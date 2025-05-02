#' Automatically process multiple statistical tests for metabolomic data
#'
#' This function performs T-test, U-test, ANOVA, Scheffe's test (ANOVA post-hoc),
#' Kruskal-Wallis test, and Dunn's test (Kruskal-Wallis post-hoc)
#' while allowing adjustment for multiple testing.
#'
#' @param Data A data.frame or data.table where the first column contains Sample IDs,
#'        the second column contains Group information, and the remaining columns
#'        are metabolites to be tested.
#' @param Adjust_p_value Logical. If TRUE (default), p-values will be adjusted for
#'        multiple testing using the method specified in Adjust_method.
#' @param Adjust_method Character string specifying the adjustment method.
#'        Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#'        "fdr", or "none". Default is "BH" (Benjamini-Hochberg).
#' @param parallel Logical. If TRUE, parallel computing will be used for calculations
#'        when the number of metabolites is ≥1000. Default is FALSE.
#'
#' @return A list containing:
#'         \item{Data}{The original input data}
#'         \item{Result}{Matrix of p-values from all tests}
#'         \item{t_test}{t-test results (if more than two groups)}
#'         \item{u_test}{U-test results (if more than two groups)}
#'         \item{Anova}{ANOVA results (if more than two groups)}
#'         \item{Anova_PostHoc}{Scheffe test results (if more than two groups)}
#'         \item{KW}{Kruskal-Wallis test results (if more than two groups)}
#'         \item{Dunn}{Dunn test results with BH adjustment (if more than two groups)}
#'
#' @import data.table
#' @importFrom matrixTests col_t_equalvar col_wilcoxon_twosample col_oneway_equalvar col_kruskalwallis
#' @importFrom PMCMRplus scheffeTest
#' @importFrom FSA dunnTest
#' @importFrom stats as.formula p.adjust
#' @importFrom utils combn
#' @importFrom parallel makeCluster detectCores clusterEvalQ parLapply stopCluster clusterExport
#'
#' @export
#'
#' @examples
#' # Load example data
#' data("Data")
#'
#' # Run statistical tests
#' library(data.table)
#' Result <- All_stats(Data)
#'
#' # View results
#' head(Result$Result)
All_stats <- function(Data,
                      Adjust_p_value = TRUE,
                      Adjust_method = "BH",
                      parallel = FALSE) {
  #### 1. Preserve original and prepare data ####
  Data_ori <- data.table::copy(Data)
  Data_renamed <- data.table::copy(Data_ori)
  data.table::setDT(Data_renamed)
  data.table::setnames(Data_renamed, 1:2, c("Sample", "Group"))
  data.table::setorder(Data_renamed, Group)
  Data_renamed[, Group := factor(Group)]
  mets <- setdiff(names(Data_renamed), c("Sample", "Group"))
  Data_renamed[, (mets) := lapply(.SD, as.numeric), .SDcols = mets]
  data.table::setkey(Data_renamed, Group)

  #### 2. Precompute group matrices ####
  groups <- levels(Data_renamed$Group)
  pairs <- utils::combn(groups, 2, simplify = FALSE)
  pair_names <- vapply(pairs, paste, character(1), collapse = "-")
  n_pairs <- length(pairs)
  mats <- setNames(lapply(groups, function(g) {
    as.matrix(Data_renamed[.(g), mets, with = FALSE])
  }), groups)

  #### 3. Parallel setup ####
  if (parallel) {
    cores <- max(1, min(5, parallel::detectCores() - 1))
    message("Parallel setup: using ", cores, " cores.")
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(matrixTests))
    parallel::clusterEvalQ(cl, library(PMCMRplus))
    parallel::clusterEvalQ(cl, library(FSA))
    parallel::clusterExport(cl, c("mats", "mets", "pairs"), envir = environment())
    apply_func <- function(X, FUN, ...) {
      parallel::parLapplyLB(cl, X, FUN, ...)
    }
  } else {
    apply_func <- function(X, FUN, ...) {
      lapply(X, FUN, ...)
    }
  }

  #### 4. t-test ####
  t_list <- apply_func(seq_len(n_pairs), function(i) {
    p <- pairs[[i]]
    matrixTests::col_t_equalvar(mats[[p[1]]], mats[[p[2]]])$pvalue
  })
  t_mat <- do.call(cbind, t_list)
  colnames(t_mat) <- paste(pair_names, "t-test", sep = "___")
  rownames(t_mat) <- mets
  message("t-test completed.")

  #### 5. u-test ####
  thresh <- 10
  max_grp_n <- Data_renamed[, .N, by = Group][which.max(N), N]
  has_ties <- any(vapply(
    mats, function(m) any(apply(m, 2, anyDuplicated) > 0),
    logical(1)
  ))
  do_exact <- max_grp_n <= thresh && !has_ties
  message(sprintf(
    "Exact test set to %s (max N = %d, ties = %s)",
    do_exact, max_grp_n, has_ties
  ))

  u_list <- apply_func(seq_len(n_pairs), function(i) {
    p <- pairs[[i]]
    matrixTests::col_wilcoxon_twosample(mats[[p[1]]], mats[[p[2]]], exact = do_exact)$pvalue
  })
  u_mat <- do.call(cbind, u_list)
  colnames(u_mat) <- paste(pair_names, "u-test", sep = "___")
  rownames(u_mat) <- mets
  message("u-test completed.")

  #### 6–9. Only if more than two groups ####
  if (length(groups) > 2) {
    #### 6. ANOVA ####
    df_anova <- data.frame(
      p_anova = matrixTests::col_oneway_equalvar(as.matrix(Data_renamed[, mets, with = FALSE]), Data_renamed$Group)$pvalue,
      row.names = mets
    )
    message("ANOVA completed.")

    #### 7. Kruskal-Wallis ####
    df_kw <- data.frame(
      p_kw = matrixTests::col_kruskalwallis(as.matrix(Data_renamed[, mets, with = FALSE]), Data_renamed$Group)$pvalue,
      row.names = mets
    )
    message("Kruskal-Wallis completed.")

    #### 8. Scheffe post-hoc ####
    get_pvec <- function(pmat) {
      vapply(pairs, function(pr) {
        pmat[pr[2], pr[1]]
      }, numeric(1))
    }
    Data_V <- as.data.frame(Data_renamed)
    safe <- paste0("V", seq_along(mets))
    colnames(Data_V)[-(1:2)] <- safe
    mmap <- setNames(safe, mets)
    sch_list <- apply_func(seq_along(mets), function(i) {
      pm <- PMCMRplus::scheffeTest(as.formula(paste(mmap[[mets[i]]], "~ Group")), data = Data_V[, c("Group", mmap[[mets[i]]])])$p.value
      get_pvec(pm)
    })
    df_scheffe <- do.call(rbind, sch_list)
    colnames(df_scheffe) <- paste(pair_names, "SCH_posthoc", sep = "___")
    rownames(df_scheffe) <- mets
    message("Scheffe post-hoc completed.")

    #### 9. Dunn post-hoc ####
    lex_pair_names <- sort(pair_names)
    canon <- function(x) {
      vapply(strsplit(gsub("\\s+", "", x), "-"), function(v) {
        paste(sort(v), collapse = "-")
      }, character(1))
    }
    dun_list <- apply_func(seq_along(mets), function(i) {
      dt <- FSA::dunnTest(as.formula(paste(mmap[[mets[i]]], "~ Group")),
        data = Data_V[, c("Group", mmap[[mets[i]]])],
        method = "none"
      )$res
      pvec <- setNames(dt$P.unadj, canon(dt$Comparison))
      pvec <- pvec[lex_pair_names]
      stats::p.adjust(pvec, method = Adjust_method)
    })
    df_dunn <- do.call(rbind, dun_list)
    colnames(df_dunn) <- paste(lex_pair_names, "DUNN_posthoc", sep = "___")
    rownames(df_dunn) <- mets
    message("Dunn post-hoc completed.")
  }

  #### 10. P-value adjustment ####
  if (Adjust_p_value) {
    df_t <- as.data.frame(apply(t_mat, 2, stats::p.adjust, method = Adjust_method))
    df_u <- as.data.frame(apply(u_mat, 2, stats::p.adjust, method = Adjust_method))
    rownames(df_t) <- rownames(df_u) <- mets
    if (length(groups) > 2) {
      df_a <- data.frame(
        p_anova = stats::p.adjust(df_anova$p_anova, method = Adjust_method),
        row.names = mets
      )
      df_k <- data.frame(
        p_kw = stats::p.adjust(df_kw$p_kw, method = Adjust_method),
        row.names = mets
      )
    }
  } else {
    df_t <- as.data.frame(t_mat)
    df_u <- as.data.frame(u_mat)
    if (length(groups) > 2) {
      df_a <- df_anova
      df_k <- df_kw
    }
  }
  message("P-value adjustment completed.")

  #### 11. Compile results ####
  Result <- cbind(df_t, df_u)
  if (length(groups) > 2) {
    Result <- cbind(Result, df_a, df_scheffe, df_k, df_dunn)
  }
  Final <- list(
    Data         = Data_ori,
    Data_renamed = as.data.frame(Data_renamed),
    Result       = as.matrix(Result),
    t_test       = df_t,
    u_test       = df_u
  )
  if (length(groups) > 2) {
    Final$Anova <- df_a
    Final$Scheffe <- df_scheffe
    Final$KW <- df_k
    Final$Dunn <- df_dunn
  }
  message("All statistical analyses completed.")
  return(Final)
}
