#' Automatically processes T-test, U-test, Anova, Scheffe(Anova Post-Hoc), Krukal Wallis, Dunn-test(BH adjusted,(Kurkal Wallis Post-Hoc)) while allowing adjustment of FDR
#'
#' @param Data csv file with Header as False First column with Sample Second column with Multilevel(Mixomics) so that it can be compatible with other multivariate statistics Third column with Group information. Rest of the columns are the metabolites to be tested.
#' @param Adjust_p_value Set True if FDR adjustments are to be made. If not set False
#' @param Adjust_method adjustment methods frequently used. "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#'
#' @return List including Result Matrix of p-values, converted datas.
#' @export
#'
#' @examples
#' data("Data")
#' Result <- All_stats(Data)
All_stats <- function(Data, Adjust_p_value = T, Adjust_method = "BH") {
  LETTERS210729 <- paste0("V", 1:500000)
  colnames(Data) <- Data[1, ]
  Data <- Data[-1, -2]
  Data <- Data %>% dplyr::arrange(Data$Group)
  Data_renamed <- Data
  nmet <- ncol(Data) - 2
  colnames(Data_renamed) <- c(colnames(Data[1:2]), LETTERS210729[1:nmet])
  rownames(Data_renamed) <- Data[, 1]
  Data_renamed_raw <- Data_renamed[, -c(1, 2)]
  Data_renamed_raw <- apply(Data_renamed_raw, 2, as.numeric)
  Data_final <- cbind(Data[, 1:2], Data_renamed_raw)
  Data_tmp <- data.table::as.data.table(Data_final)
  Data_tmp1 <- Data_tmp
  Data_tmp1$Group <- as.factor(Data_tmp1$Group)
  
  
  groups_split <- split(Data_tmp, Data_tmp$Group)
  #### ttest ####
  split_t_test <- function(x, i) {
    if (var(groups_split[[x[1]]][[i]]) == 0 & var(groups_split[[x[2]]][[i]]) == 0) {
      1
    } else {
      t.test(
        groups_split[[x[1]]][[i]],
        groups_split[[x[2]]][[i]]
      )[["p.value"]]
    }
  }
  
  
  res_ttest <- lapply(
    (seq_len(ncol(Data_tmp) - 2) + 2),
    function(i) as.list(combn(names(groups_split), 2, split_t_test, i = i))
  )
  
  
  df_ttest <- data.table::rbindlist(res_ttest)
  df_ttest <- data.frame(df_ttest)
  
  print("T-test has finished")
  
  #### utest ####
  split_u_test <- function(x, i) {
    if (var(groups_split[[x[1]]][[i]]) == 0 & var(groups_split[[x[2]]][[i]]) == 0) {
      if (groups_split[[x[1]]][[i]] == groups_split[[x[2]]][[i]]) {
        1
      } else {
        wilcox.test(
          groups_split[[x[1]]][[i]],
          groups_split[[x[2]]][[i]]
        )[["p.value"]]
      }
    } else {
      wilcox.test(
        groups_split[[x[1]]][[i]],
        groups_split[[x[2]]][[i]]
      )[["p.value"]]
    }
  }
  res_utest <- lapply(
    (seq_len(ncol(Data_tmp) - 2) + 2),
    function(i) as.list(combn(names(groups_split), 2, split_u_test, i = i))
  )
  df_utest <- data.table::rbindlist(res_utest)
  df_utest <- data.frame(df_utest)
  print("U-test has finished")
  
  
  
  if (length(unique(Data_final$Group)) > 2) {
    groups_split <- split(Data_tmp, Data_tmp$Group)
    
    
    #### ANOVAPostHoc ####
    formula_anova <- lapply(colnames(Data_tmp)[3:ncol(Data_tmp)], function(x) as.formula(paste0(x, " ~ Group")))
    res_anova <- lapply(formula_anova, function(x) summary(aov(x, data = Data_tmp)))
    names(res_anova) <- format(formula_anova)
    p_anova <- unlist(lapply(res_anova, function(x) x[[1]]$"Pr(>F)"[1]))
    
    df_anova <- data.table::data.table(p_anova)
    df_anova <- data.frame(df_anova)
    
    anovapost_name <- lapply(colnames(Data_tmp)[3:ncol(Data_tmp)], function(x) as.formula(paste0(x, " ~ Group")))
    res_anovapost <- lapply(anovapost_name, function(x) DescTools::PostHocTest(aov(x, data = Data_tmp), method = "scheffe"))
    names(res_anovapost) <- format(anovapost_name)
    post_anova <- lapply(res_anovapost, function(x) x[["Group"]][, 4])
    
    df_anova_post <- t(data.frame(post_anova))
    
    
    print("Anova & PostHoc has finished")
    
    
    #### KruskalWallisPostHoc ####
    formula_kw <- lapply(colnames(Data_tmp)[3:ncol(Data_tmp)], function(x) as.formula(paste0(x, " ~ Group")))
    res_kw <- lapply(formula_anova, function(x) kruskal.test(x, data = Data_tmp)[["p.value"]])
    names(res_kw) <- format(formula_kw)
    p_kw <- unlist(res_kw)
    
    df_kw <- data.table::data.table(p_kw)
    df_kw <- data.frame(df_kw)
    
    kwpost_name <- lapply(colnames(Data_tmp)[3:ncol(Data_tmp)], function(x) as.formula(paste0(x, " ~ Group")))
    res_kwpost <- lapply(kwpost_name, function(x) FSA::dunnTest(x, data = Data_tmp1, method = "none"))
    names(res_kwpost) <- format(kwpost_name)
    post_kw <- lapply(res_kwpost, function(x) x[["res"]][["P.unadj"]])
    post_kw <- lapply(post_kw, function(x) p.adjust(x, method = "BH"))
    
    df_kw_post <- t(data.table::data.table(data.frame(post_kw)))
    df_kw_post <- data.frame(df_kw_post)
    
    print("Kruskal Wallis & PostHoc has finished")
    
    #### finalization ####
    
    Names <- NULL
    for (i in 1:choose(length(groups_split), 2)) {
      Names <- rbind(Names, paste(combn(names(groups_split), 2)[1, i],
                                  combn(names(groups_split), 2)[2, i],
                                  sep = "-"
      ))
    }
    
    AN_Post_names <- colnames(df_anova_post)
    DU_post_names <- res_kwpost[[1]]$res$Comparison
    
    rownamechange <- colnames(Data)[3:ncol(Data_final)]
    rownames(df_ttest) <- rownamechange
    rownames(df_utest) <- rownamechange
    rownames(df_anova) <- rownamechange
    rownames(df_anova_post) <- rownamechange
    rownames(df_kw) <- rownamechange
    rownames(df_kw_post) <- rownamechange
    
    colnames(df_ttest) <- paste(Names[,1], "t-test", sep = "___")
    colnames(df_utest) <- paste(Names[,1], "u-test", sep = "___")
    colnames(df_anova) <- "Anova"
    colnames(df_anova_post) <- paste(AN_Post_names, "ANO_posthoc", sep = "___")
    colnames(df_kw) <- "Kruskal_Wallis"
    colnames(df_kw_post) <- paste(DU_post_names, "Kru_posthoc(Dunn)", sep = "___")
    
    if (Adjust_p_value == T) {
      print("###########################################")
      print(paste0(
        "adjusted according to the ",
        Adjust_method, " method"
      ))
      print("###########################################")

        adj_func <- function(x) {
        p.adjust(x, method = Adjust_method)
      }
      
      df_ttest <- apply(df_ttest, 2, adj_func)
      df_utest <- apply(df_utest, 2, adj_func)
      df_anova <- apply(df_anova, 2, adj_func)
      df_kw <- apply(df_kw, 2, adj_func)
    } else {
      print("###########################################")
      print("p_value not adjusted")
      print("###########################################")
    }
    
    Result <- cbind(df_ttest, df_utest, df_anova, df_anova_post, df_kw, df_kw_post)
    Final <- list()
    Final$Data <- Data
    Final$Data_renamed <- Data_renamed
    Final$Result <- Result
    Final$Anova <- df_anova
    Final$Anova_PostHoc <- df_anova_post
    Final$KW <- df_kw
    Final$Dunn <- df_kw_post
    Final$t_test <- df_ttest
    Final$u_test <- df_utest
    Final
  }
}
