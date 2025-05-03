#' Create and save dotplots for all metabolites with significance indicators
#'
#' This function automatically generates and saves dotplots for all metabolites in the dataset,
#' with asterisks indicating statistical significance between groups based on the selected
#' statistical test. Plots are saved as PNG files in a 'dotplot' directory.
#'
#' @param data List inheriting from Allstats containing metabolite data and statistical test results
#' @param asterisk Statistical test to use for significance indicators: "Dunn", "Scheffe", "u_test", or "t_test"
#' @param significant_variable_only Logical; if TRUE, only plots significant variables
#' @param colour Vector of color codes used for different groups in the dotplots
#' @param legend_position Position of the legend: "none", "left", "right", "bottom", or "top"
#' @param order Optional vector specifying the order of groups (not currently used in function)
#' @param width Width of the dotplots
#' @param size Line size for the dotplot outlines
#' @param label_size Size of significance labels (asterisks)
#' @param tip Length of the significance bracket tips
#' @param step Vertical step increase between significance brackets
#' @param X_text Size of X-axis text
#' @param Y_text Size of Y-axis title
#' @param Y_lab Size of Y-axis labels
#' @param T_size Size of plot title
#' @param fig_width Width of saved figure in inches
#' @param fig_height Height of saved figure in inches
#' @param sig Vector of two significance thresholds: first for displaying any significance,
#'        second for displaying higher significance (** vs *)
#' @param progress_bar Logical; if TRUE, displays a progress bar during plot generation
#'
#' @importFrom foreach %dopar%
#' @importFrom stats setNames
#' @return Invisibly returns NULL, with plots saved to the 'dotplot' directory
#' @export
#'
#' @examples
#' \dontrun{
#' data(Data)
#' Test <- All_stats(Data)
#' Dotplot_new(Test,
#'   asterisk = "Dunn",
#'   significant_variable_only = FALSE,
#'   colour = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", "#660099")
#' )
#' }
Dotplot_new <- function(data,
                        asterisk = "t_test",
                        # "Dunn","Scheffe","u_test","t_test"
                        significant_variable_only = FALSE,
                        colour = c(
                          "#FF3300",
                          "#FF6600",
                          "#FFCC00",
                          "#99CC00",
                          "#0066CC",
                          "#660099"
                        ),
                        legend_position = "none",
                        order = NULL,
                        width = .3,
                        size = .5,
                        label_size = 2.8,
                        tip = .01,
                        step = .05,
                        X_text = 10,
                        Y_text = 12,
                        Y_lab = 10,
                        T_size = 15,
                        fig_width = 6,
                        fig_height = 6,
                        sig = c(.05, .01),
                        progress_bar = TRUE) {
  `%>%` <- magrittr::`%>%`
  `%dopar%` <- foreach::`%dopar%`
  sanitize <- \(x) gsub("[\\\\/:*?\"<>|]", "_", x)

  # 0) 필수 패키지 ----------------------------------------------------------
  need <- c(
    "ggplot2",
    "ggpubr",
    "dplyr",
    "tidyr",
    "gtools",
    "foreach",
    "scales",
    "ragg"
  )
  miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(miss)) {
    stop("Install packages: ", paste(miss, collapse = ", "))
  }

  # 1) p-value 행렬 ---------------------------------------------------------
  pick <- switch(tolower(asterisk),
    dunn    = "Dunn",
    scheffe = "Scheffe",
    u_test  = "u_test",
    t_test  = "t_test",
    stop("asterisk must be Dunn, Scheffe, u_test, or t_test")
  )
  p <- data[[pick]]
  p <- as.matrix(if (is.list(p)) {
    unlist(p)
  } else {
    p
  })
  p[is.nan(p)] <- 1

  # 2) 데이터 + 이름 안전화 --------------------------------------------------
  raw <- colnames(data$Data)[-c(1, 2)]
  safe <- make.names(raw, unique = TRUE)
  map <- setNames(raw, safe)

  dat <- data$Data_renamed %>%
    dplyr::transmute(dplyr::across(-c(1, 2), as.numeric), Group = .[[2]])
  colnames(dat)[-ncol(dat)] <- safe

  if (!is.null(order)) {
    # 유효성 검사
    bad <- setdiff(order, unique(dat$Group))
    if (length(bad))
      stop("Invalid group names in 'order': ",
           paste(bad, collapse = ", "))
    dat$Group <- factor(dat$Group, levels = order)
  } else {
    dat$Group <- factor(dat$Group)
  }
  
  grp <- levels(dat$Group)
  pal <- colour[seq_along(grp)]
  Comb <- gtools::combinations(length(grp), 2, grp)

  # 3) stat-table 미리 계산 --------------------------------------------------
  tbl <- lapply(seq_along(safe), function(i) {
    if (length(grp) == 2) {
      tibble::tibble(
        group1 = Comb[, 1],
        group2 = Comb[, 2],
        p = as.numeric(p[i, ])
      ) %>%
        dplyr::filter(p <= sig[1]) %>%
        dplyr::mutate(p.adj.signif = ifelse(p <= sig[2], "**", "*"))
    } else {
      tibble::tibble(
        pair = sub("___.*", "", colnames(p)),
        p = as.numeric(p[i, ])
      ) %>%
        tidyr::separate(pair,
          c("group1", "group2"),
          sep = " - |-",
          remove = FALSE
        ) %>%
        dplyr::filter(p <= sig[1]) %>%
        dplyr::mutate(p.adj.signif = ifelse(p <= sig[2], "**", "*")) %>%
        dplyr::select(group1, group2, p, p.adj.signif)
    }
  })

  # 4) 공통 테마 ------------------------------------------------------------
  base <- ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = legend_position,
      plot.title = ggplot2::element_text(
        size = T_size,
        face = "bold",
        hjust = .5
      ),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = Y_text),
      axis.text.x = ggplot2::element_text(size = X_text, vjust = .5),
      axis.text.y = ggplot2::element_text(size = Y_lab)
    )

  dir.create("dotplot", showWarnings = FALSE)

  # 5) 플롯 + 저장 함수 -----------------------------------------------------
  draw <- function(i) {
    st <- tbl[[i]]
    show <- nrow(st) > 0
    if (!show && significant_variable_only) {
      return()
    }

    p <- ggplot2::ggplot(dat, ggplot2::aes(Group, .data[[safe[i]]], colour = Group)) +
      ggforce::geom_sina(maxwidth = width, size = 1) +
      ggplot2::stat_summary(
        fun.y = median,
        fun.ymin = median,
        fun.ymax = median,
        geom = "crossbar",
        width = 0.4,
        color = "Black",
        size = 0.2
      ) +
      ggplot2::scale_colour_manual(values = pal) +
      ggplot2::scale_y_continuous(labels = scales::scientific_format()) +
      ggplot2::labs(title = map[[safe[i]]], y = "Intensity") +
      base

    if (show) {
      p <- p + ggpubr::stat_pvalue_manual(
        st,
        y.position = 1.05 * max(dat[[safe[i]]], na.rm = TRUE),
        step.increase = step,
        tip.length = tip,
        label.size = label_size,
        label = "p.adj.signif",
        size = 3.5,
        vjust = .05
      )
    }

    suppressMessages(
      ggplot2::ggsave(
        file.path("dotplot", sprintf("%s_dot.png", sanitize(map[[safe[i]]]))),
        plot = p,
        device = ragg::agg_png,
        width = fig_width,
        height = fig_height,
        dpi = 600
      )
    )
  }

  # 6) 실행 (병렬 ≤ 5 코어) -------------------------------------------------
  N <- length(safe)
  cores <- min(5, parallel::detectCores())

  if (requireNamespace("doSNOW", quietly = TRUE)) {
    cl <- parallel::makeCluster(min(cores, N))
    doSNOW::registerDoSNOW(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    if (progress_bar) {
      pb <- utils::txtProgressBar(max = N, style = 3)
      progress <- function(n) {
        utils::setTxtProgressBar(pb, n)
      }
      opts <- list(progress = progress)
      on.exit(close(pb), add = TRUE)
    } else {
      opts <- list()
    }

    foreach::foreach(
      i = seq_len(N),
      .options.snow = opts,
      .packages = c("ggplot2", "ggpubr", "dplyr", "tidyr", "scales", "ragg")
    ) %dopar% draw(i)
  } else {
    # 순차
    if (progress_bar) {
      pb <- utils::txtProgressBar(
        min = 0,
        max = N,
        style = 3
      )
      on.exit(close(pb), add = TRUE)
    }
    for (i in seq_len(N)) {
      draw(i)
      if (progress_bar) {
        utils::setTxtProgressBar(pb, i)
      }
    }
  }
  invisible(NULL)
}
