#' Create and save barplots for all metabolites with significance indicators
#'
#' This function automatically generates and saves barplots for all metabolites in the dataset,
#' with asterisks indicating statistical significance between groups based on the selected
#' statistical test. Plots are saved as PNG files in a 'barplot' directory.
#'
#' @param data List inheriting from Allstats containing metabolite data and statistical test results
#' @param asterisk Statistical test to use for significance indicators: "Dunn", "Scheffe", "u_test", or "t_test"
#' @param significant_variable_only Logical; if TRUE, only plots significant variables
#' @param colour Vector of color codes used for different groups in the barplots
#' @param legend_position Position of the legend: "none", "left", "right", "bottom", or "top"
#' @param order Optional vector specifying the order of groups (not currently used in function)
#' @param width Width of the barplots
#' @param size Line size for the bar outlines and error bars
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
#' @importFrom stats setNames sd
#' @return Invisibly returns NULL, with plots saved to the 'barplot' directory
#' @export
#'
#' @examples
#' \dontrun{
#' data(Data)
#' Test <- All_stats(Data)
#' Barplot_new(Test,
#'   asterisk = "t_test",
#'   significant_variable_only = FALSE,
#'   colour = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", "#660099")
#' )
#' }
Barplot_new <- function(data,
                        asterisk = "t_test",
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
                        width = 0.6,
                        size = 0.5,
                        label_size = 2.8,
                        tip = 0.01,
                        step = 0.05,
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

  ## 0) 필수 패키지 ----------------------------------------------------------
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

  ## 1) p-value 행렬 ---------------------------------------------------------
  pick <- switch(tolower(asterisk),
    dunn = "Dunn",
    scheffe = "Scheffe",
    u_test = "u_test",
    t_test = "t_test",
    stop("asterisk must be Dunn, Scheffe, u_test, or t_test")
  )
  p <- as.matrix(if (is.list(data[[pick]])) {
    unlist(data[[pick]])
  } else {
    data[[pick]]
  })
  p[is.nan(p)] <- 1

  ## 2) 데이터 & 안전 변수명 --------------------------------------------------
  raw <- colnames(data$Data)[-c(1, 2)]
  safe <- make.names(raw, unique = TRUE)
  map <- setNames(raw, safe)

  dat <- data$Data_renamed %>%
    dplyr::transmute(dplyr::across(-c(1, 2), as.numeric), Group = .[[2]])
  colnames(dat)[-ncol(dat)] <- safe
  if (!is.null(order)) {
    dat$Group <- factor(dat$Group, levels = order)
  }

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

  ## 3) stat-table -----------------------------------------------------------
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

  ## 4) 공통 테마 ------------------------------------------------------------
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

  dir.create("barplot", showWarnings = FALSE)

  ## 5) draw 함수 ------------------------------------------------------------
  draw <- function(i) {
    st <- tbl[[i]]
    show <- !is.null(st) && is.data.frame(st) && nrow(st) > 0
    if (!show && significant_variable_only) {
      return(invisible(NULL))
    }

    sumdat <- dat %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(
        mean = mean(.data[[safe[i]]], na.rm = TRUE),
        se = sd(.data[[safe[i]]], na.rm = TRUE) /
          sqrt(sum(!is.na(.data[[safe[i]]]))),
        .groups = "drop"
      ) %>%
      tidyr::complete(Group = levels(dat$Group)) %>%
      dplyr::mutate(
        Group = factor(Group, levels = levels(dat$Group)),
        se = ifelse(is.na(se) | !is.finite(se), 0, se)
      )

    p <- ggplot2::ggplot(sumdat, ggplot2::aes(x = Group, y = mean, fill = Group)) +
      ggplot2::geom_col(
        width = width,
        colour = "black",
        size = size
      ) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se),
        width = width / 2,
        size = size
      ) +
      ggplot2::scale_fill_manual(
        values = pal,
        limits = levels(dat$Group),
        drop = FALSE
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::scientific_format(),
        expand = ggplot2::expansion(mult = c(0, 0.1))
      ) +
      ggplot2::labs(title = map[[safe[i]]], y = "Intensity") +
      base

    if (show) {
      y_pos_sig <- 1.05 * max(sumdat$mean + sumdat$se, na.rm = TRUE)
      if (!is.finite(y_pos_sig) ||
        y_pos_sig <= 0) {
        y_pos_sig <- max(sumdat$mean, na.rm = TRUE) * 1.1
      }

      st_filtered <- st %>%
        dplyr::filter(group1 %in% levels(dat$Group) &
          group2 %in% levels(dat$Group))

      if (nrow(st_filtered) > 0) {
        p <- p +
          ggpubr::stat_pvalue_manual(
            data = st_filtered,
            y.position = y_pos_sig,
            label = "p.adj.signif",
            tip.length = tip,
            step.increase = step,
            label.size = label_size,
            # dotplot 코드 참고: size = 3.5 추가, bracket.size 제거
            size = 3.5,
            # dotplot 코드 참고: vjust = 0.05 사용
            vjust = .05,
            # dotplot 코드 참고: hide.ns 제거 (st_filtered에서 이미 필터링 가정)
            # dotplot 코드 참고: mapping 인수 제거
            inherit.aes = FALSE # 명시적으로 FALSE 유지
          )
      }
    }

    safe_filename <- sanitize(map[[safe[i]]])
    save_path <- file.path("barplot")
    save_filename <- file.path(save_path, sprintf("%s_bar.png", safe_filename))

    suppressMessages(
      ggplot2::ggsave(
        filename = save_filename,
        plot     = p,
        device   = ragg::agg_png,
        width    = fig_width,
        height   = fig_height,
        units    = "in",
        dpi      = 600
      )
    )

    return(invisible(NULL))
  }

  ## 6) 실행 (병렬 ≤ 5 코어) -------------------------------------------------
  N <- length(safe)
  cores <- min(5, parallel::detectCores())

  vars_to_export <- c("asterisk")

  if (requireNamespace("doSNOW", quietly = TRUE)) {
    cl <- parallel::makeCluster(min(cores, N))
    doSNOW::registerDoSNOW(cl)

    # --- 수정된 on.exit 설정 ---
    # on.exit 호출을 한 번만 사용하여 모든 정리 작업 수행
    pb <- NULL # pb 초기화
    if (progress_bar) {
      pb <- utils::txtProgressBar(max = N, style = 3)
      prg <- function(n) {
        utils::setTxtProgressBar(pb, n)
      }
      opts <- list(progress = prg)
      # on.exit: 프로그레스 바 닫고 클러스터 중지
      on.exit(
        {
          if (!is.null(pb)) {
            close(pb)
          }
          parallel::stopCluster(cl)
        },
        add = TRUE
      )
    } else {
      opts <- list()
      # on.exit: 프로그레스 바 없을 때는 클러스터만 중지
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }
    # --- 수정 끝 ---

    foreach::foreach(
      i = seq_len(N),
      .options.snow = opts,
      .packages = c(
        "ggplot2",
        "ggpubr",
        "dplyr",
        "tidyr",
        "scales",
        "ragg",
        "stats",
        "magrittr"
      ),
      .export = vars_to_export # 경고가 나오더라도 일단 유지
    ) %dopar% {
      tryCatch(
        {
          draw(i)
        },
        error = function(e) {
          warning("Parallel task failed for index ",
            i,
            " (",
            map[[safe[i]]],
            "): ",
            e$message,
            call. = FALSE
          )
          NULL
        }
      )
    }

    # 루프 정상 종료 시 프로그레스 바 닫기 (on.exit 에서도 처리되지만 명시적 추가)
    # if (progress_bar && !is.null(pb)) close(pb) # on.exit에서 처리하므로 중복 불필요
  } else {
    pb <- NULL
    if (progress_bar) {
      pb <- utils::txtProgressBar(
        min = 0,
        max = N,
        style = 3
      )
      on.exit(
        {
          if (!is.null(pb)) {
            close(pb)
          }
        },
        add = TRUE
      )
    }
    for (i in seq_len(N)) {
      tryCatch(
        {
          draw(i)
        },
        error = function(e) {
          warning("Sequential task failed for index ",
            i,
            " (",
            map[[safe[i]]],
            "): ",
            e$message,
            call. = FALSE
          )
        }
      )
      if (progress_bar &&
        !is.null(pb)) {
        utils::setTxtProgressBar(pb, i)
      }
    }
    if (progress_bar && !is.null(pb)) {
      close(pb)
    }
  }
  invisible(NULL)
}
