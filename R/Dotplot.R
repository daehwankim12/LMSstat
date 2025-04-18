#' Automatically save dotplots for all metabolites including asterisks for significance.
#'
#' @param data List inheriting from Allstats
#' @param asterisk Choose asterisk to plot ("Dunn","Scheffe","u_test","t_test")
#' @param significant_variable_only show significant variable only (T,F)
#' @param color colors used for ggplots.
#' @param order order of the groups; order = c("A","B","C")
#' @param legend_position legend position "none","left","right","bottom","top"
#' @param tip_length significance tip length
#' @param label_size significance label size
#' @param step_increase significance step increase
#' @param size dot size
#' @param fig_width figure size
#' @param fig_height figure size
#' @param Y_text Y axis title size
#' @param X_text X axis text size
#' @param Y_lab y axis text size
#' @param T_size Title size
#' @param sig_int significance parameter
#'
#' @importFrom foreach %dopar%
#' @return ggdotplot
#' @export
#'
#' @examples data(Data)
#' Test <- All_stats(Data)
#' Dotplot(Test,
#'   asterisk = "Dunn", significant_variable_only = F,
#'   color = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", "#660099")
#' )
Dotplot <- function(data,
                    asterisk = "t_test",
                    significant_variable_only = F,
                    color = c(
                      "#FF3300",
                      "#FF6600",
                      "#FFCC00",
                      "#99CC00",
                      "#0066CC",
                      "#660099"
                    ),
                    legend_position = "none",
                    order = NULL,
                    tip_length = 0.01,
                    label_size = 2.88,
                    step_increase = 0.05,
                    size = NULL,
                    fig_width = NA,
                    fig_height = NA,
                    X_text = 10,
                    Y_text = 12,
                    Y_lab = 10,
                    T_size = 15,
                    sig_int = c(0.05, 0.01)) {
  {
    # Summary
    ### Plot_data_prep###
    ifelse(!dir.exists(file.path(getwd(), "dotplot")), dir.create(file.path(getwd(), "dotplot")), FALSE)
    data[["Data_renamed"]] <-
      data[["Data_renamed"]] %>% mutate(ZZZZ = data[["Data_renamed"]][, 2])
    data[["Data_renamed"]] <- data[["Data_renamed"]][, c(-1, -2)]
    data[["Data_renamed"]][, 1:(ncol(data[["Data_renamed"]]) - 1)] <-
      sapply(
        data[["Data_renamed"]][, 1:(ncol(data[["Data_renamed"]]) - 1)],
        function(x) {
          as.numeric(x)
        }
      )
    NAMES <- colnames(data[["Data"]][3:ncol(data[["Data"]])])
    data[["Data_renamed"]]$ZZZZ <-
      as.factor(data[["Data_renamed"]]$ZZZZ)
    colnames(data[["Data_renamed"]])[ncol(data[["Data_renamed"]])] <-
      "Group"
    Comb <-
      gtools::combinations(length(unique(data[["Data"]]$Group)), 2, unique(data[["Data"]]$Group))
    ckey <- color[1:length(unique(data[["Data"]]$Group))]

    if (asterisk == "Scheffe") {
      p_val_data <- data[["Anova_PostHoc"]]
    } else if (asterisk == "t_test") {
      p_val_data <- data[["t_test"]]
    } else if (asterisk == "u_test") {
      p_val_data <- data[["u_test"]]
    } else if (asterisk == "Dunn") {
      p_val_data <- data[["Dunn"]]
    } else {
      print("Wrong asterisk input must be one of (Dunn,Scheffe,u_test,t_test)")
    }
    if (length(unique(data[["Data"]]$Group)) != 2) {
      for (a in 1:nrow(p_val_data)) {
        for (b in 1:ncol(p_val_data)) {
          if (is.nan(p_val_data[a, b]) == T) {
            p_val_data[a, b] <- 1
          }
        }
      }
    } else if (length(unique(data[["Data"]]$Group)) == 2) {
      for (b in 1:length(p_val_data)) {
        if (is.nan(p_val_data[b]) == T) {
          p_val_data[b] <- 1
        }
      }
    }
  }
  ### Plots###
  group_not_two <-
    length(unique(data[["Data_renamed"]][["Group"]])) != 2
  tryCatch(expr = {
    num_core <- parallel::detectCores()
    final_cores <- min(5, num_core)
    cl <- parallel::makeCluster(final_cores)
    doSNOW::registerDoSNOW(cl)
    if (group_not_two) {
      pb <- txtProgressBar(max = nrow(p_val_data), style = 3)
    } else {
      pb <- txtProgressBar(max = length(p_val_data), style = 3)
    }
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)
    suppressWarnings(if (group_not_two) {
      foreach::foreach(
        number = 1:nrow(p_val_data),
        .options.snow = opts,
        .packages = c("LMSstat", "dplyr", "plyr")
      ) %dopar% {
        stat.test <- as.data.frame(matrix(
          data = NA,
          nrow = nrow(Comb),
          ncol = 2
        ))
        stat.test[, 1:2] <- Comb
        if (asterisk == "Dunn") {
          rownames(stat.test) <- paste0(stat.test[, 1], " - ", stat.test[, 2])
        } else if (asterisk == "Scheffe") {
          rownames(stat.test) <- paste0(stat.test[, 2], "-", stat.test[, 1])
        } else {
          rownames(stat.test) <- paste0(stat.test[, 1], "-", stat.test[, 2])
        }
        colnames(p_val_data) <-
          as.data.frame(strsplit(colnames(p_val_data), "___"))[1, ]
        t1 <-
          as.data.frame(t(as.data.frame(p_val_data)[number, ]))
        t2 <- stat.test
        t3 <- merge(t2, t1, by = 0)
        stat.test <- t3[, 2:4]
        colnames(stat.test) <- c(
          "group1", "group2",
          "p"
        )
        stat.test <- stat.test %>% plyr::mutate(
          p.adj.signif = case_when(
            p >
              sig_int[1] ~ "NS",
            p <= sig_int[1] & p > sig_int[2] ~ "*",
            p <= sig_int[2] ~ "**"
          )
        )
        stat.test <- stat.test[stat.test$p.adj.signif != "NS", ]
        if (length(stat.test > 4)) {
          ggpubr::ggstripchart(
            data[["Data_renamed"]],
            x = "Group",
            y = colnames(data[["Data_renamed"]])[number],
            color = "Group",
            palette = ckey,
            order = order,
            size = size
          ) +
            ggplot2::stat_summary(
              fun.y = median,
              fun.ymin = median,
              fun.ymax = median,
              geom = "crossbar",
              width = 0.4,
              color = "Black",
              size = 0.2
            ) +
            ggplot2::scale_y_continuous(label = scientific_10x) +
            ggplot2::labs(
              title = NAMES[number],
              x = NULL,
              y = "Intensity"
            ) +
            ggplot2::theme_classic() + ggplot2::theme(
              plot.title = ggplot2::element_text(
                size = T_size,
                face = "bold",
                color = "Black",
                hjust = 0.5,
                lineheight = 1.2
              ),
              # title
              plot.subtitle = ggplot2::element_text(
                size = 15,
                hjust = 0.5
              ),
              # subtitle
              plot.caption = ggplot2::element_text(size = 15),
              # caption
              axis.title.x = ggplot2::element_text(
                vjust = 10,
                size = 15
              ),
              # X axis title
              axis.title.y = ggplot2::element_text(size = Y_text),
              # Y axis title
              axis.text.x = ggplot2::element_text(
                size = X_text,
                vjust = .5
              ),
              # X axis text
              axis.text.y = ggplot2::element_text(size = Y_lab),
              legend.position = legend_position
            ) +
            ggpubr::stat_pvalue_manual(
              stat.test,
              y.position = 1.05 * max(data[["Data_renamed"]][, number]),
              step.increase = step_increase,
              label.size = label_size,
              tip.length = tip_length,
              label = "p.adj.signif",
              size = 3.5,
              vjust = 0.05
            )
          ggplot2::ggsave(
            filename = paste(NAMES[number], "dotplot.png", collapse = ""),
            path = paste0(getwd(), "/dotplot"),
            width = fig_width,
            height = fig_height
          )
        } else if (significant_variable_only == F) {
          ggpubr::ggstripchart(
            data[["Data_renamed"]],
            x = "Group",
            y = colnames(data[["Data_renamed"]])[number],
            color = "Group",
            palette = ckey,
            order = order,
            size = size
          ) +
            ggplot2::stat_summary(
              fun.y = median,
              fun.ymin = median,
              fun.ymax = median,
              geom = "crossbar",
              width = 0.4,
              color = "Black",
              size = 0.2
            ) +
            ggplot2::scale_y_continuous(label = scientific_10x) +
            ggplot2::labs(
              title = NAMES[number],
              x = NULL,
              y = "Intensity"
            ) +
            ggplot2::theme_classic() + ggplot2::theme(
              plot.title = ggplot2::element_text(
                size = T_size,
                face = "bold",
                color = "Black",
                hjust = 0.5,
                lineheight = 1.2
              ),
              # title
              plot.subtitle = ggplot2::element_text(
                size = 15,
                hjust = 0.5
              ),
              # subtitle
              plot.caption = ggplot2::element_text(size = 15),
              # caption
              axis.title.x = ggplot2::element_text(
                vjust = 10,
                size = 15
              ),
              # X axis title
              axis.title.y = ggplot2::element_text(size = Y_text),
              # Y axis title
              axis.text.x = ggplot2::element_text(
                size = X_text,
                vjust = .5
              ),
              # X axis text
              axis.text.y = ggplot2::element_text(size = Y_lab),
              legend.position = legend_position
            )

          ggplot2::ggsave(
            filename = paste(NAMES[number], "dotplot.png", collapse = ""),
            path = paste0(getwd(), "/dotplot"),
            width = fig_width,
            height = fig_height
          )
        }
      }
    } else {
      p_val_data <- as.data.frame(p_val_data)
      if (asterisk == "t_test") {
        colnames(p_val_data) <- colnames(data[["Result"]])[1]
      } else if (asterisk == "u_test") {
        colnames(p_val_data) <- colnames(data[["Result"]])[2]
      }
      foreach::foreach(
        number = 1:nrow(p_val_data),
        .options.snow = opts,
        .packages = c("LMSstat", "dplyr", "plyr")
      ) %dopar% {
        df <- data_summary(data[["Data_renamed"]],
          varname = colnames(data[["Data_renamed"]])[number],
          groupnames = c("Group")
        )
        df <- df %>% plyr::mutate(ymax = len + SEM)
        colnames(df)[1] <- "Group"
        stat.test <- as.data.frame(matrix(
          data = NA,
          nrow = nrow(Comb),
          ncol = 2
        ))
        stat.test[, 1:2] <- Comb
        if (asterisk == "Dunn") {
          rownames(stat.test) <- paste0(
            stat.test[, 1],
            " - ", stat.test[, 2]
          )
        } else if (asterisk == "Scheffe") {
          rownames(stat.test) <- paste0(
            stat.test[, 2],
            "-", stat.test[, 1]
          )
        } else {
          rownames(stat.test) <- paste0(
            stat.test[, 1],
            "-", stat.test[, 2]
          )
        }
        colnames(p_val_data) <-
          as.data.frame(strsplit(
            colnames(p_val_data),
            "___"
          ))[1, ]
        t1 <-
          as.data.frame(t(as.data.frame(p_val_data[number, ])))
        rownames(t1) <- colnames(p_val_data)

        t2 <- stat.test
        t3 <- merge(t2, t1, by = 0)
        stat.test <- t3[, 2:4]
        colnames(stat.test) <- c(
          "group1", "group2",
          "p"
        )
        stat.test <-
          stat.test %>% plyr::mutate(
            p.adj.signif = case_when(
              p >
                sig_int[1] ~ "NS",
              p <= sig_int[1] & p > sig_int[2] ~ "*",
              p <= sig_int[2] ~ "**"
            )
          )
        stat.test <- stat.test[stat.test$p.adj.signif !=
          "NS", ]
        if (length(stat.test > 4)) {
          ggpubr::ggstripchart(
            data[["Data_renamed"]],
            x = "Group",
            y = colnames(data[["Data_renamed"]])[number],
            color = "Group",
            palette = ckey,
            order = order,
            size = size
          ) +
            ggplot2::stat_summary(
              fun.y = median,
              fun.ymin = median,
              fun.ymax = median,
              geom = "crossbar",
              width = 0.4,
              color = "Black",
              size = 0.2
            ) +
            ggplot2::scale_y_continuous(label = scientific_10x) +
            ggplot2::labs(
              title = NAMES[number],
              x = NULL,
              y = "Intensity"
            ) +
            ggplot2::theme_classic() + ggplot2::theme(
              plot.title = ggplot2::element_text(
                size = T_size,
                face = "bold",
                color = "Black",
                hjust = 0.5,
                lineheight = 1.2
              ),
              # title
              plot.subtitle = ggplot2::element_text(
                size = 15,
                hjust = 0.5
              ),
              # subtitle
              plot.caption = ggplot2::element_text(size = 15),
              # caption
              axis.title.x = ggplot2::element_text(
                vjust = 10,
                size = 15
              ),
              # X axis title
              axis.title.y = ggplot2::element_text(size = Y_text),
              # Y axis title
              axis.text.x = ggplot2::element_text(
                size = X_text,
                vjust = .5
              ),
              # X axis text
              axis.text.y = ggplot2::element_text(size = Y_lab),
              legend.position = legend_position
            ) +
            ggpubr::stat_pvalue_manual(
              stat.test,
              y.position = 1.05 * max(data[["Data_renamed"]][, number]),
              step.increase = step_increase,
              label.size = label_size,
              tip.length = tip_length,
              label = "p.adj.signif",
              size = 3.5,
              vjust = 0.05
            )
          ggplot2::ggsave(
            filename = paste(NAMES[number], "dotplot.png", collapse = ""),
            path = paste0(getwd(), "/dotplot"),
            width = fig_width,
            height = fig_height
          )
        } else if (significant_variable_only == F) {
          ggpubr::ggstripchart(
            data[["Data_renamed"]],
            x = "Group",
            y = colnames(data[["Data_renamed"]])[number],
            color = "Group",
            palette = ckey,
            order = order,
            size = size
          ) +
            ggplot2::stat_summary(
              fun.y = median,
              fun.ymin = median,
              fun.ymax = median,
              geom = "crossbar",
              width = 0.4,
              color = "Black",
              size = 0.2
            ) +
            ggplot2::scale_y_continuous(label = scientific_10x) +
            ggplot2::labs(
              title = NAMES[number],
              x = NULL,
              y = "Intensity"
            ) +
            ggplot2::theme_classic() + ggplot2::theme(
              plot.title = ggplot2::element_text(
                size = T_size,
                face = "bold",
                color = "Black",
                hjust = 0.5,
                lineheight = 1.2
              ),
              # title
              plot.subtitle = ggplot2::element_text(
                size = 15,
                hjust = 0.5
              ),
              # subtitle
              plot.caption = ggplot2::element_text(size = 15),
              # caption
              axis.title.x = ggplot2::element_text(
                vjust = 10,
                size = 15
              ),
              # X axis title
              axis.title.y = ggplot2::element_text(size = Y_text),
              # Y axis title
              axis.text.x = ggplot2::element_text(
                size = X_text,
                vjust = .5
              ),
              # X axis text
              axis.text.y = ggplot2::element_text(size = Y_lab),
              legend.position = legend_position
            )
          ggplot2::ggsave(
            filename = paste(NAMES[number], "dotplot.png", collapse = ""),
            path = paste0(getwd(), "/dotplot"),
            width = fig_width,
            height = fig_height
          )
        }
      }
    })
    close(pb)
  }, finally = {
    parallel::stopCluster(cl)
  })
}
