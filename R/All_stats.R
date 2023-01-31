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
All_stats <-
  function(Data,
           Adjust_p_value = TRUE,
           Adjust_method = "BH") {
    # If the dataset has a column named "group", rename it to "Group"
    if ("group" %in% colnames(Data)) {
      Data <- dplyr::rename(Data, "Group" = "group")
    }

    # If the dataset has a column named "sample", rename it to "Sample"
    if ("sample" %in% colnames(Data)) {
      Data <- dplyr::rename(Data, "Sample" = "sample")
    }

    # Convert the "Group" column to character type
    Data$Group <- as.character(Data$Group)

    # Sort the data by "Group"
    Data <- Data %>% dplyr::arrange(Data$Group)

    # Create a new data frame with modified column names
    Data_renamed <- Data
    nmet <- ncol(Data) - 2 # number of variables (metabolites)
    nmet_seq <-
      seq_len(nmet) + 2 # sequence of column indices for variables
    colnames(Data_renamed) <-
      c(colnames(Data[1:2]), paste0("V", 1:nmet))
    rownames(Data_renamed) <-
      Data[, 1] # set row names to sample names

    # Extract the variables (metabolites) as a separate data frame
    Data_renamed_raw <- Data_renamed[, -c(1, 2)]

    # Convert the variables to numeric type
    Data_renamed_raw <- apply(Data_renamed_raw, 2, as.numeric)

    # Combine the variables with sample names and group information
    Data_final <- cbind(Data[, 1:2], Data_renamed_raw)

    # Convert the combined data frame to a data table
    Data_final <- data.table::as.data.table(Data_final)

    # Convert the "Group" column to a factor
    Data_final$Group <- as.factor(as.character(Data_final$Group))

    # Split the data table by group
    groups_split <-
      split(Data_final, Data_final$Group)

    # Check if there are more than two groups
    group_nottwo <-
      length(unique(Data_final$Group)) > 2

    # Remove the intermediate data frame
    rm(Data_renamed_raw)

    #### ttest ####

    # Function to perform t-test on two groups for a given variable (i)
    split_t_test <- function(x, i) {
      # If both groups have zero variance, return p-value of 1
      if (var(groups_split[[x[1]]][[i]]) == 0 &&
          var(groups_split[[x[2]]][[i]]) == 0) {
        1
      } else {
        # Otherwise, perform t-test and return p-value
        t.test(groups_split[[x[1]]][[i]], groups_split[[x[2]]][[i]])[["p.value"]]
      }
    }

    # Perform t-test on all pairs of groups for each variable
    res_ttest <-
      lapply((nmet_seq), function(i) {
        as.list(combn(names(groups_split), 2, split_t_test, i = i))
      })

    # Bind the results into a single data frame
    df_ttest <-
      data.table::rbindlist(res_ttest)
    df_ttest <- data.frame(df_ttest)

    # Print message
    print("T-test has finished")

    #### utest ####

    # Function to perform u-test (Wilcoxon rank sum test) on two groups for a given variable (i)
    split_u_test <- function(x, i) {
      # If both groups have zero variance, return p-value of 1
      if (var(groups_split[[x[1]]][[i]]) == 0 &&
          var(groups_split[[x[2]]][[i]]) == 0) {
        1
      } else {
        # Otherwise, perform u-test and return p-value
        wilcox.test(groups_split[[x[1]]][[i]], groups_split[[x[2]]][[i]])[["p.value"]]
      }
    }

    # Perform u-test on all pairs of groups for each variable
    res_utest <-
      lapply((nmet_seq), function(i) {
        as.list(combn(names(groups_split), 2, split_u_test, i = i))
      })

    # Bind the results into a single data frame
    df_utest <-
      data.table::rbindlist(res_utest)
    df_utest <- data.frame(df_utest)

    # Print message
    print("U-test has finished")

    if (group_nottwo) {
      #### ANOVAPostHoc ####
      # Perform ANOVA and Scheffe's post-hoc test
      formula_anova <-
        lapply(colnames(Data_final)[nmet_seq], function(x) {
          # Create a formula for each metabolite with the metabolite as the response variable and group as the predictor
          as.formula(paste0(x, " ~ Group"))
        })
      res_anova <-
        lapply(formula_anova, function(x) {
          # Apply ANOVA to each formula
          summary(aov(x, data = Data_final))
        })
      names(res_anova) <-
        format(formula_anova)
      p_anova <-
        unlist(lapply(res_anova, function(x) {
          # Extract the p-value for each ANOVA test
          x[[1]]$"Pr(>F)"[1]
        }))

      # Create a data frame to store the p-values
      df_anova <-
        data.table::data.table(p_anova)
      df_anova <- data.frame(df_anova)

      anovapost_name <-
        lapply(colnames(Data_final)[nmet_seq], function(x) {
          # Create a formula for each metabolite with the metabolite as the response variable and group as the predictor
          as.formula(paste0(x, " ~ Group"))
        })
      res_anovapost <-
        lapply(anovapost_name, function(x) {
          # Perform Scheffe's post-hoc test on each formula
          DescTools::PostHocTest(aov(x, data = Data_final), method = "scheffe")
        })
      names(res_anovapost) <-
        format(anovapost_name)
      post_anova <-
        lapply(res_anovapost, function(x) {
          # Extract the p-values for each post-hoc test
          x[["Group"]][, 4]
        })

      # Create a data frame to store the p-values for the post-hoc tests
      df_anova_post <-
        t(data.frame(post_anova))

      # Print message indicating ANOVA and post-hoc tests are finished
      print("Anova & PostHoc has finished")


      #### KruskalWallisPostHoc ####
      # Perform Kruskal-Wallis test and post-hoc analysis using the Dunn test

      # Create formulas for each metabolite to be used in the Kruskal-Wallis test
      formula_kw <-
        lapply(colnames(Data_final)[nmet_seq], function(x) {
          as.formula(paste0(x, " ~ Group"))
        })

      # Perform Kruskal-Wallis test for each metabolite and store p-values in a list
      res_kw <-
        lapply(formula_anova, function(x) {
          kruskal.test(x, data = Data_final)[["p.value"]]
        })

      # Name the list elements with the formulas used in the Kruskal-Wallis test
      names(res_kw) <-
        format(formula_kw)

      # Unlist the p-values and store them in a data frame
      p_kw <- unlist(res_kw)
      df_kw <-
        data.table::data.table(p_kw)
      df_kw <- data.frame(df_kw)

      # Perform post-hoc analysis using the Dunn test for each metabolite
      # and store the results in a list
      kwpost_name <-
        lapply(colnames(Data_final)[nmet_seq], function(x) {
          as.formula(paste0(x, " ~ Group"))
        })
      res_kwpost <-
        lapply(kwpost_name, function(x) {
          FSA::dunnTest(x, data = Data_final, method = "none")
        })

      # Name the list elements with the formulas used in the Dunn test
      names(res_kwpost) <-
        format(kwpost_name)

      # Extract the unadjusted p-values for each pairwise comparison and store them in a list
      post_kw <-
        lapply(res_kwpost, function(x) {
          x[["res"]][["P.unadj"]]
        })

      # Adjust the p-values using the Benjamini-Hochberg procedure and store them in a data frame
      post_kw <-
        lapply(post_kw, function(x) {
          p.adjust(x, method = "BH")
        })
      df_kw_post <-
        t(data.table::data.table(data.frame(post_kw)))
      df_kw_post <-
        data.frame(df_kw_post)

      # Print a message indicating that the Kruskal-Wallis test and post-hoc analysis have finished
      print("Kruskal Wallis & PostHoc has finished")
    }

    #### finalization ####
    # Finalize the results of the statistical tests and rename the rows and columns of the data frames

    # Initialize the Names variable to NULL
    Names <- NULL

    # In R 4.2.2, remove if-else. (Legacy)
    # If there are more than two groups, create all possible pairs of group names
    # and store them in the Names matrix
    # if (group_nottwo) {
    #   for (i in seq_len(choose(length(groups_split), 2))) {
    #     Names <- rbind(Names, paste(combn(names(groups_split), 2)[1, i],
    #                                 combn(names(groups_split), 2)[2, i],
    #                                 sep = "-"
    #     ))
    #   }
    #
    #   # If there are only two groups, create a single pair of group names
    #   # and store it in the Names matrix
    # } else {
    #   Names <- rbind(Names, paste(combn(names(groups_split), 2)[1],
    #                               combn(names(groups_split), 2)[2],
    #                               sep = "-"
    #   ))
    # }

    # Change the row names of the data frames containing the results of the t-test and U-test
    # to the names of the metabolites
    rownamechange <-
      colnames(Data)[nmet_seq]
    rownames(df_ttest) <- rownamechange
    rownames(df_utest) <- rownamechange

    # Change the column names of the data frames containing the results of the t-test and U-test
    # to include the names of the groups and the type of test
    colnames(df_ttest) <-
      paste(Names[, 1], "t-test", sep = "___")
    colnames(df_utest) <-
      paste(Names[, 1], "u-test", sep = "___")
    # If there are more than two groups, change the row names of the data frames containing
    # the results of the ANOVA, ANOVA post-hoc, Kruskal-Wallis, and Kruskal-Wallis post-hoc tests
    # to the names of the metabolites
    if (group_nottwo) {
      AN_Post_names <- colnames(df_anova_post)
      DU_post_names <-
        res_kwpost[[1]]$res$Comparison

      rownames(df_anova) <-
        rownamechange
      rownames(df_anova_post) <-
        rownamechange
      rownames(df_kw) <- rownamechange
      rownames(df_kw_post) <-
        rownamechange

      # Change the column names of the data frames containing the results of the ANOVA and ANOVA post-hoc tests
      # to include the names of the groups and the type of test
      colnames(df_anova) <- "Anova"
      colnames(df_anova_post) <-
        paste(AN_Post_names, "ANO_posthoc", sep = "___")

      # Change the column names of the data frames containing the results of the Kruskal-Wallis and Kruskal-Wallis post-hoc tests
      # to include the names of the groups and the type of test
      colnames(df_kw) <-
        "Kruskal_Wallis"
      colnames(df_kw_post) <-
        paste(DU_post_names, "Kru_posthoc(Dunn)", sep = "___")
    }

    # Remove the nmet_seq variable
    rm(nmet_seq)

    # If p-value adjustment is enabled, adjust the p-values in the data frames
    # containing the results of the t-test, U-test, ANOVA, and Kruskal-Wallis tests
    # using the specified method
    if (Adjust_p_value == TRUE) {
      print("###########################################")
      print(paste0("adjusted according to the ", Adjust_method, " method"))
      print("###########################################")

      # Define a function that adjusts p-values using the specified method
      adj_func <- function(x) {
        p.adjust(x, method = Adjust_method)
      }

      # Adjust the p-values in the data frames using the adj_func function
      df_ttest <-
        apply(df_ttest, 2, adj_func)
      df_utest <-
        apply(df_utest, 2, adj_func)
      if (group_nottwo) {
        df_anova <- apply(df_anova, 2, adj_func)
        df_kw <-
          apply(df_kw, 2, adj_func)
      }

      # If p-value adjustment is not enabled, print a message indicating that the p-values are not adjusted
    } else {
      print("###########################################")
      print("p_value not adjusted")
      print("###########################################")
    }

    # Print a message indicating that the statistical tests are complete
    print("statistical test has finished")

    # Combine the data frames containing the results of the t-test and U-test into a single data frame
    Result <- cbind(df_ttest, df_utest)

    # If there are more than two groups, combine the data frames containing the results of the ANOVA, ANOVA post-hoc,
    # Kruskal-Wallis, and Kruskal-Wallis post-hoc tests into the Result data frame
    if (group_nottwo) {
      Result <- cbind(Result, df_anova, df_anova_post, df_kw, df_kw_post)
    }
    # Convert the Result data frame to a matrix
    Result <- as.matrix(Result)

    # If there are more than two groups, convert the data frames containing the results of the ANOVA,
    # Kruskal-Wallis, and Kruskal-Wallis post-hoc tests to data frames or matrices
    if (group_nottwo) {
      df_anova <- as.data.frame(df_anova)
      df_kw <- as.data.frame(df_kw)
      df_kw_post <-
        as.matrix(df_kw_post)
    }

    # Create a list containing the input data, the data with renamed columns, the results of the statistical tests,
    # and the results of the ANOVA, ANOVA post-hoc, Kruskal-Wallis, and Kruskal-Wallis post-hoc tests (if applicable)
    Final <- list()
    Final$Data <- Data
    Final$Data_renamed <- Data_renamed
    Final$Result <- Result
    if (group_nottwo) {
      Final$Anova <- df_anova
      Final$Anova_PostHoc <-
        df_anova_post
      Final$KW <- df_kw
      Final$Dunn <- df_kw_post
      Final$t_test <- df_ttest
      Final$u_test <- df_utest
    } else {
      # If there are only two groups, include the results of the t-test and U-test in the Final list
      Final$t_test <- df_ttest[, 1]
      Final$u_test <- df_utest[, 1]
    }

    # Return the Final list
    Final
  }
