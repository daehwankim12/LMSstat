#' Automatically processes T-test, U-test, Anova, Scheffe(Anova Post-Hoc),
#' Krukal Wallis, Dunn-test(BH adjusted,(Kurkal Wallis Post-Hoc)) while allowing adjustment of FDR
#'
#' @param Data csv file with Header as False First column with Sample Second column
#' with Multilevel(Mixomics) so that it can be compatible with other multivariate
#' statistics Third column with Group information. Rest of the columns are the metabolites to be tested.
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
    # Rename the first and second columns to "Group" and "Sample"
    colnames(Data)[1:2] <- c("Sample", "Group")

    # Convert the "Group" column to character type
    Data$Group <- as.character(Data$Group)

    # Sort the data by "Group"
    Data <- dplyr::arrange(Data, Data$Group)
    Data_ori <- Data

    # Modify column names and convert variables to numeric
    nmet <- ncol(Data) - 2 # number of variables (metabolites)
    nmet_seq <-
      seq_len(nmet) # sequence of column indices for variables


    # Convert relevant columns to numeric using lapply (vectorization)
    cols_to_convert <- 3:ncol(Data)
    Data[, cols_to_convert] <-
      lapply(Data[, cols_to_convert], as.numeric)

    # Convert the data frame to a data table
    Data_final <- data.table::as.data.table(Data)
    Data_final_raw <- Data_final[, -c("Sample", "Group")]

    # Convert the "Group" column to a factor
    Data_final[, Group := as.factor(as.character(Group))]

    # Split the data table by group
    groups_split <- split(Data_final_raw, Data_final$Group)

    # Check if there are more than two groups
    group_nottwo <- length(unique(Data_final$Group)) > 2

    variances <-
      lapply(groups_split, function(group) {
        lapply(group, var)
      })
    Data_final <- as.data.frame(Data_final)

    #### ttest ####
    # Generate a matrix for each group, where each column is a variable (assuming groups_split is a list of data.frames)
    group_matrices <- lapply(groups_split, as.matrix)

    # Helper function to perform vectorized t-test between two matrices with unequal variances
    vectorized_t_test <- function(mat1, mat2) {
      n1 <- nrow(mat1)
      n2 <- nrow(mat2)

      means1 <- colMeans(mat1)
      means2 <- colMeans(mat2)
      vars1 <- apply(mat1, 2, var)
      vars2 <- apply(mat2, 2, var)

      # Using formula for t-statistic under the assumption of unequal variances (Welch's t-test)
      t_stats <- (means1 - means2) / sqrt(vars1 / n1 + vars2 / n2)

      # Degrees of freedom for Welch’s t-test (unequal variances)
      df <-
        (((vars1 / n1) + (vars2 / n2))^2) / ((vars1 / n1)^2 / (n1 - 1) + (vars2 /
          n2)^
          2 / (n2 - 1))
      p_values <- 2 * (1 - pt(abs(t_stats), df = df))

      return(p_values)
    }

    # Compute all combinations of groups
    group_combinations <-
      combn(names(groups_split), 2, simplify = FALSE)

    # Perform vectorized t-tests for each combination of groups
    result_list_t <- lapply(group_combinations, function(combo) {
      mat1 <- group_matrices[[combo[1]]]
      mat2 <- group_matrices[[combo[2]]]
      vectorized_t_test(mat1, mat2)
    })


    # Convert results into desired format (data frame)
    df_ttest <- t(data.frame(do.call(rbind, result_list_t)))

    print("T-test has finished")

    #### utest ####
    # Helper function to perform vectorized Wilcoxon test between two matrices
    vectorized_u_test <- function(mat1, mat2) {
      result <- matrixTests::col_wilcoxon_twosample(mat1, mat2)
      return(result$pvalue)
    }

    # Perform vectorized Wilcoxon tests for each combination of groups
    result_list_u <- lapply(group_combinations, function(combo) {
      mat1 <- group_matrices[[combo[1]]]
      mat2 <- group_matrices[[combo[2]]]
      vectorized_u_test(mat1, mat2)
    })

    # Convert Wilcoxon test results into desired format (data frame)
    df_utest <- t(data.frame(do.call(rbind, result_list_u)))

    print("U-test has finished")


    if (group_nottwo) {
      metabolite_names <- colnames(Data_final)[nmet_seq + 2]

      # Setting up the cluster
      num_core <- parallel::detectCores() - 1
      final_cores <- min(5, num_core)
      cl <- parallel::makeCluster(final_cores)
      on.exit(parallel::stopCluster(cl))

      perform_anova_tests <- function(colname) {
        data <- Data_final[, c("Group", colname)]
        formula_str <- as.formula(paste0(colname, " ~ Group"))
        model <- aov(formula_str, data = data)
        anova_res <- summary(model)[[1]]$"Pr(>F)"[1]
        posthoc_res <-
          DescTools::PostHocTest(model, method = "scheffe")$Group[, 4]
        list(anova_res = anova_res, posthoc_res = posthoc_res)
      }

      results <-
        parallel::parLapply(cl, metabolite_names, perform_anova_tests)

      p_anova <- sapply(results, function(x) {
        x$anova_res
      })
      df_anova <- data.frame(p_anova)
      df_anova_post <- do.call(rbind, lapply(results, function(x) {
        x$posthoc_res
      }))

      print("Anova & PostHoc has finished")

      #### KruskalWallisPostHoc ####
      # Use matrixTests to perform Kruskal-Wallis tests for each column
      kw_results <-
        matrixTests::col_kruskalwallis(as.matrix(Data_final[metabolite_names]), Data_final$Group)

      perform_posthoc_tests <- function(colname) {
        formula <- as.formula(paste0(colname, " ~ Group"))
        data_subset <- Data_final[, c("Group", colname)]
        dunn_res <-
          FSA::dunnTest(formula, data = data_subset, method = "none")
        post_pvals <- dunn_res[["res"]][["P.unadj"]]
        adjusted_pvals <- p.adjust(post_pvals, method = "BH")
        list(adjusted_pvals = adjusted_pvals, dunn_res = dunn_res)
      }

      # Use parLapply to perform the post-hoc tests in parallel
      results_kw_posthoc <-
        parallel::parLapply(cl, metabolite_names, perform_posthoc_tests)

      # Extract the Kruskal-Wallis results and put them into a dataframe
      df_kw <- data.frame(p_kw = kw_results$pvalue)

      # Extract the adjusted p-values from the posthoc tests and combine them into a dataframe
      df_kw_post <- do.call(
        rbind,
        lapply(results_kw_posthoc, function(x) {
          x$adjusted_pvals
        })
      )

      print("Kruskal Wallis & PostHoc has finished")
    }


    #### finalization ####
    # Finalize the results of the statistical tests and rename the rows and columns of the data frames

    # Initialize the Names variable to NULL
    Names <- NULL

    for (i in seq_len(choose(length(groups_split), 2))) {
      Names <- rbind(Names, paste(combn(names(groups_split), 2)[1, i],
        combn(names(groups_split), 2)[2, i],
        sep = "-"
      ))
    }

    # Change the row names of the data frames containing the results of the t-test and U-test
    # to the names of the metabolites
    rownamechange <-
      colnames(Data)[nmet_seq + 2]
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
      AN_post_names <- colnames(df_anova_post)
      DU_post_names <-
        results_kw_posthoc[[1]]$dunn_res$res$Comparison

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
        paste(AN_post_names, "ANO_posthoc", sep = "___")

      # Change the column names of the data frames containing the results of the Kruskal-Wallis and Kruskal-Wallis post-hoc tests
      # to include the names of the groups and the type of test
      colnames(df_kw) <-
        "Kruskal_Wallis"
      colnames(df_kw_post) <-
        paste(DU_post_names, "Kru_posthoc(Dunn)", sep = "___")
    }

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
    Final$Data <- Data_ori
    Final$Result <- Result
    if (group_nottwo) {
      Final$Anova <- df_anova
      Final$Anova_PostHoc <- df_anova_post
      Final$KW <- df_kw
      Final$Dunn <- df_kw_post
      Final$t_test <- as.data.frame(df_ttest)
      Final$u_test <- as.data.frame(df_utest)
    } else {
      # If there are only two groups, include the results of the t-test and U-test in the Final list
      Final$t_test <- df_ttest[, 1]
      Final$u_test <- df_utest[, 1]
    }

    # Return the Final list
    Final
  }
