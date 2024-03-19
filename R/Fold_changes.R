#' Compute Fold Changes and P-values for Grouped Data
#'
#' Given a dataset with groupings, this function computes the fold changes
#' and p-values between groups of data starting from the third column onward.
#'
#' @param data A data frame where the second column represents the grouping
#' variable, and subsequent columns represent numeric data to be analyzed.
#'
#' @return A data frame containing the original data, appended with columns
#' representing the fold change and p-values for each group combination.
#'
#' @examples
#' # data(Data)
#' # result <- Fold_changes(Data)
#' @export
Fold_changes <- function(data) {
  # Ensure the second column is used to order the data and derive groups
  data_col_names <- colnames(data)

  sorted_indices <- order(data[[2]])
  data <- data[sorted_indices, ]

  group <- as.factor(as.character(data[[2]]))

  # Extract numeric data starting from the third column
  raw_data <- data.frame(lapply(data[3:ncol(data)], as.numeric))

  # Split raw_data by group
  groups_split <- split(raw_data, group)

  # Generate group combinations
  group_combinations <-
    combn(names(groups_split), 2, simplify = FALSE)

  # Compute fold change and p-values
  compute_fc_and_pvals <- function(x) {
    means1 <- colMeans(groups_split[[x[1]]], na.rm = TRUE)
    means2 <- colMeans(groups_split[[x[2]]], na.rm = TRUE)
    fc <- means1 / means2
    pval <-
      matrixTests::col_t_equalvar(groups_split[[x[1]]], groups_split[[x[2]]])$pvalue
    list(fc = fc, pval = pval)
  }

  results <- lapply(group_combinations, compute_fc_and_pvals)

  fc_matrix <- do.call(cbind, lapply(results, "[[", "fc"))
  pval_matrix <- do.call(cbind, lapply(results, "[[", "pval"))

  # Naming helper function
  rename_columns <- function(mat, separator, suffix) {
    new_names <-
      apply(do.call(rbind, group_combinations), 1, function(row) {
        paste(row, collapse = separator)
      })
    colnames(mat) <- paste(new_names, suffix)
    return(mat)
  }

  fc_matrix <- t(rename_columns(fc_matrix, " / ", " fold change"))
  pval_matrix <- t(rename_columns(pval_matrix, " - ", " p-value"))

  fc_matrix <- cbind(row.names(fc_matrix), row.names(fc_matrix), fc_matrix)
  pval_matrix <- cbind(row.names(pval_matrix), row.names(pval_matrix), pval_matrix)

  colnames(fc_matrix) <- data_col_names
  colnames(pval_matrix) <- data_col_names

  # Bind the results together
  result <- rbind(data, fc_matrix, pval_matrix)

  return(result)
}
