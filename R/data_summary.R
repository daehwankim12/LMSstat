#' Group-wise numeric summary
#'
#' @param data       Data frame containing the variables.
#' @param varname    Name (string) of the numeric variable to summarise.
#' @param groupnames Character vector with the grouping variables.
#' @return A data frame with N, mean, and SD for each group.
#' @importFrom stats sd
#' @export
data_summary <- function(data, varname, groupnames) {
  summary_func <- function(x, col) {
    c(
      mean = mean(x[[col]], na.rm = TRUE),
      SEM = sd(x[[col]] / sqrt(length(x[[col]])), na.rm = TRUE)
    )
  }
  data_sum <- plyr::ddply(data, groupnames,
    .fun = summary_func,
    varname
  )
  data_sum <- plyr::rename(data_sum, c("mean" = "len"))
  return(data_sum)
}
