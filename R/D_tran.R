#' Scale_Transformation of data file
#'
#' @param Data uploaded data
#' @param param "None","Auto","log10","Pareto"
#' @param save T,F
#'
#' @return scaled, transformed CSV file
#' @export
#'
#' @examples data(Data)
#' D_tran(Data, param = "Pareto")
D_tran <- function(data, param = "None", save = FALSE) {
  data <- data %>% dplyr::arrange(data$Group)
  rownames(data) <- data[, 1]
  data <- data[, -1]
  data[, 2:ncol(data)] <-
    sapply(data[, 2:ncol(data)], function(x) {
      as.numeric(x)
    })

  if (param == "None") {
    data <- data
  } else if (param == "Auto") {
    data[, 2:ncol(data)] <- scale(data[, 2:ncol(data)], scale = T)
  } else if (param == "log10") {
    data[, 2:ncol(data)] <- log10(data[, 2:ncol(data)])
  } else if (param == "Pareto") {
    data[, 2:ncol(data)] <-
      scale(data[, 2:ncol(data)], scale = sqrt(sapply(data[, 2:ncol(data)], sd)))
  }
  if (save) {
    write.csv(data, paste0("datafile_scaled_to_", param, ".csv"))
  }
  return(data)
}
