scientific_10x <- function(values) {
  s_values <- sprintf("%.1e", values)
  coefs <- sub("e[+|-]?\\d+$", "", s_values)
  exps <- gsub("^.*e([+|-]?\\d+)$", "\\1", s_values)
  exps <- sub("^\\+", "", exps)
  formatted <- mapply(function(coef, exp) {
    coef_decimal <- ifelse(!grepl("\\.", coef), paste0(coef, ".0"), coef)
    parse(text = paste0(coef_decimal, " %*% 10^", exp))
  }, coefs, exps)

  return(formatted)
}
