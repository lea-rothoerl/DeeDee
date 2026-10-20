#' DeeDee
#'
#' `DeeDee` is a package that allows for the visual comparison of two or more
#'   DEA results.
#'
#' `DeeDee` contains functions for the preprocessing of raw DEA result data and
#'   multiple functions that create graphical output and give an impression of
#'   similarities and differences in the expression profiles. The functions are
#'   contained individually as well as included in an interactive Shiny Web App.
#'
#' @name DeeDee-pkg
"_PACKAGE"
NULL

utils::globalVariables(
  c("logFC1", "logFC2", "rowname", "hover_text", "contrast", "direction", "n")
)
