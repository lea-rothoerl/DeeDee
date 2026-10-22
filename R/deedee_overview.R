#' DeeDee Overview Table
#'
#' @description `deedee_overview` creates a table containing information on the
#'              DEA contrasts in a DeeDeeExperiment object: The package that the
#'              DEA was run with, the number of genes, differentially expressed
#'              genes and number of up- as well as down-regulated genes for each
#'              contrast in the DDE.
#'
#' @param data instance of the DeeDeeExperiment class
#' @param pthresh threshold for adjusted p-values to be in-/excluded
#'                (default = 0.05)
#'
#' @return data.frame
#'
#' @examples
#'
#' data(dde_macrophage, package = "DeeDee")
#' deedee_overview(dde_macrophage)
#'
#' @export
#'

deedee_overview <- function(data, pthresh = 0.05) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_class(data, "DeeDeeExperiment")
  checkmate::assert_number(pthresh, lower = 0, upper = 1)

  # --------------------------------- summary ---------------------------------
  dea_info <- DeeDeeExperiment::getDEAInfo(data)
  dea_list <- deedee_from_dde(data)
  contrasts <- DeeDeeExperiment::getDEANames(data)

  dea_package <- vapply(contrasts, function(name) {
    pkg <- dea_info[[name]]$package
    if (is.null(pkg) || length(pkg) != 1 || is.na(pkg)) {
      NA_character_
    } else {
      as.character(pkg)
    }
  }, character(1))

  res <- data.frame(
    contrast = contrasts,
    package = unname(dea_package),
    genes = unname(vapply(dea_list, nrow, integer(1))),
    de_genes = unname(vapply(dea_list, function(d) {
      sum(d$padj < pthresh, na.rm = TRUE)
    }, integer(1))),
    up = unname(vapply(dea_list, function(d) {
      sum(d$padj < pthresh & d$logFC > 0, na.rm = TRUE)
    }, integer(1))),
    down = unname(vapply(dea_list, function(d) {
      sum(d$padj < pthresh & d$logFC < 0, na.rm = TRUE)
    }, integer(1))),
    stringsAsFactors = FALSE
  )

  # --------------------------------- return ----------------------------------
  return(res)
}
