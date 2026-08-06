#' DeeDee Prepare
#'
#' @description `deedee_from_dde` prepares an instance of the DeeDeeExperiment
#' class for the usage in other DeeDee functions.
#'
#' @param dde Instance of the DeeDeeExperiment class
#'
#' @return list of DeeDee tables, to be used as part of the input for the other
#' DeeDee functions
#'

deedee_from_dde <- function(dde) {
  dea_names <- names(dde@dea)
  symbols <- SingleCellExperiment::rowData(dde)$SYMBOL

  deedee_list <- lapply(dea_names, function(name) {
    res <- DeeDeeExperiment::getDEA(dde, name)
    df <- data.frame(
      logFC = res[[1]],
      pval  = res[[2]],
      row.names = rownames(res)
    )
    if (!is.null(symbols)) {
      df$symbol <- symbols[rownames(res)]
    }
    df$logFC <- as.numeric(df$logFC)
    df$pval  <- as.numeric(df$pval)
    df
  })

  names(deedee_list) <- dea_names
  deedee_list
}
