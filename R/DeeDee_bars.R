#' Title
#'
#' @param dde todo
#'
#' @return todo
#' @export
#'
#' @examples
#' data("de_named_list", package = "DeeDee")
#' library(SummarizedExperiment)
#'
#' rd_macrophage <- DataFrame(
#'   gene_id = rownames(del$ifng_vs_naive))
#' rownames(rd_macrophage) <- rownames(del$ifng_vs_naive)
#' se_macrophage_noassays <- SummarizedExperiment(
#'   assays = SimpleList(),
#'   rowData = rd_macrophage
#' )
#' dde <- DeeDeeExperiment(
#'   se_macrophage_noassays,
#'   del
#' )
#' dde
ddedde_bars <- function(dde) {
  # extract the summary infos


  # put them together in a proper df


  # plot via gg



  # alternatively: return the df itself
}

