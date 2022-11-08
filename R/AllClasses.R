#' DeeDeeExperiment
#'
#' @rdname DeeDeeExperiment
#'
#' @exportClass DeeDeeExperiment
#'
#' @slot dea todo
#'
setClass("DeeDeeExperiment",
         contains = "RangedSummarizedExperiment",
         slots = representation(
           dea = "list"
         )
)

# TODOs in general: provide some corollary info as messages for most functions?
