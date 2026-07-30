#' A sample `DeeDeeExperiment` object
#'
#' A sample `DeeDeeExperiment` object, built around a `DESeqDataSet` from the
#' `macrophage` package, containing four DESeq2 contrasts as differential
#' expression results.
#'
#' @details The underlying `DESeqDataSet` compares macrophages under four
#' conditions: naive, IFNg-stimulated, Salmonella-infected, and both
#' IFNg-stimulated and Salmonella-infected. The `dea` slot of this object
#' contains four contrasts:
#' \itemize{
#'   \item `IFNg_naive`: IFNg vs. naive
#'   \item `IFNg_both`: IFNg + Salmonella vs. IFNg
#'   \item `Salm_naive`: Salmonella vs. naive
#'   \item `Salm_both`: IFNg + Salmonella vs. Salmonella
#' }
#'
#' The code to create this object can be found in `/inst/scripts/generate_data.R`
#' in the DeeDee package.
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name dde_macrophage
#' @docType data
NULL
