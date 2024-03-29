#' @title DeeDee Heatmap
#'
#' @description `deedee_heatmap` creates a heatmap depicting the logFC as a
#' measure of the differential expression of the first genes in the given
#' datasets.
#'
#' @param data named list of results from deedee_prepare()
#' @param show_first indicating the number of genes depicted (default = 25)
#' @param show_gene_names boolean, show row names next to heatmap
#'                        (default = FALSE)
#' @param dist select the distance measure (`euclidean`, `manhattan`, `pearson`,
#'             `spearman`)
#' @param clust select the clustering method (`single`, `complete`, `average`,
#'              `centroid`)
#' @param show_na boolean, include genes with NA values in heatmap
#'                (default = FALSE)
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#'
#' @return Heatmap object (plottable with show()/print()). The resulting heatmap
#' (`res`) can be opened in a Shiny App with interactive functionality by
#' running the command `ht_shiny(res)` (see package `InteractiveComplexHeatmap`).
#'
#' @examples
#'
#' data(DE_results_IFNg_naive, package = "DeeDee")
#' IFNg_naive <- deedee_prepare(IFNg_naive, "DESeq2")
#'
#' data(DE_results_IFNg_both, package = "DeeDee")
#' IFNg_both <- deedee_prepare(IFNg_both, "DESeq2")
#'
#' data(DE_results_Salm_naive, package = "DeeDee")
#' Salm_naive <- deedee_prepare(Salm_naive, "DESeq2")
#'
#' data(DE_results_Salm_both, package = "DeeDee")
#' Salm_both <- deedee_prepare(Salm_both, "DESeq2")
#'
#' dd_list <- list(
#'   IFNg_naive = IFNg_naive, IFNg_both = IFNg_both,
#'   Salm_naive = Salm_naive, Salm_both = Salm_both
#' )
#'
#' deedee_heatmap(dd_list,
#'   pthresh = 0.05, show_first = 25,
#'   show_gene_names = FALSE, dist = "euclidean",
#'   clust = "average", show_na = FALSE
#' )
#' @export
#'

deedee_heatmap <- function(data,
                           show_first = 25,
                           show_gene_names = FALSE,
                           dist = "euclidean",
                           clust = "average",
                           show_na = FALSE,
                           pthresh = 0.05) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  checkmate::assert_number(show_first, lower = 1)
  checkmate::assert_logical(show_gene_names)
  choices1 <- c("euclidean", "manhattan", "pearson", "spearman")
  checkmate::assert_choice(dist, choices1)
  choices2 <- c("single", "complete", "average", "centroid")
  checkmate::assert_choice(clust, choices2)
  checkmate::assert_logical(show_na)

  # ---------------------------- data preparation -----------------------------
  for (i in 1:length(data)) {
    data[i][[1]] <- subset(data[i][[1]], data[i][[1]]$pval < pthresh) # pthresh
    data[i][[1]] <- data[i][[1]]["logFC"] # removing p-value column
    colnames(data[i][[1]]) <- names(data)[i] # creating unique colnames
    data[i][[1]] <- tibble::rownames_to_column(data[i][[1]])
  }

  # ----- creation of a comp matrix with the modified results from above ------
  comp <- dplyr::full_join(data[1][[1]],
    data[2][[1]],
    by = "rowname",
    copy = FALSE
  )
  if (length(data) > 2) {
    for (i in 3:length(data)) {
      comp <- dplyr::full_join(comp, data[i][[1]], by = "rowname", copy = FALSE)
    }
  }

  row.names(comp) <- comp$rowname
  comp <- subset(comp, select = -c(rowname)) # removing column with rownames
  comp <- comp[rowSums(!is.na(comp)) >= floor(length(comp) / 2) + 1, ]
  if (show_na == FALSE) {
    comp <- comp[stats::complete.cases(comp[colnames(comp)]), ]
  }
  comp <- as.matrix(comp)

  if (length(comp[, 1]) == 0) {
    return(NULL)
  }

  # ------------------- creation of the resulting heatmap ---------------------
  if (show_gene_names == FALSE) {
    rownames(comp) <- c()
  }

  col <- viridis::viridis(n = 15, option = "magma")

  res <- ComplexHeatmap::Heatmap(comp[1:min(
    show_first,
    length(comp[, 1])
  ), ],
  name = "logFC",
  col = col,
  clustering_distance_rows = dist,
  clustering_method_rows = clust
  )

  # ComplexHeatmap::draw(res)

  # --------------------------------- return ----------------------------------
  return(res)
}
