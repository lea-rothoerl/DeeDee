#' deedee_heatmap
#'
#' Construct a logFC heatmap on a subset of features
#'
#' @param dde A [DeeDeeExperiment] object.
#' @param show_first Numeric value, specifying the number of features to include.
#' The ranking is done according to ... TODO
#' @param show_gene_names Logical value, whether to display the gene names on the
#' heatmap's side.
#' @param dist Character value, specifying the distance type to use in the call
#' to `ComplexHeatmap`. Accordingly, it can be a pre-defined character
#' ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",
#' "pearson", "spearman", "kendall").
#' @param clust Character value. Defines the method to perform hierarchical
#' clustering, passed to hclust, as used in `ComplexHeatmap`.
#' @param show_na Logical value, whether to include features that have NA value
#' for the log fold change.
#' @param pthresh Numeric value, corresponding to the p-value to use as a
#' threshold to subset the features to include.
#'
#' @return  A plot object generated via the `ComplexHeatmap` package.
#' @export
#'
#' @examples
#' data("de_named_list", package = "DeeDee")
#' library("SummarizedExperiment")
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
#'   de_results = del
#' )
#' dde
#'
#' deedee_heatmap(dde,
#'                pthresh = 0.05, show_first = 25,
#'                show_gene_names = FALSE, dist = "euclidean",
#'                clust = "average", show_na = FALSE
#' )
#'
deedee_heatmap <- function(dde,
                           show_first = 25,
                           show_gene_names = FALSE,
                           dist = "euclidean",
                           clust = "average",
                           show_na = FALSE,
                           pthresh = 0.05) {
  if (length(dea(dde)) < 2)
    stop("Please provide a dde object including at least two DE analyses")

  if (!is.numeric(show_first) | !(show_first > 0))
    stop("Please specify a positive integer for the show_first parameter")

  if (!is.numeric(pthresh) | !(pthresh > 0 & pthresh <= 1))
    stop("Please specify a valid p-value threshold")

  stopifnot(is.logical(show_gene_names))
  stopifnot(is.logical(show_na))

  stopifnot(is.character(dist))
  stopifnot(is.character(clust))


  dea_list <- get_dea_list(dde)
  # and then proceed as for the old implementation

  for (i in seq_len(length(dea_list))) {
    dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$padj < pthresh) # pthresh
    dea_list[[i]] <- dea_list[[i]]["log2FoldChange"] # removing p-value column
    colnames(dea_list[[i]]) <- names(dea_list)[[i]] # creating unique colnames
    dea_list[[i]] <- tibble::rownames_to_column(dea_list[[i]])
  }

  comp <- dplyr::full_join(dea_list[[1]],
                           dea_list[[2]],
                           by = "rowname",
                           copy = FALSE
  )
  if (length(dea_list) > 2) {
    for (i in 3:length(dea_list)) {
      comp <- dplyr::full_join(comp, dea_list[[i]], by = "rowname", copy = FALSE)
    }
  }


  row.names(comp) <- comp$rowname

  # comp <- subset(comp, select = -c(rowname)) # removing column with rownames
  comp <- comp[, !(names(comp) %in% "rowname")]

  comp <- comp[rowSums(!is.na(comp)) >= floor(length(comp) / 2) + 1, ]
  if (show_na == FALSE) {
    comp <- comp[stats::complete.cases(comp[colnames(comp)]), ]
  }
  comp <- as.matrix(comp)

  if (length(comp[, 1]) == 0) {
    return(NULL)
  }

  if (show_gene_names == FALSE) {
    rownames(comp) <- c()
  }

  col <- viridis::viridis(n = 15, option = "magma")

  res <- ComplexHeatmap::Heatmap(
    comp[seq_len(min(show_first, length(comp[, 1]))), ],
  name = "logFC",
  col = col,
  clustering_distance_rows = dist,
  clustering_method_rows = clust
  )

  # ComplexHeatmap::draw(res)

  return(res)

}



#' deedee_heatmap(dd_list_original,
#'                pthresh = 0.05, show_first = 25,
#'                show_gene_names = FALSE, dist = "euclidean",
#'                clust = "average", show_na = FALSE
#' )

# deedee_heatmap(dde = dde,
#                pthresh = 0.05, show_first = 25,
#                show_gene_names = FALSE, dist = "euclidean",
#                clust = "average", show_na = FALSE
# )
#
# deedee_heatmap(dde = dde,
#                pthresh = 0.05, show_first = 25,
#                show_gene_names = TRUE, dist = "euclidean",
#                clust = "average", show_na = TRUE
# )
# TODOs: better scale, blue to red
# TODOs: do this by overriding the selection of genes
# TODOs: add option to return the data frame instead of the plotted object


