
#' deedee_heatmap(dd_list_original,
#'                pthresh = 0.05, show_first = 25,
#'                show_gene_names = FALSE, dist = "euclidean",
#'                clust = "average", show_na = FALSE
#' )


#' Title
#'
#' @param dde todo
#' @param show_first todo
#' @param show_gene_names todo
#' @param dist todo
#' @param clust todo
#' @param show_na todo
#' @param pthresh todo
#'
#' @return todo
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
ddedde_heatmap <- function(dde,
                           show_first = 25,
                           show_gene_names = FALSE,
                           dist = "euclidean",
                           clust = "average",
                           show_na = FALSE,
                           pthresh = 0.05) {

  # checkmate::assert_list(data, type = "data.frame", min.len = 2)
  # for (i in 1:length(data)) {
  #   checkmate::assert_data_frame(data[[i]], type = "numeric")
  # }
  # checkmate::assert_number(pthresh, lower = 0, upper = 1)
  # checkmate::assert_number(show_first, lower = 1)
  # checkmate::assert_logical(show_gene_names)
  # choices1 <- c("euclidean", "manhattan", "pearson", "spearman")
  # checkmate::assert_choice(dist, choices1)
  # choices2 <- c("single", "complete", "average", "centroid")
  # checkmate::assert_choice(clust, choices2)
  # checkmate::assert_logical(show_na)

  dea_list <- get_deas_list(dde)
  # and then proceed as for the old implementation

  for (i in 1:length(dea_list)) {
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

  return(res)

}



#' deedee_heatmap(dd_list_original,
#'                pthresh = 0.05, show_first = 25,
#'                show_gene_names = FALSE, dist = "euclidean",
#'                clust = "average", show_na = FALSE
#' )

# ddedde_heatmap(dde = dde,
#                pthresh = 0.05, show_first = 25,
#                show_gene_names = FALSE, dist = "euclidean",
#                clust = "average", show_na = FALSE
# )
#
# ddedde_heatmap(dde = dde,
#                pthresh = 0.05, show_first = 25,
#                show_gene_names = TRUE, dist = "euclidean",
#                clust = "average", show_na = TRUE
# )
# TODOs: better scale, blue to red
# TODOs: do this by overriding the selection of genes
# TODOs: add option to return the data frame instead of the plotted object


