#' @title DeeDee Heatmap
#'
#' @description "deedee_heatmap" creates a heatmap depicting the logFC as a
#' measure of the differential expression of the first genes in the given
#' datasets.
#'
#' @param data (named) list of 2-5 results from deedee_format
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#' @param show_first indicating the number of genes depicted (default = 500)
#' @param show_gene_names boolean, show rownames next to heatmap
#'                        (default = FALSE)
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' @export
#'

deedee_heatmap <- function(data,
                           pthresh = 0.05,
                           show_first = 25,
                           show_gene_names = FALSE) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0)
  checkmate::assert_number(show_first, lower = 1,
                           upper = max(length(data[[i]]$logFC)))
  checkmate::assert_logical(show_gene_names)

  # ---------------------------- data preparation -----------------------------
  for(i in 1:length(data)){
    data[i][[1]] <- subset(data[i][[1]], data[i][[1]]$pval < pthresh)  # pthresh
    data[i][[1]] <- data[i][[1]]["logFC"]   # removing p-value column
    colnames(data[i][[1]]) <- names(data)[i]  # creating unique colnames
    data[i][[1]] <- tibble::rownames_to_column(data[i][[1]])
  }

  # ----- creation of a comp matrix with the modified results from above ------
  comp <- dplyr::inner_join(data[1][[1]],
                            data[2][[1]],
                            by = "rowname",
                            copy = FALSE)
  for (i in 3:length(data)) {
    comp <- dplyr::inner_join(comp, data[i][[1]], by = "rowname", copy = FALSE)
  }

  row.names(comp) <- comp$rowname
  comp <- subset(comp, select = -c(rowname))  # removing column with rownames
  comp <- comp[stats::complete.cases(comp[colnames(comp)]),]
  comp <- as.matrix(comp)

  # ------------------- creation of the resulting heatmap ---------------------
  if (show_gene_names == FALSE) {
   rownames(comp) <- c()
  }

  col <- viridis::viridis(n = 15, option = "magma")

  res <- ComplexHeatmap::Heatmap(comp[1:min(show_first,
                                                length(comp[,1])),],
                                 col = col)

  res <- ComplexHeatmap::draw(res)

  # res <- ggplotify::as.ggplot(function() gplots::heatmap.2(comp[1:min(show_first,
  #                                                length(comp[,1])),],
  #                                     col = col,
  #                                     margins=c(9,8),
  #                                     key = TRUE,
  #                                     density.info = "density",
  #                                     key.title = NA,
  #                                     key.xlab = NA,
  #                                     key.ylab = NA,
  #                                     trace = "none",
  #                                     distfun = stats::dist))

  # --------------------------------- return ----------------------------------
  return(res)
}
