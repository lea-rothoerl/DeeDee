#' DeeDee logFC Scatterplot
#'
#' @description "deedee_scatter" creates a scatterplot of the genes in two input
#' datasets based on their logFC values.
#'
#' @param data list of results from deedee_prepare
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#' @param select1 index of first data-list element to be used (default = 1)
#' @param select2 index of second data-list element to be used (default = 2)
#' @param color_by indicates which set of values the output should be colored by
#'                 (possible values = "pval1" (default), "pval2", "idr") (blue
#'                 -> lower values, red -> higher values)
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' @export
#'

deedee_scatter <- function(data,
                           pthresh = 0.05,
                           select1 = 1,
                           select2 = 2,
                           color_by = "pval1") {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0)
  checkmate::assert_number(select1, lower = 1, upper = length(data))
  checkmate::assert_number(select2, lower = 1, upper = length(data))
  choices <- c("pval1", "pval2", "idr")
  checkmate::assert_choice(color_by, choices)

  # ---------------------------- data preparation -----------------------------
  data_red <- list(data[[select1]], data[[select2]]) # selected samples

  for(i in 1:length(data_red)){
    data_red[i][[1]] <- subset(data_red[i][[1]],
                               data_red[i][[1]]$pval < pthresh)
    data_red[i][[1]] <- tibble::rownames_to_column(data_red[i][[1]])
  }
  comp <- dplyr::inner_join(data_red[1][[1]],
                     data_red[2][[1]],
                     by = "rowname",
                     copy = FALSE)
  comp <- comp[stats::complete.cases(comp[colnames(comp)]), ]

  names(comp) <- c("rowname", "logFC1", "pval1", "logFC2", "pval2")

  if (length(comp[,1]) == 0) {
    return(NULL)
  }

  # ----------------- creation of the resulting scatter plot ------------------
  res <- ggplot2::ggplot(data = comp, ggplot2::aes(logFC1, logFC2,
                                                   col = get(color_by))) +
      ggplot2::geom_point() +
      viridis::scale_color_viridis(option = "magma") +
      ggplot2::xlab(names(data)[select1]) +
      ggplot2::ylab(names(data)[select2]) +
      ggplot2::labs(color=color_by)

  # --------------------------------- return ----------------------------------
  return(res)
}
