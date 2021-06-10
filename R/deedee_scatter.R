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
#'                 (possible values = "pval1" (default), "pval2", "pval_mean",
#'                 "idr") (blue -> lower values, red -> higher values)
#' @param ggplot output as ggplot (default = FALSE)
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
                           color_by = "pval1",
                           ggplot = FALSE) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0)
  checkmate::assert_number(select1, lower = 1, upper = length(data))
  checkmate::assert_number(select2, lower = 1, upper = length(data))
  choices <- c("pval1", "pval2", "pval_mean", "idr")
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
  comp <- comp[stats::complete.cases(comp[colnames(comp)]), ]  # remove rows with NAs

  if (color_by == "pval_mean") {
    comp$col <- (comp$pval.x + comp$pval.y)/2  # mean of p-values
  }
    else if (color_by == "pval1") {
      comp$col <- comp$pval.x
    }
    else if (color_by == "pval2") {
      comp$col <- comp$pval.y
    }
    # else if (color_by == "idr") {
      # source("~/Development/DeeDee_wip/deedee_idr.R")
      # idr <- deedee_idr(data = data, pthresh = pthresh)
      # idr <- rownames_to_column(idr)
      # colnames(idr) <- c("rowname", "idr")
      # comp <- full_join(comp, idr, by = "rowname", copy = TRUE)
      # names(comp)[names(comp) == 'idr'] <- 'col'
    # }

  comp <- subset(comp, select = -c(comp$pval.y, comp$pval.x)) # removing single p-values
  # comp <- tibble::column_to_rownames(comp, "rowname")

  # -------------------------------- coloring ---------------------------------
  comp$col <- viridis::viridis(1000, option = "magma")[as.numeric(cut(comp$col,
                                                             breaks = 1000))]

  # ----------------- creation of the resulting scatter plot ------------------
  if (ggplot == FALSE) {
    res <- ggplotify::as.ggplot(function() (plot(comp$logFC.x,
              comp$logFC.y,
              xlab = names(data)[select1],
              ylab = names(data)[select2],
              col = comp$col,
              pch = 20,
              xlim=c(min(comp$logFC.x), max(comp$logFC.x)),
              ylim=c(min(comp$logFC.y), max(comp$logFC.y)))))
  }

  else {
    res <- ggplot2::ggplot(data = comp, aes(logFC.x, logFC.y)) +
      geom_point(colour = comp$col)+
      xlab(names(data)[select1]) +
      ylab(names(data)[select2])
  }

  # --------------------------------- return ----------------------------------
  return(res)
}
