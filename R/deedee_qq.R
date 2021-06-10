#' DeeDee QQ Plot
#'
#' @description "deedee_qq" creates a QQ-plot comparing the statistical
#' distribution of the logFC of the genes in two input datasets.
#'
#' @param data (named) list of results from deedee_prepare
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#' @param select1 index of first data-list element to be used (default = 1)
#' @param select2 index of second data-list element to be used (default = 2)
#' @param color_by indicates which set of values the output should be colored by
#'                 (possible values = "logFC1", "logFC2", "pval1" (default),
#'                 "pval2")
#' @param ggplot output as ggplot (default = FALSE)
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' @export
#'

deedee_qq <- function(data,
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
  choices <- c("pval1", "pval2", "logFC1", "logFC2")
  checkmate::assert_choice(color_by, choices)

  # ---------------------------- data preparation -----------------------------
  data_red <- list(data[[select1]], data[[select2]])

  for(i in 1:length(data_red)) {
    data_red[i][[1]] <- subset(data_red[i][[1]], data_red[i][[1]]$pval < pthresh)
    colnames(data_red[i][[1]]) <- c(paste("logFC", i, sep=""),
                                    paste("pval", i, sep=""))
  }

  # -------------------------------- coloring ---------------------------------
  pal <- viridis::viridis(1000, option = "magma")

  if (color_by == "logFC1") {
    data_red$col <- pal[as.numeric(cut(data_red[[1]]$logFC1,
                                            breaks = 1000))]
  }
    else if (color_by == "logFC2") {
      data_red$col <- pal[as.numeric(cut(data_red[[2]]$logFC2,
                                            breaks = 1000))]
   }
    else if (color_by == "pval1") {
      data_red$col <- pal[as.numeric(cut(data_red[[1]]$pval1,
                                            breaks = 1000))]
   }
    else if (color_by == "pval2") {
       data_red$col <- pal[as.numeric(cut(data_red[[2]]$pval2,
                                            breaks = 1000))]
   }

  # ------------------- creation of the resulting qq plot ---------------------
  if (ggplot == FALSE) {
    res <- ggplotify::as.ggplot(function() (stats::qqplot(data_red[1][[1]]$logFC,
                data_red[2][[1]]$logFC,
                plot.it = TRUE,
                xlab = names(data)[select1],
                ylab = names(data)[select2],
                pch = 20,
                xlim = c(min(data_red[1][[1]]$logFC),
                         max(data_red[1][[1]]$logFC)),
                ylim = c(min(data_red[2][[1]]$logFC),
                         max(data_red[2][[1]]$logFC)),
                col = data_red$col)))
  }

  else {
    qq <- as.data.frame(stats::qqplot(data_red[1][[1]]$logFC,
                               data_red[2][[1]]$logFC,
                               plot.it=FALSE,
                               col = data_red$col))
    res <- ggplot2::ggplot(qq, ggplot2::aes(x, y)) +
      ggplot2::geom_point() +
      ggplot2::xlab(names(data)[select1]) +
      ggplot2::ylab(names(data)[select2])
  }

  # --------------------------------- return ----------------------------------
  return(res)
}
