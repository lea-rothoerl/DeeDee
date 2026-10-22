#' DeeDee Overview Plot
#'
#' @description `deedee_overview_plot` creates a plot depicting the number of
#'              differentially expressed genes in each contrast of a DDE. The
#'              shares of up- vs. down-regulated genes are plotted as
#'              color-coded bars for each contrast in a collective ggplot.
#'
#' @param data instance of the DeeDeeExperiment class
#' @param pthresh threshold for adjusted p-values to be in-/excluded
#'                (default = 0.05)
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' data(dde_macrophage, package = "DeeDee")
#' deedee_overview_plot(dde_macrophage)
#'
#' @export
#'

deedee_overview_plot <- function(data, pthresh = 0.05) {
  overview <- deedee_overview(data, pthresh)

  long <- data.frame(
    contrast = rep(overview$contrast, times = 2),
    direction = rep(c("up", "down"), each = nrow(overview)),
    n = c(overview$up, overview$down),
    stringsAsFactors = FALSE)

  long$direction <- factor(long$direction, levels = c("down", "up"))
  long$contrast <- factor(long$contrast, levels = rev(overview$contrast))

  col_up <- viridis::viridis(n = 1, begin = 0.4, option = "magma")
  col_down <- viridis::viridis(n = 1, begin = 0.9, option = "magma")

  ggplot2::ggplot(long, ggplot2::aes(x = contrast, y = n, fill = direction )) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::scale_fill_manual(values = c(up = col_up, down = col_down)) +
    ggplot2::labs(x = NULL, y = paste0("DE genes (p < ", pthresh, ")"), fill = NULL) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14))
}
