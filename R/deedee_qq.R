#' DeeDee QQ Plot with Multiple Contrasts
#'
#' @description `deedee_qq` creates a plot containing Q-Q-lines comparing
#' the statistical distribution of the logFC of the genes in each input datasets
#' to a contrast chosen from the input data.
#'
#' @param data instance of the DeeDeeExperiment class
#' @param ref index of the contrast in data to be used as reference contrast
#'            (default = 1)
#' @param pthresh threshold for adjusted p-values to be in-/excluded
#'                (default = 0.05)
#' @param ref_line boolean, draws a line on the angle bisector as reference for
#'                 optimal agreement (default = TRUE)
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' data(dde_macrophage, package = "DeeDee")
#' deedee_qq(dde_macrophage)
#'
#' @export
#'

deedee_qq <- function(data,
                          ref = 1,
                          pthresh = 0.05,
                          ref_line = TRUE) {

  # ----------------------------- argument check ------------------------------
  checkmate::assertClass(data, "DeeDeeExperiment")
  data_list <- deedee_from_dde(data)
  contrasts <- names(data_list)
  checkmate::assert(length(contrasts) >= 2)
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  checkmate::assert_number(ref, lower = 1, upper = length(contrasts))
  checkmate::assert_logical(ref_line)

  # ------------------- creation of the resulting qq plot ---------------------
  ref <- as.integer(ref)

  output <- list()

  ref_data <- data_list[[ref]]
  ref_data <- ref_data[!is.na(ref_data$logFC), , drop = FALSE]

  for (i in seq_along(contrasts)) {

    if (i == ref) {
      next
    }

    other_data <- data_list[[i]]
    other_data <- other_data[!is.na(other_data$logFC), , drop = FALSE]

    # Sort logFC values independently, as in deedee_qq()
    x <- sort(ref_data$logFC)
    y <- sort(other_data$logFC)

    padj1 <- ref_data$padj[order(ref_data$logFC)]
    padj2 <- other_data$padj[order(other_data$logFC)]

    lenx <- length(x)
    leny <- length(y)

    # Handle contrasts with different numbers of genes
    if (lenx == 0 || leny == 0) {
      next
    }

    if (leny < lenx) {
      padj1 <- stats::approx(seq_along(padj1), padj1,
                             n = leny)$y
      x <- stats::approx(seq_along(x), x, n = leny)$y
    }

    if (leny > lenx) {
      padj2 <- stats::approx(seq_along(padj2), padj2,
                             n = lenx)$y
      y <- stats::approx(seq_along(y), y, n = lenx)$y
    }

    # Retain points significant in either contrast
    qq <- data.frame(
      x = x,
      y = y,
      padj1 = padj1,
      padj2 = padj2
    )

    qq <- qq[!is.na(qq$padj1) & !is.na(qq$padj2) &
               (qq$padj1 <= pthresh | qq$padj2 <= pthresh), ,
             drop = FALSE]

    if (nrow(qq) == 0) {
      next
    }

    output[[contrasts[i]]] <- qq[, c("x", "y"), drop = FALSE]
  }

  if (length(output) == 0) {
    stop("No Q-Q lines to plot: no contrasts contain points passing pthresh.",
         call. = FALSE)
  }

  qq_data <- dplyr::bind_rows(output, .id = "contrast")

  axis_range <- range(c(qq_data$x, qq_data$y), finite = TRUE)
  axis_pad <- max(diff(axis_range) * 0.05, 0.05)
  axis_range <- axis_range + c(-axis_pad, axis_pad)

  res <- ggplot2::ggplot(
    qq_data,
    ggplot2::aes_string("x", "y", colour = "contrast")
  ) +
    ggplot2::xlab(paste0("reference ", contrasts[ref], " (logFC)")) +
    ggplot2::ylab("logFC (other contrasts, see legend)") +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(limits = axis_range) +
    ggplot2::scale_y_continuous(limits = axis_range) +
    ggplot2::theme_light() +
    viridis::scale_color_viridis(
      option = "magma", discrete = TRUE,
      begin = 0, end = 0.9
  )
  if (ref_line) {
     res <- res +
      ggplot2::geom_abline(intercept = 0,
                           slope = 1,
                           linetype = "dashed",
                           color = "grey")
  }


  # --------------------------------- return ----------------------------------
  return(res)
}
