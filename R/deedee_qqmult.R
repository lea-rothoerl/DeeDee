#' DeeDee QQ Plot with Multiple Contrasts
#'
#' @description TODO
#'
#' @param data TODO
#'
#' @return ggplot object (plottable with show()/print())
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
#' deedee_qqmult(TODO)
#' @export
#'

deedee_qq <- function(data,
                      pthresh = 0.05,
                      ref = 1) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  checkmate::assert_number(ref, lower = 1, upper = length(data))

  # ---------------------------- data preparation -----------------------------
  for (i in length(data)-1) {
    data_red <- list(data[[i]], data[[ref]])

    for (i in 1:length(data_red)) {
      colnames(data_red[i][[1]]) <- c(
        paste("logFC", i, sep = ""),
        paste("pval", i, sep = "")
      )
    }

    x <- data_red[1][[1]]$logFC
    y <- data_red[2][[1]]$logFC
    pval1 <- data_red[1][[1]]$pval
    pval2 <- data_red[2][[1]]$pval
    names(x) <- row.names(data_red[1][[1]])
    names(y) <- row.names(data_red[2][[1]])

    sx_idx <- order(x)
    sy_idx <- order(y)

    sx <- x[sx_idx]
    sy <- y[sy_idx]
    pval1 <- pval1[sx_idx]
    pval2 <- pval2[sy_idx]

    lenx <- length(sx)
    leny <- length(sy)

    if (leny < lenx) {
      sx <- stats::approx(1L:lenx, sx, n = leny)$y
      pval1 <- stats::approx(1L:lenx, pval1, n = leny)$y
    }
    if (leny > lenx) {
      sy <- stats::approx(1L:leny, sy, n = lenx)$y
      pval2 <- stats::approx(1L:leny, pval2, n = lenx)$y
    }

    qq <- data.frame(x = sx, y = sy, pval1 = pval1, pval2 = pval2)

    qq_f <- qq[qq$pval1 <= pthresh | qq$pval2 <= pthresh, ]

    if (length(qq_f[[1]]) == 0 || length(qq_f[[2]]) == 0) {
      return(NULL)
    }
  }

  res <- ggplot2::ggplot(qq_f, ggplot2::aes(x, y, col = -log10(get(color_by)))) +
    ggplot2::geom_point() +
    viridis::scale_color_viridis(option = "magma") +
    ggplot2::xlab(names(data)[select1]) +
    ggplot2::ylab(names(data)[select2]) +
    ggplot2::labs(color = paste("-log10(", color_by, ")", sep = "")) +
    ggplot2::theme_light()

  # --------------------------------- return ----------------------------------
  return(res)
}
