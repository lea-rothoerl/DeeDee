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
#'                 (possible values = "pval1" (default),
#'                 "pval2")
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
                      color_by = "pval1") {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  checkmate::assert_number(select1, lower = 1, upper = length(data))
  checkmate::assert_number(select2, lower = 1, upper = length(data))
  choices <- c("pval1", "pval2")
  checkmate::assert_choice(color_by, choices)

  # ---------------------------- data preparation -----------------------------
  data_red <- list(data[[select1]], data[[select2]])

  for(i in 1:length(data_red)) {
    colnames(data_red[i][[1]]) <- c(paste("logFC", i, sep=""),
                                    paste("pval", i, sep=""))
  }

  # ------------------- creation of the resulting qq plot ---------------------
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
    sx <- approx(1L:lenx, sx, n = leny)$y
    pval1 <- approx(1L:lenx, pval1, n = leny)$y
  }
  if (leny > lenx) {
    sy <- approx(1L:leny, sy, n = lenx)$y
    pval2 <- approx(1L:leny, pval2, n = lenx)$y
  }

  qq <- data.frame(x = sx, y = sy, pval1 = pval1, pval2 = pval2)

  qq_f <- qq[qq$pval1 <= pthresh | qq$pval2 <= pthresh, ]

  if (length(qq_f[[1]]) == 0 || length(qq_f[[2]]) == 0) {
    return(NULL)
  }

  res <- ggplot2::ggplot(qq_f, ggplot2::aes(x, y, col = -log10(get(color_by)))) +
    ggplot2::geom_point() +
    viridis::scale_color_viridis(option = "magma") +
    ggplot2::xlab(names(data)[select1]) +
    ggplot2::ylab(names(data)[select2]) +
    ggplot2::labs(color= paste("-log10(", color_by, ")", sep="")) +
    ggplot2::theme_light()

    # --------------------------------- return ----------------------------------
  return(res)
}
