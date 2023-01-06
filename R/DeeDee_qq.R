#' Title
#'
#' @param dde todo
#' @param select1 todo
#' @param select2 todo
#' @param color_by todo
#' @param as_line todo
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
#'
#' deedee_qq(dde)
#'
deedee_qq <- function(dde,
                      select1 = 1,
                      select2 = 2,
                      color_by = "pval1",
                      as_line = FALSE,
                      pthresh = 0.05) {

  # checkmate::assert_list(data, type = "data.frame", min.len = 2)
  # for (i in 1:length(data)) {
  #   checkmate::assert_data_frame(data[[i]], type = "numeric")
  # }
  # checkmate::assert_number(pthresh, lower = 0, upper = 1)
  # checkmate::assert_number(select1, lower = 1, upper = length(data))
  # checkmate::assert_number(select2, lower = 1, upper = length(data))
  # choices <- c("pval1", "pval2")
  # checkmate::assert_choice(color_by, choices)


  dea_list <- get_dea_list(dde)

  data_red <- list(dea_list[[select1]], dea_list[[select2]])

  for (i in seq_len(length(data_red))) {
    colnames(data_red[[i]]) <- c(
      paste("log2FoldChange", i, sep = ""),
      paste("pval", i, sep = ""),
      paste("padj", i, sep = "")
    )
  }

  x <- data_red[[1]]$log2FoldChange
  y <- data_red[[2]]$log2FoldChange
  pval1 <- data_red[[1]]$padj
  pval2 <- data_red[[2]]$padj
  names(x) <- row.names(data_red[[1]])
  names(y) <- row.names(data_red[[2]])

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

  if (as_line == FALSE) {
    res <- ggplot2::ggplot(qq_f, ggplot2::aes(x, y, col = -log10(get(color_by)))) +
      ggplot2::geom_point() +
      viridis::scale_color_viridis(option = "magma") +
      ggplot2::xlab(names(dea(dde))[select1]) +
      ggplot2::ylab(names(dea(dde))[select2]) +
      ggplot2::labs(color = paste("-log10(", color_by, ")", sep = "")) +
      ggplot2::theme_light()
  } else {
    res <- ggplot2::ggplot(qq_f, ggplot2::aes(x, y, col = -log10(get(color_by)))) +
      ggplot2::geom_line() +
      viridis::scale_color_viridis(option = "magma") +
      ggplot2::xlab(names(data)[select1]) +
      ggplot2::ylab(names(data)[select2]) +
      ggplot2::labs(color = paste("-log10(", color_by, ")", sep = "")) +
      ggplot2::theme_light()
  }

  return(res)
}

# deedee_qq(dde,
#           pthresh = 0.05,
#           select = 1,
#           select2 = 2,
#           color_by = "pval1")

# compare the outputs




