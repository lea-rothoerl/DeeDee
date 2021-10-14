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
    deedee_qq(data = data,
              select1 = ref,
              select2 = i,
              as_line = TRUE)
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
