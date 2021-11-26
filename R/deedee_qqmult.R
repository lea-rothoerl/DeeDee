#' DeeDee QQ Plot with Multiple Contrasts
#'
#' @description `deedee_qqmult` creates a plot containing Q-Q-lines comparing
#' the statistical distribution of the logFC of the genes in each input datasets
#' to a contrast chosen from the input data. For a Q-Q plot comparing two
#' contrasts in a p-value-colored manner, have a look at `deedee_qq`.
#'
#' @param data named list of results from deedee_prepare()
#' @param ref index of the contrast in data to be used as reference contrast
#'            (default = 1)
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
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
#' deedee_qqmult(dd_list, pthresh = 0.05, ref = 1)
#' @export
#'

deedee_qqmult <- function(data,
                          ref = 1,
                          pthresh = 0.05) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  checkmate::assert_number(ref, lower = 1, upper = length(data))

  # ------------------- creation of the resulting qq plot ---------------------
  output <- list()
  nm <- c()

  for (i in 1:length(data)) {
    if (i != ref) {
      output[[i]] <- data.frame(ggplot2::ggplot_build(deedee_qq(
        data = data,
        select1 = ref,
        select2 = i,
        as_line = TRUE
      ))$plot$data)
      nm[[i]] <- names(data[i])
    }
  }

  names(output) <- nm

  res <- ggplot2::ggplot(
    dplyr::bind_rows(output, .id = "contrast"),
    ggplot2::aes_string("x", "y", colour = "contrast")
  ) +
    ggplot2::xlab(names(data)[ref]) +
    ggplot2::ylab("see legend") +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    viridis::scale_color_viridis(
      option = "magma", discrete = TRUE,
      begin = 0, end = 0.9
    )


  # --------------------------------- return ----------------------------------
  return(res)
}
