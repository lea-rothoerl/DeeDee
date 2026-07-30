#' DeeDee QQ Plot with Multiple Contrasts
#'
#' @description `deedee_qqmult` creates a plot containing Q-Q-lines comparing
#' the statistical distribution of the logFC of the genes in each input datasets
#' to a contrast chosen from the input data.
#'
#' @param data instance of the DeeDeeExperiment class
#' @param ref index of the contrast in data to be used as reference contrast
#'            (default = 1)
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' TODO
#'
#' @export
#'

deedee_qqmult <- function(data,
                          ref = 1,
                          pthresh = 0.05) {

  # ----------------------------- argument check ------------------------------
  checkmate::assertClass(data, "DeeDeeExperiment")
  data_list <- deedee_from_dde(data)
  contrasts <- names(data_list)
  checkmate::assert(length(contrasts) >= 2)
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  checkmate::assert_number(ref, lower = 1, upper = length(contrasts))

  # ------------------- creation of the resulting qq plot ---------------------
  output <- list()
  nm <- c()

  for (i in seq_along(contrasts)) {
    if (i != ref) {
      output[[i]] <- data.frame(ggplot2::ggplot_build(deedee_qq(
        data = data,
        select1 = ref,
        select2 = i,
        as_line = TRUE
      ))$plot$data)
      nm[[i]] <- contrasts[i]
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
