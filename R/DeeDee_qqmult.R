#' Title
#'
#' @param dde todo
#' @param ref todo
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
#' ddedde_qqmult(dde, pthresh = 0.05, ref = 1)
#'
ddedde_qqmult <- function(dde,
                          ref = 1,
                          pthresh = 0.05) {

  # checkmate::assert_list(data, type = "data.frame", min.len = 2)
  # for (i in 1:length(data)) {
  #   checkmate::assert_data_frame(data[[i]], type = "numeric")
  # }
  # checkmate::assert_number(pthresh, lower = 0, upper = 1)
  # checkmate::assert_number(ref, lower = 1, upper = length(data))

  dea_list <- get_dea_list(dde)

  output <- list()
  nm <- c()

  for (i in 1:length(dea_list)) {
    if (i != ref) {
      output[[i]] <- data.frame(
        ggplot2::ggplot_build(ddedde_qq(
          dde = dde,
          select1 = ref,
          select2 = i,
          as_line = TRUE
        ))$plot$data)

      nm[[i]] <- names(dea_list[i])
    }
  }

  names(output) <- nm

  res <- ggplot2::ggplot(
    dplyr::bind_rows(output, .id = "contrast"),
    ggplot2::aes(.data$x, .data$y, colour = .data$contrast)
  ) +
    ggplot2::xlab(names(dea_list)[ref]) +
    ggplot2::ylab("Contrasts") +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    viridis::scale_color_viridis(
      option = "magma", discrete = TRUE,
      begin = 0, end = 0.9
    )

  return(res)
}

#' ddedde_qqmult(dde,
#'               pthresh = 0.05,
#'               ref = 1)
#'
#' deedee_qqmult(dd_list_original, pthresh = 0.05, ref = 2)
#' ddedde_qqmult(dde, pthresh = 0.05, ref = 2)
#'
#' deedee_qqmult(dd_list_original, pthresh = 0.05, ref = 3)
#' ddedde_qqmult(dde, pthresh = 0.05, ref = 3)
#'
#' deedee_qqmult(dd_list_original, pthresh = 0.05, ref = 4)
#' ddedde_qqmult(dde, pthresh = 0.05, ref = 4)




