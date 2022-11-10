#' Title
#'
#' @param dde todo
#' @param select1 todo
#' @param select2 todo
#' @param color_by todo
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
#' ddedde_scatter(dde, 1, 3, color_by = "pval2")
#' ddedde_scatter(dde, 1, 2)
#' ddedde_scatter(dde, 1, 4)
#' ddedde_scatter(dde, 3, 2)
ddedde_scatter <- function(dde,
                           select1 = 1,
                           select2 = 2,
                           color_by = "pval1",
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

  dea_names <- names(dea(dde))
  names1 <- dea_names[select1]
  names2 <- dea_names[select2]

  de_1 <- get_dea_df(dde, names1)
  de_2 <- get_dea_df(dde, names2)
  data_red <- list(
    tibble::rownames_to_column(as.data.frame(de_1)),
    tibble::rownames_to_column(as.data.frame(de_2))) # selected samples

  # 928
  # 2692
  # 185

  comp <- dplyr::inner_join(data_red[[1]],
                            data_red[[2]],
                            by = "rowname",
                            copy = FALSE
  )

  # for (i in 1:length(data_red)) {
  #   data_red[i][[1]] <- subset(
  #     data_red[i][[1]],
  #     data_red[i][[1]]$pval < pthresh
  #   )
  #   data_red[i][[1]] <- tibble::rownames_to_column(data_red[i][[1]])
  # }
  # comp <- dplyr::inner_join(data_red[1][[1]],
  #                           data_red[2][[1]],
  #                           by = "rowname",
  #                           copy = FALSE
  # )
  # comp <- comp[stats::complete.cases(comp[colnames(comp)]), ]

  names(comp) <- c("gene_id",
                   "logFC1", "punadj1", "pval1",
                   "logFC2", "punadj2", "pval2")


  comp_sel <- comp[(comp$pval1 <= pthresh) & (comp$pval2 <= pthresh), ]

  comp_sel$info <- paste0(
    comp_sel$gene_id,": ", comp_sel$logFC1, ", ", comp_sel$logFC2, "\n something else\nanother line"
  )

  if (length(comp[, 1]) == 0) {
    return(NULL)
  }

  res <- ggplot2::ggplot(
    data = comp_sel,
    aes(.data[["logFC1"]], .data[["logFC2"]],
        col = .data[[color_by]],
        label = .data[["gene_id"]],
        info = .data[["info"]]
    )) +
    ggplot2::geom_point() +
    viridis::scale_color_viridis(option = "magma") +
    ggplot2::xlab(names1) +
    ggplot2::ylab(names2) +
    ggplot2::labs(color = color_by) +
    ggplot2::theme_light()

  return(res)
}
