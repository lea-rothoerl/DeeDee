#' DeeDee logFC Scatter Plot
#'
#' @description `deedee_scatter` creates a scatterplot of the genes in two input
#' datasets based on their logFC values.
#'
#' @param data instance of the DeeDeeExperiment class
#' @param select1 index of first data-list element to be used (default = 1)
#' @param select2 index of second data-list element to be used (default = 2)
#' @param color_by indicates which set of values the output should be colored by
#'                 (possible values = `pval1` (default), `pval2`)
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#' @param show_symbols boolean whether to show gene symbols instead of Ensembl
#'.                    IDs, `TRUE` by default
#' @param species species the data comes from (for clickable IDs),
#'                `Homo_sapiens` by default
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' data(dde_macrophage, package = "DeeDee")
#' deedee_scatter(dde_macrophage, select1 = 1, select2 = 2)
#'
#' @export
#'

deedee_scatter <- function(data,
                           select1 = 1,
                           select2 = 2,
                           color_by = "pval1",
                           pthresh = 0.05,
                           show_symbols = TRUE,
                           species = "Homo_sapiens") {

  # ----------------------------- argument check ------------------------------
  checkmate::assert(length(names(data@dea)) >= 2)
  data <- deedee_from_dde(data)
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  checkmate::assert_number(select1, lower = 1, upper = length(data))
  checkmate::assert_number(select2, lower = 1, upper = length(data))
  choices <- c("pval1", "pval2")
  checkmate::assert_choice(color_by, choices)
  checkmate::assert_logical(show_symbols)
  checkmate::assert_choice(species,c("Homo_sapiens", "Mus_musculus"))

  # ---------------------------- data preparation -----------------------------
  data_red <- list(data[[select1]], data[[select2]]) # selected samples

  for (i in 1:length(data_red)) {
    data_red[i][[1]] <- subset(
      data_red[i][[1]],
      data_red[i][[1]]$pval < pthresh
    )
    data_red[i][[1]] <- tibble::rownames_to_column(data_red[i][[1]])
  }
  comp <- dplyr::inner_join(data_red[1][[1]],
    data_red[2][[1]],
    by = "rowname",
    copy = FALSE
  )
  comp <- comp[stats::complete.cases(comp[colnames(comp)]), ]

  names(comp) <- c("rowname",
                   "logFC1", "pval1", "symbol1",
                   "logFC2", "pval2", "symbol2")

  if (length(comp[, 1]) == 0) {
    return(NULL)
  }

  if (show_symbols) {
    comp$label <- mosdef::create_link_NCBI(comp$symbol1)
  } else {
    comp$label <- mosdef::create_link_ENSEMBL(
      comp$rowname,
      species = species)
  }

  # additional hover text
  contrast1_name <- names(data)[select1]
  contrast2_name <- names(data)[select2]

  comp$hover_text <- paste0(
    "Symbol: ", comp$symbol1, "<br>",
    "Ensembl ID: ", comp$rowname, "<br>",
    contrast1_name, ": logFC = ", round(comp$logFC1, 2),
    ", p = ", formatC(comp$pval1, format = "e", digits = 2), "<br>",
    contrast2_name, ": logFC = ", round(comp$logFC2, 2),
    ", p = ", formatC(comp$pval2, format = "e", digits = 2), "<br>",
    comp$label
  )

  axis_range <- range(c(comp$logFC1, comp$logFC2), finite = TRUE)
  axis_pad <- max(diff(axis_range) * 0.05, 0.05)
  axis_range <- axis_range + c(-axis_pad, axis_pad)

  # ----------------- creation of the resulting scatter plot ------------------
  res <- ggplot2::ggplot(data = comp, ggplot2::aes(logFC1, logFC2,
    col = get(color_by)
  )) +
    ggplot2::geom_point(ggplot2::aes(text = label)) +
    viridis::scale_color_viridis(option = "magma") +
    ggplot2::xlab(names(data)[select1]) +
    ggplot2::ylab(names(data)[select2]) +
    ggplot2::scale_x_continuous(limits = axis_range) +
    ggplot2::scale_y_continuous(limits = axis_range) +
    ggplot2::labs(color = color_by) +
    ggplot2::theme_light()


  res <- plotly::ggplotly(res, tooltip = "text")

  # --------------------------------- return ----------------------------------
  return(res)
}
