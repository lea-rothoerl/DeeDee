#' DeeDee Venn Diagram
#'
#' @description `deedee_venn` creates a Venn diagram depicting the overlaps of
#' differentially expressed genes in the input datasets.
#'
#' @param data instance of the DeeDeeExperiment class
#' @param mode show all overlapping DE genes (`both`, default),
#'             only conjointly up-regulated (`up`)
#'             or only conjointly down-regulated (`down`) genes
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#' @param select vector of DEA slot indexes to use for plotting, must be 2 to 4
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' data(dde_macrophage, package = "DeeDee")
#' deedee_venn(dde_macrophage, mode = "both", pthresh = 0.05)
#'
#' @export
#'

deedee_venn <- function(data,
                        mode = "both",
                        pthresh = 0.05,
                        select = NULL) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert(length(names(data@dea)) >= 2)
  data <- deedee_from_dde(data)
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  choices <- c("up", "down", "both")
  checkmate::assert_choice(mode, choices)

  n_contrasts <- length(data)

  # --------------------------- selection handling ----------------------------
  if (is.null(select)) {

    if (n_contrasts >= 4) {
      select <- seq_len(4)
    } else {
      select <- seq_len(n_contrasts)
    }

  } else {

    checkmate::assert_integerish(select, lower = 1)
    checkmate::assert_true(length(select) >= 2 && length(select) <= 4)
    checkmate::assert_true(length(unique(select)) == length(select))
    checkmate::assert_true(all(select <= n_contrasts))
  }

  data <- data[select]

  # ---------------------------- data preparation -----------------------------
  for (i in 1:length(data)) {
    data[i][[1]] <- subset(
      data[i][[1]],
      data[i][[1]]$pval < pthresh
    )

    if (length(data[i][[1]][[1]]) == 0) {
      return(NULL)
    }

    if (mode == "up") {
      data[i][[1]] <- subset(
        data[i][[1]],
        data[i][[1]]$logFC > 0
      )
    }
    if (mode == "down") {
      data[i][[1]] <- subset(
        data[i][[1]],
        data[i][[1]]$logFC < 0
      )
    }
    data[i][[1]] <- data[i][[1]]["logFC"] # removing p-value column
    colnames(data[i][[1]]) <- c(paste("logFC", i, sep = ""))
    data[i][[1]] <- as.matrix(data[i][[1]]) # conversion to matrix
    names <- rownames(data[i][[1]])
    data[i][[1]] <- as.vector(data[i][[1]]) # conversion to vector
    names(data[i][[1]]) <- names
    data[i][[1]] <- sort(data[i][[1]], decreasing = TRUE)
    data[i][[1]] <- names(data[i][[1]])
  }

  # ----------------- creation of the resulting venn diagram ------------------
  pal <- c(viridis::viridis(length(data), option = "magma"))

  res <- ggvenn::ggvenn(data,
    fill_alpha = 0.2,
    fill_color = pal,
    show_percentage = FALSE,
    stroke_color = "grey80",
    set_name_size = 4
  )

  # --------------------------------- return ----------------------------------
  return(res)
}
