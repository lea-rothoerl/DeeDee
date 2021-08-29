#' DeeDee Concordance At the Top Plot
#'
#' @description `deedee_cat` creates a Concordance at the Top plot depicting
#' the concordance of genes in the first n elements of the logFC-ranked lists
#' for the given (non-reference) input contrasts against a chosen reference
#' contrast from the input list.
#'
#' @param data named list of results from deedee_prepare()
#' @param ref index of the contrast in data to be used as reference contrast
#'            (default = 1)
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05
#' @param mode sort by highest logFC (`up`, default), lowest logFC (`down`) or
#'             greatest deviation from zero (`both`)
#' @param maxrank highest rank that should be displayed (default = 1000)
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
#' deedee_cat(dd_list, ref = 1, maxrank = 1000, mode = "up", pthresh = 0.05)
#' @export
#'

deedee_cat <- function(data,
                       ref = 1,
                       maxrank = 1000,
                       mode = "up",
                       pthresh = 0.05) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(ref, lower = 1, upper = length(data))
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  checkmate::assert_number(maxrank, lower = 1)
  choices <- c("up", "down", "both")
  checkmate::assert_choice(mode, choices)

  # ---------------------------- data preparation -----------------------------
  for (i in 1:length(data)) {
    data[i][[1]] <- subset(
      data[i][[1]],
      data[i][[1]]$pval < pthresh
    )

    if (length(data[i][[1]][[1]]) == 0) {
      return(NULL)
    }

    data[i][[1]] <- data[i][[1]]["logFC"] # remove p-value column
    colnames(data[i][[1]]) <- c(paste("logFC", i, sep = ""))
    data[i][[1]] <- as.matrix(data[i][[1]]) # conversion to matrix
    names <- rownames(data[i][[1]])
    data[i][[1]] <- as.vector(data[i][[1]]) # conversion to vector
    names(data[i][[1]]) <- names
    if (mode == "up") {
      data[i][[1]] <- sort(data[i][[1]], decreasing = TRUE)
    } else if (mode == "down") {
      data[i][[1]] <- sort(data[i][[1]], decreasing = FALSE)
    } else if (mode == "both") {
      data[i][[1]] <- abs(data[i][[1]])
      data[i][[1]] <- sort(data[i][[1]], decreasing = TRUE)
    }

    data[i][[1]] <- names(data[i][[1]])
  }

  # ----------------------- calculation of concordance ------------------------
  output <- list()
  nm <- c()

  for (i in 1:length(data)) {
    if (i != ref) {
      output[[i]] <- data.frame(
        rank = 1:min(maxrank, length(data[i][[1]])),
        concordance = NA
      )

      for (j in 1:nrow(output[[i]])) {
        intsec <- intersect(data[ref][[1]][1:j], data[i][[1]][1:j])
        output[[i]][[j, "concordance"]] <- length(intsec) / j
      }
      nm[[i]] <- names(data[i])
    }
  }

  names(output) <- nm

  # ------------------- creation of the resulting CAT plot --------------------
  res <- ggplot2::ggplot(
    dplyr::bind_rows(output, .id = "contrast"),
    ggplot2::aes_string("rank", "concordance", colour = "contrast")
  ) +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    viridis::scale_color_viridis(
      option = "magma", discrete = TRUE,
      begin = 0, end = 0.9
    ) +
    ggplot2::annotate("text",
      label = paste("reference: ",
        names(data)[ref],
        sep = ""
      ),
      x = maxrank * 0.8,
      y = max(output[[i]][["concordance"]]) * 1.1
    )

  # --------------------------------- return ----------------------------------
  return(res)
}
