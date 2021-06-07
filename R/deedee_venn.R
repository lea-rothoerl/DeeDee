#' DeeDee Venn Diagram
#'
#' @description "deedee_venn" creates a Venn diagram depicting the overlaps of
#' differentially expressed genes in the input datasets.
#'
#' @param data (named) list of 2-4 results from deedee_format
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#' @param mode show all overlapping DE genes ("both", default),
#'             only conjointly up-regulated ("up")
#'             or only conjointly down-regulated ("down") genes
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' @export
#'

deedee_venn <- function(data,
                        pthresh = 0.05,
                        mode = "both") {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2, max.len = 4)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0)
  choices <- c("up", "down", "both")
  checkmate::assert_choice(mode, choices)

  # ---------------------------- data preparation -----------------------------
  for(i in 1:length(data)) {
    data[i][[1]] <- subset(data[i][[1]],
                           data[i][[1]]$pval < pthresh)
    if (mode == "up") {
      data[i][[1]] <- subset(data[i][[1]],
                             data[i][[1]]$logFC > 0)
    }
    if (mode == "down") {
      data[i][[1]] <- subset(data[i][[1]],
                             data[i][[1]]$logFC < 0)
    }
    data[i][[1]] <- data[i][[1]]["logFC"]   # removing p-value column
    colnames(data[i][[1]]) <- c(paste("logFC", i, sep=""))
    data[i][[1]] <- as.matrix(data[i][[1]]) # conversion to matrix
    names <- rownames(data[i][[1]])
    data[i][[1]] <- as.vector(data[i][[1]]) # conversion to vector
    names(data[i][[1]]) <- names
    data[i][[1]] <- sort(data[i][[1]], decreasing = TRUE)
    data[i][[1]] <- names(data[i][[1]])
  }

  # ----------------- creation of the resulting venn diagram ------------------
  pal = c(viridis::viridis(length(data), option = "magma"))

  res <- ggvenn::ggvenn(data,
                fill_alpha = 0.2,
                fill_color = pal,
                show_percentage = FALSE,
                stroke_color = "grey80")

  # --------------------------------- return ----------------------------------
  return(res)
}

