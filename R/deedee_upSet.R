#' DeeDee UpSet Plot
#'
#' @description "deedee_upSet" creates an UpSet plot depicting the overlaps of
#' differentially expressed genes in the input datasets.
#'
#' @param data list of results from deedee_prepare
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#' @param mode show all overlapping DE genes ("both", default),
#'             all overlapping genes colored by DE direction ("both_colored",
#'             red -> all up, blue -> all down),
#'             only conjointly up-regulated ("up")
#'             or only conjointly down-regulated ("down") genes
#'
#' @return upset element (plottable with show()/print())
#'
#' @examples
#'
#' @export
#'

deedee_upSet <- function(data,
                         pthresh = 0.05,
                         mode = "both",
                         min_setsize = 10) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  choices <- c("up", "down", "both", "both_colored")
  checkmate::assert_choice(mode, choices)
  checkmate::assert_number(min_setsize, lower = 0)

  # ---------------------------- data preparation -----------------------------
  for(i in 1:length(data)) {
    data[i][[1]] <- subset(data[i][[1]], data[i][[1]]$pval < pthresh) # pthresh

    if (length(data[i][[1]][[1]]) == 0) {
      return(NULL)
    }

    data[i][[1]] <- data[i][[1]]["logFC"]   # removing p-value column
    data[i][[1]]$logFC[data[i][[1]]$logFC < 0] <- -1  # 1/0/-1 distinction
    data[i][[1]]$logFC[data[i][[1]]$logFC > 0] <- 1

    if (mode == "both" || mode == "both_colored") {
      data[i][[1]] <- subset(data[i][[1]], data[i][[1]]$logFC != 0)
    }

    if (mode == "up") {
      data[i][[1]] <- subset(data[i][[1]], data[i][[1]]$logFC == 1)
    }

    if (mode == "down") {
      data[i][[1]] <- subset(data[i][[1]], data[i][[1]]$logFC == -1)
    }
    data[i][[1]] <- tibble::rownames_to_column(data[i][[1]])
    names(data[i][[1]]) <- c("rowname", names(data[i]))
  }

  comp <- dplyr::full_join(data[1][[1]], data[2][[1]], by = "rowname",
                             copy = FALSE)

  if (length(data) > 2) {
    for (i in 3:length(data)) {
      comp <- dplyr::full_join(comp, data[i][[1]], by = "rowname", copy = FALSE)
    }
  }
  comp <- tibble::column_to_rownames(comp)
  contrasts <- colnames(comp)

  # -------------------------------- coloring ---------------------------------
  if (mode == "both_colored") {
    count <- vector(mode = "numeric", length = length(comp[[1]]))
    comp <- cbind(comp, count)
    comp["count"] <- rowSums(comp, na.rm = TRUE)

    na_count <- vector(mode = "numeric", length = length(comp[[1]]))
    comp <- cbind(comp, na_count)
    comp["na_count"] <- rowSums(is.na(comp))

    dr <- vector(mode = "character", length = length(comp[[1]]))
    comp <- cbind(comp, dr)

    nsets <- length(data)

    for (i in 1:length(comp[[1]])) {
      if ((comp[i, "count"] + comp[i, "na_count"]) == nsets) {
        comp[i, "dr"] <- "up"
      }
      else if ((comp[i, "count"] - comp[i, "na_count"]) == -nsets) {
        comp[i, "dr"] <- "down"
      }
      else {
        comp[i, "dr"] <- "different"
      }
    }

    comp <- subset(comp, select = -c(count, na_count))
    comp[is.na(comp)] <- FALSE
    comp[comp == -1 | comp == 1] <- TRUE
  }

  # ------------------ creation of the resulting UpSet plot -------------------
  if (mode == "both_colored") {
    col_up = viridis::viridis(n = 1, begin = 0.4, option="magma")
    col_down = viridis::viridis(n = 1, begin = 0.9, option="magma")
    res <- ComplexUpset::upset(comp,
            contrasts,
            guide = NULL,
            base_annotations = list(
              "Intersection of genes" = ComplexUpset::intersection_size(
                counts = FALSE,
                mapping = ggplot2::aes(fill = dr))
              + ggplot2::scale_fill_manual(values = c("different"= "dark grey",
                                                      "down" = col_down,
                                                      "up" = col_up),
                                           name = "Regulation direction")
                               ),
            width_ratio = 0.1,
            min_size = min_setsize)
  }
  else {
    comp <- (comp != 0 & !is.na(comp))
    comp <- as.data.frame(comp)
    res <- ComplexUpset::upset(comp,
                               contrasts,
                               guide = NULL)
  }

  # --------------------------------- return ----------------------------------
  return(res)
}
