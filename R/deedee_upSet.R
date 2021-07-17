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
                         mode = "both") {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0)
  choices <- c("up", "down", "both", "both_colored")
  checkmate::assert_choice(mode, choices)

  # ---------------------------- data preparation -----------------------------
  for(i in 1:length(data)){
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
    if (mode != "both_colored") {
      data[i][[1]] <- subset(data[i][[1]], select = "rowname")
      colnames(data[i][[1]]) <- names(data)[i]
      data[i] <- data[i][[1]]
    }
  }
  # -------------------------------- coloring ---------------------------------
  if (mode == "both_colored") {
    comp <- dplyr::full_join(data[1][[1]], data[2][[1]], by = "rowname",
                             copy = FALSE)
    if (length(data) > 2) {
      for (i in 3:length(data)) {
        comp <- dplyr::full_join(comp, data[i][[1]], by = "rowname", copy = FALSE)
      }
    }

    comp <- tibble::column_to_rownames(comp, "rowname")
    count <- vector(mode = "numeric", length = length(comp["logFC.x"]))
    comp <- cbind(comp, count)
    comp["count"] <- rowSums(comp, na.rm = TRUE)

    na_count <- vector(mode = "numeric", length = length(comp[, "logFC.x"]))
    comp <- cbind(comp, na_count)
    comp["na_count"] <- rowSums(is.na(comp))

    all_up <- vector(mode = "logical", length = length(comp[, "logFC.x"]))
    comp <- cbind(comp, all_up)
    all_down <- vector(mode = "logical", length = length(comp[, "logFC.x"]))
    comp <- cbind(comp, all_down)

    for (i in 1:length(comp[, "logFC.x"])) {
      if ((comp[i, "count"] + comp[i, "na_count"]) == 4) {
        comp[i, "all_up"] <- TRUE
      }
      else if ((comp[i, "count"] - comp[i, "na_count"]) == -4) {
        comp[i, "all_down"] <- TRUE
      }
    }

    comp <- subset(comp, select = -c(count, na_count))
    comp[is.na(comp)] <- 0
    comp[comp == -1] <- 1
    colnames(comp) <- c(names(data), "all_up", "all_down")

    same_dir <- vector(mode = "logical", length = length(comp[, "all_up"]))
    comp <- cbind(comp, same_dir)

    for (i in 1:length(comp[, "all_up"])) {
      if (comp[i, "all_up"] == TRUE || comp[i, "all_down"] == TRUE) {
        comp[i, "same_dir"] <- TRUE
      }
    }

    for (i in 1:(length(data))) {
      data[i][[1]] <- subset(data[i][[1]], select = "rowname")
      colnames(data[i][[1]]) <- names(data)[i]
      data[i] <- data[i][[1]]
    }
  }

  # ------------------ creation of the resulting UpSet plot -------------------
  col_up = viridis::viridis(n = 1, begin = 0.4, option="magma")
  col_down = viridis::viridis(n = 1, begin = 0.9, option="magma")

  if (mode == "both_colored") {
    res <- UpSetR::upset(comp,
                 order.by = "freq",
                 nsets = length(data),
                 queries = list(list(query = UpSetR::elements,
                                     params = c("same_dir", TRUE),
                                     color = col_up, active = TRUE,
                                     query.name = "only positive logFCs"),
                                list(query = UpSetR::elements,
                                     params = list("all_down", TRUE),
                                     color = col_down, active = TRUE,
                                     query.name = "only negative logFCs")),
                 query.legend = "top")

  }
  else {
    res <- UpSetR::upset(UpSetR::fromList(data),
                         order.by = "freq",
                         nsets = length(data))
  }

  # --------------------------------- return ----------------------------------
  return(res)
}
