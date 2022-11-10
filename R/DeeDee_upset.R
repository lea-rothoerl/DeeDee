#' Title
#'
#' @param dde todo
#' @param mode todo
#' @param min_setsize todo
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
#' ddedde_upset(dde, pthresh = 0.05, mode = "both_colored", min_setsize = 10)
ddedde_upset <- function(dde,
                         mode = "both_colored",
                         min_setsize = 10,
                         pthresh = 0.05) {

  # checkmate::assert_list(data, type = "data.frame", min.len = 2)
  # for (i in 1:length(data)) {
  #   checkmate::assert_data_frame(data[[i]], type = "numeric")
  # }
  # checkmate::assert_number(pthresh, lower = 0, upper = 1)
  # choices <- c("up", "down", "both", "both_colored")
  # checkmate::assert_choice(mode, choices)
  # checkmate::assert_number(min_setsize, lower = 0)


  dea_list <- get_dea_list(dde)
  for (i in 1:length(dea_list)) {
    dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$padj < pthresh) # pthresh

    if (length(dea_list[[i]][[1]]) == 0) {
      return(NULL)
    }

    dea_list[[i]] <- dea_list[[i]]["log2FoldChange"] # removing p-value column
    dea_list[[i]]$log2FoldChange[dea_list[[i]]$log2FoldChange < 0] <- -1 # 1/0/-1 distinction
    dea_list[[i]]$log2FoldChange[dea_list[[i]]$log2FoldChange > 0] <- 1

    if (mode == "both" || mode == "both_colored") {
      dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$log2FoldChange != 0)
    }

    if (mode == "up") {
      dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$log2FoldChange == 1)
    }

    if (mode == "down") {
      dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$log2FoldChange == -1)
    }
    dea_list[[i]] <- tibble::rownames_to_column(dea_list[[i]])
    names(dea_list[[i]]) <- c("rowname", names(dea_list[i]))
  }

  comp <- dplyr::full_join(dea_list[[1]], dea_list[[2]],
                           by = "rowname",
                           copy = FALSE
  )

  if (length(dea_list) > 2) {
    for (i in 3:length(dea_list)) {
      comp <- dplyr::full_join(comp, dea_list[[i]], by = "rowname", copy = FALSE)
    }
  }
  comp <- tibble::column_to_rownames(comp)
  contrasts <- colnames(comp)

  if (mode == "both_colored") {
    count <- vector(mode = "numeric", length = length(comp[[1]]))
    comp <- cbind(comp, count)
    comp["count"] <- rowSums(comp, na.rm = TRUE)

    na_count <- vector(mode = "numeric", length = length(comp[[1]]))
    comp <- cbind(comp, na_count)
    comp["na_count"] <- rowSums(is.na(comp))

    dr <- vector(mode = "character", length = length(comp[[1]]))
    comp <- cbind(comp, dr)

    nsets <- length(dea_list)

    for (i in 1:length(comp[[1]])) {
      if ((comp[i, "count"] + comp[i, "na_count"]) == nsets) {
        comp[i, "dr"] <- "up"
      } else if ((comp[i, "count"] - comp[i, "na_count"]) == -nsets) {
        comp[i, "dr"] <- "down"
      } else {
        comp[i, "dr"] <- "different"
      }
    }

    comp <- subset(comp, select = -c(count, na_count))
    comp[is.na(comp)] <- FALSE
    comp[comp == -1 | comp == 1] <- TRUE
  }

  if (mode == "both_colored") {
    col_up <- viridis::viridis(n = 1, begin = 0.4, option = "magma")
    col_down <- viridis::viridis(n = 1, begin = 0.9, option = "magma")
    res <- ComplexUpset::upset(comp,
                               contrasts,
                               guide = NULL,
                               base_annotations = list(
                                 "Intersection of genes" = ComplexUpset::intersection_size(
                                   counts = FALSE,
                                   mapping = ggplot2::aes(fill = dr)
                                 )
                                 + ggplot2::scale_fill_manual(
                                   values = c(
                                     "different" = "dark grey",
                                     "down" = col_down,
                                     "up" = col_up
                                   ),
                                   name = "Regulation direction"
                                 )
                               ),
                               width_ratio = 0.1,
                               min_size = min_setsize
    )
  } else {
    comp <- (comp != 0 & !is.na(comp))
    comp <- as.data.frame(comp)
    res <- ComplexUpset::upset(comp,
                               contrasts,
                               guide = NULL,
                               width_ratio = 0.1,
                               min_size = min_setsize
    )
  }

  return(res)
}


# ddedde_upset(dde,
#              pthresh = 0.05,
#              mode = "both_colored",
#              min_setsize = 10)
