
#' deedee_cat(dd_list_original, ref = 1, maxrank = 1000, mode = "up", pthresh = 0.05)

#' Title
#'
#' @param dde todo
#' @param ref todo
#' @param maxrank todo
#' @param mode todo
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
#'   del
#' )
#' dde
ddedde_cat <- function(dde,
                       ref = 1,
                       maxrank = 1000,
                       mode = "up",
                       pthresh = 0.05) {

  # checkmate::assert_list(data, type = "data.frame", min.len = 2)
  # for (i in 1:length(data)) {
  #   checkmate::assert_data_frame(data[[i]], type = "numeric")
  # }
  # checkmate::assert_number(ref, lower = 1, upper = length(data))
  # checkmate::assert_number(pthresh, lower = 0, upper = 1)
  # checkmate::assert_number(maxrank, lower = 1)
  # choices <- c("up", "down", "both")
  # checkmate::assert_choice(mode, choices)

  dea_list <- get_deas_list(dde)

  for (i in 1:length(dea_list)) {
    dea_list[i][[1]] <- subset(
      dea_list[i][[1]],
      dea_list[i][[1]]$padj < pthresh
    )

    if (length(dea_list[i][[1]][[1]]) == 0) {
      return(NULL)
    }

    dea_list[i][[1]] <- dea_list[i][[1]]["log2FoldChange"] # remove p-value column
    colnames(dea_list[i][[1]]) <- c(paste("log2FoldChange", i, sep = ""))
    dea_list[i][[1]] <- as.matrix(dea_list[i][[1]]) # conversion to matrix
    names <- rownames(dea_list[i][[1]])
    dea_list[i][[1]] <- as.vector(dea_list[i][[1]]) # conversion to vector
    names(dea_list[i][[1]]) <- names
    if (mode == "up") {
      dea_list[i][[1]] <- sort(dea_list[i][[1]], decreasing = TRUE)
    } else if (mode == "down") {
      dea_list[i][[1]] <- sort(dea_list[i][[1]], decreasing = FALSE)
    } else if (mode == "both") {
      dea_list[i][[1]] <- abs(dea_list[i][[1]])
      dea_list[i][[1]] <- sort(dea_list[i][[1]], decreasing = TRUE)
    }

    dea_list[i][[1]] <- names(dea_list[i][[1]])
  }

  output <- list()
  nm <- c()

  for (i in 1:length(dea_list)) {
    if (i != ref) {
      output[[i]] <- data.frame(
        rank = 1:min(maxrank, length(dea_list[i][[1]])),
        concordance = NA
      )

      for (j in 1:nrow(output[[i]])) {
        intsec <- intersect(dea_list[ref][[1]][1:j], dea_list[i][[1]][1:j])
        output[[i]][[j, "concordance"]] <- length(intsec) / j
      }
      nm[[i]] <- names(dea_list[i])
    }
  }

  names(output) <- nm

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
                                    names(dea_list)[ref],
                                    sep = ""
                      ),
                      x = maxrank * 0.8,
                      y = max(output[[i]][["concordance"]]) * 1.1
    )

  return(res)
}

# ddedde_cat(dde,
#            ref = 1,
#            maxrank = 1000,
#            mode = "up",
#            pthresh = 0.05)



