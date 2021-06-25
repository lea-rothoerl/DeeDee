#' Title
#'
#' @description "deedee_cat" creates a Concordance at the Top plot depicting
#' the concordance of genes in the first n elements of the logFC-ranked lists
#' for the given input datasets.
#'
#' @param data (named) list of results from deedee_prepare
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05
#' @param maxrank highest rank that should be displayed (default = 1000)
#'
#' @return ggplot object (plottable with show()/print())
#'
#' @examples
#'
#' @export
#'

deedee_cat <- function(data,
                       pthresh = 0.05,
                       maxrank = 1000) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(data, type = "data.frame", min.len = 2)
  for (i in 1:length(data)) {
    checkmate::assert_data_frame(data[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0)
  checkmate::assert_number(maxrank, lower = 1)

  # ---------------------------- data preparation -----------------------------
  for(i in 1:length(data)) {
    data[i][[1]] <- subset(data[i][[1]],
                               data[i][[1]]$pval < pthresh)
    data[i][[1]] <- data[i][[1]]["logFC"]   # removing p-value column
    colnames(data[i][[1]]) <- c(paste("logFC", i, sep=""))
    data[i][[1]] <- as.matrix(data[i][[1]]) # conversion to matrix
    names <- rownames(data[i][[1]])
    data[i][[1]] <- as.vector(data[i][[1]]) # conversion to vector
    names(data[i][[1]]) <- names
    data[i][[1]] <- sort(data[i][[1]], decreasing = TRUE)
    data[i][[1]] <- names(data[i][[1]])
  }

  # ----------------------- calculation of concordance ------------------------
  output <- data.frame(rank=1:min(maxrank, min(length(data[i][[1]]))),
                       concordance=NA)

  for (i in 1:nrow(output)){
    intsec <- intersect(data[1][[1]][1:i], data[2][[1]][1:i])
    dat_new <- rlist::list.append(data[1][[1]][1:i], data[2][[1]][1:i])
    if (length(data) > 2) {
      for (j in 3:length(data)) {
        intsec <- rlist::list.append(intersect(dat_new, data[j][[1]][1:i]))
        dat_new <- rlist::list.append(dat_new, data[j][[1]][1:i])
      }
    }
    output[i,"concordance"] <- length(intsec)/i
  }

  # ------------------------- calculation of the AUC --------------------------
  auc <- DescTools::AUC(output$rank, output$concordance)
  auc <- auc/(min(maxrank, min(length(data[i][[1]]))))
  # auc <- mean(output$concordance)

  # ------------------- creation of the resulting CAT plot --------------------
  res <- ggplot2::ggplot(data = output, ggplot2::aes(rank, concordance)) +
      ggplot2::geom_line()
  res <- res + ggplot2::annotate("text",
                        x = (7/8*length(output$rank)),
                        y = 0.9,
                        label = paste("AUC =", round(auc, 2), sep = " "))

  # --------------------------------- return ----------------------------------
  # print(paste("Area under curve: ", round(auc, 2), sep=""))
  return(res)
}
