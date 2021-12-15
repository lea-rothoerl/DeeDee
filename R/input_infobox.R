
input_infobox <- function(deedee_obj,
                          sets) {

  ext <- c()
  filename <- c()
  res <- list()
  type <- c()
  contrast <- c()
  genes <- c()
  count <- 0

  if (!is.null(deedee_obj)) {
    for (j in 1:length(deedee_obj@DeeDeeList)) {
      count <- count + 1
      type[count] <- "DeeDee object"
      filename[count] <- "input as argument"
      contrast[count] <- names(deedee_obj@DeeDeeList)[j]
      genes[count] <- length(deedee_obj@DeeDeeList[j]
                             [[contrast[count]]][["logFC"]])
    }
  }

  # reading out input files
  if (length(sets[, 1] > 0)) {
    for (i in 1:(length(sets[, 1]))) {
      ext[i] <- tools::file_ext(sets[i, "datapath"])

      if (ext[[i]] == "rds" || ext[[i]] == "RDS") {
        res[[i]] <- readRDS(sets[[i, "datapath"]])
      } else if (ext[[i]] == "xlsx") {
        sheets <- readxl::excel_sheets(sets[[i, "datapath"]])
        res[[i]] <- lapply(sheets,
                           readxl::read_excel,
                           path = sets[[i, "datapath"]]
        )
        names(res[[i]]) <- sheets
        for (j in 1:length(sheets)) {
          res[[i]][[sheets[j]]] <- as.data.frame(res[[i]][[sheets[j]]])
          res[[i]][[sheets[j]]] <- tibble::column_to_rownames(
            res[[i]][[sheets[j]]], "rowname"
          )
        }
      } else if (ext[[i]] == "txt") {
        temp <- utils::read.table(sets[[i, "datapath"]])
        res[[i]] <- list(temp)
        names(res[[i]]) <- unlist(strsplit(sets[i, "name"],
                                           split = ".",
                                           fixed = TRUE
        ))[1]
      }

      if (class(res[[i]]) == "DESeqResults" ||
          class(res[[i]]) == "DGEExact" ||
          length(names(res[[i]])) == 6 ||
          class(res[[i]]) == "DeeDeeObject") {
        count <- count + 1
        type[count] <- class(res[[i]])
        filename[count] <- sets[i, "name"]
        contrast[count] <- unlist(strsplit(filename[count],
                                           split = ".",
                                           fixed = TRUE
        ))[1]
        genes[count] <- length(mydata()@DeeDeeList
                               [[contrast[count]]][["logFC"]])
      } else {
        if (class(res[[i]]) == "data.frame") {
          res[[i]] <- list(res[[i]])
          names(res[[i]]) <- unlist(strsplit(sets[i, "name"],
                                             split = ".",
                                             fixed = TRUE
          ))[1]
        }
        for (j in 1:length(res[[i]])) {
          count <- count + 1
          type[count] <- "DeeDee object"
          filename[count] <- sets[i, "name"]
          contrast[count] <- names(res[[i]])[j]
          genes[count] <- length(res[[i]][j]
                                 [[contrast[count]]][["logFC"]])
        }
      }
    }
  }

  df <- data.frame(filename, type, contrast, genes)
  return(df)
}
