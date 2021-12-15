
rds_input <- function(obj,
                      nm) {

  if (class(obj) == "DeeDeeObject") {
    obj <- obj@DeeDeeList
  }
  if (class(obj) == "DESeqResults") {
    obj <- deedee_prepare(obj, "DESeq2")
    obj <- list(obj)
    names(obj) <- unlist(strsplit(nm,
                                       split = ".",
                                       fixed = TRUE
    ))[1]
  } else if (class(obj) == "DGEExact") {
    obj <- deedee_prepare(obj, "edgeR")
    obj <- list(obj)
    names(obj) <- unlist(strsplit(nm,
                                       split = ".",
                                       fixed = TRUE
    ))[1]
  } else if (class(obj) == "list") {
    for (j in length(obj)) {
      if (checkmate::test_subset(
        names(obj[[j]]),
        c("logFC", "pval")
      ) == FALSE) {
        return(NULL)
      }
    }
  } else if (class(obj) == "data.frame") {
    if (length(obj) == 2) {
      if (checkmate::test_subset(
        names(obj),
        c("logFC", "pval")
      ) == FALSE) {
        return(NULL)
      }
      obj <- list(obj)
      names(obj) <- unlist(strsplit(nm,
                                         split = ".",
                                         fixed = TRUE
      ))[1]
    } else if (length(obj) == 6) {
      if (checkmate::test_subset(names(obj), c(
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B"
      )) == FALSE) {
        return(NULL)
      }
      obj <- deedee_prepare(obj, "limma")
      obj <- list(obj)
      names(obj) <- unlist(strsplit(nm,
                                         split = ".",
                                         fixed = TRUE
      ))[1]
    } else {
      return(NULL)
    }
  }

  return(obj)
}
