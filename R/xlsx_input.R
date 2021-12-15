
xlsx_input <- function(obj,
                       path,
                       nm) {

  if (length(obj) > 1) {
    out <- lapply(obj,
                       readxl::read_excel,
                       path = path
    )
    names(out) <- obj
    for (j in 1:length(obj)) {
      out[[obj[j]]] <- as.data.frame(out[[obj[j]]])
      out[[obj[j]]] <- tibble::column_to_rownames(
        out[[obj[j]]], "rowname"
      )
      if (checkmate::test_subset(
        names(out[[j]]) == FALSE,
        c("logFC", "pval")
      )) {
        return(NULL)
      }
    }
  } else {
    out <- readxl::read_excel(path = path)
    out <- tibble::column_to_rownames(out, "rowname")
    out <- list(out)
    names(out) <- unlist(strsplit(nm,
                                       split = ".",
                                       fixed = TRUE
    ))[1]
  }
  return(out)
}

