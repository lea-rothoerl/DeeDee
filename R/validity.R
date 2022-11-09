validDeeDeeExperiment <- function(object) {
  msg <- NULL

  if (!is(object@dea, "list")) {
    msg <- c(msg, "'dea' must be a list")
  }

  if (length(dea(object)) > 0) {
    if (any(is.null(names(object@dea)))) {
      msg <- c(msg, "'dea' must be a named list")
    }

    dea_names <- names(object@dea)

    required_rowdata <- unlist(
      lapply(
        dea_names, function(arg)
          c(
            paste0(arg, "_log2FoldChange"),
            paste0(arg, "_pvalue"),
            paste0(arg, "_padj")
          )
      )
    )
    if (!all(required_rowdata %in% colnames(rowData(object)))) {
      msg <- c(msg, "some required columns were not found in the rowData")
    }
  }

  if (is.null(msg)) {
    TRUE
  } else msg


}

setValidity("DeeDeeExperiment", validDeeDeeExperiment)

