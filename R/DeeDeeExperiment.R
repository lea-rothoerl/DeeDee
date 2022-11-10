
#' se: one SE object
#' de_results: a list of DE results, whatever format these are
#' will be handled by other function to split/convert/make uniform
#' has to have some names, otherwise will override with some defaults - ideally preferring names
#'
#' dde <- DeeDeeExperiment(se_macrophage, de_results = del)
#'
#' @param se SumExp
#' @param de_results ssss
#'
#' @return tODO
#' @export
#'
#' @examples
#' # todo
DeeDeeExperiment <- function(se,
                             de_results = NULL) {

  # se <- SummarizedExperiment(...)

  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }

  # TODO: if no SE is really provided, instantiate some rownames, at least directly
  # from the rownames of the result objects
  # TODO: the row names are taken from the FIRST object in the de results then - or
  # from the union of all of them?

  if (is.null(de_results)) {
    object <- new("DeeDeeExperiment",
                  se,
                  dea = list()
    )

    # stash the package version
    metadata(object)[["version"]] <- packageVersion("DeeDee")

    return(object)
  }

  # TODO: does not have to relate to an SE which has all the slots and all
  # ...


  # TODO: additional checks
  dde <- se

  # here is where I add the names in the rowData to make all info matched
  # checks on the names
  names(de_results)
  # if not there, "force add"
  # TODO

  dde_ids <- rownames(dde)

  dea_contrasts <- list()

  for (i in names(de_results)){
    this_de <- de_results[[i]]

    # do different things according to what these objects are
    if(is(this_de, "DESeqResults")) {
      matched_ids <- match(rownames(dde), dde_ids)

      # if not tested, add NA - everywhere? -> pre-fill?
      rowData(dde)[[paste0(i,"_log2FoldChange")]] <- NA
      rowData(dde)[[paste0(i,"_pvalue")]] <- NA
      rowData(dde)[[paste0(i,"_padj")]] <- NA

      rowData(dde)[[paste0(i,"_log2FoldChange")]][matched_ids] <- this_de$log2FoldChange
      rowData(dde)[[paste0(i,"_pvalue")]][matched_ids] <- this_de$pvalue
      rowData(dde)[[paste0(i,"_padj")]][matched_ids] <- this_de$padj

      dea_contrasts[[i]] <- list(
        alpha = metadata(this_de)$alpha,
        lfcThreshold = metadata(this_de)$lfcThreshold,
        metainfo_logFC = mcols(this_de)$description[colnames(this_de) == "log2FoldChange"],
        metainfo_pvalue = mcols(this_de)$description[colnames(this_de) == "pvalue"],
        original_object = this_de,
        package = "DESeq2"
      )
    }
  }

  # rowData(dde)[["new_rd"]] <- de_name

  object <- new("DeeDeeExperiment",
                dde,
                dea = dea_contrasts
  )

  # stash the package version
  metadata(object)[["version"]] <- packageVersion("DeeDee")

  return(object)
}






#' Title
#'
#' @param del todo
#'
#' @return todo
#' @export
#'
#' @examples
#' # todo
de_list_checker <- function(del) {
  # check that all things are as expected/expectable
}


.importDE_DESeq2 <- function(x) {
  # checks TODO:
  # correct object format
  # contain the right columns
  # contain the feature ids
  # p value respect the 0-1 interval

}


.importDE_edgeR <- function(x) {

}

.importDE_limma <- function(x) {

}

.importDE_custom <- function(x) {

}




deedee_import <- function(x) {

  # legacy code:

  # # ----------------------------- argument check ------------------------------
  # choices <- c("DESeq2", "edgeR", "limma")
  # checkmate::assertChoice(input_type, choices)
  #
  # if (input_type == "DESeq2") {
  #   checkmate::assertClass(data, "DESeqResults")
  #   logFC <- data$log2FoldChange
  #   pval <- data$padj
  #   input <- data.frame(logFC, pval)
  #   rownames(input) <- data@rownames
  # } else if (input_type == "edgeR") {
  #   checkmate::assertClass(data, "DGEExact")
  #   logFC <- data[["table"]][["logFC"]]
  #   pval <- data[["table"]][["PValue"]]
  #   input <- data.frame(logFC, pval)
  #   rownames(input) <- data[["genes"]][["genes"]]
  # } else if (input_type == "limma") {
  #   checkmate::assertDataFrame(data, types = "numeric")
  #   logFC <- data$logFC
  #   pval <- data$adj.P.Val
  #   input <- data.frame(logFC, pval)
  #   rownames(input) <- rownames(data)
  # }
  # return(input)
}


