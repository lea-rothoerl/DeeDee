
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
DeeDeeExperiment <- function(se = NULL,
                             de_results = NULL) {

  # old <- S4Vectors:::disableValidity()
  # if (!isTRUE(old)) {
  #   S4Vectors:::disableValidity(TRUE)
  #   on.exit(S4Vectors:::disableValidity(old))
  # }

  if (!is.null(de_results))
    .check_de_results(de_results)

  # se <- SummarizedExperiment(...)
  if(!is.null(se)) {
    if (!is(se, "RangedSummarizedExperiment")) {
      if (is(se, "SummarizedExperiment")) {
        se <- as(se, "RangedSummarizedExperiment")
      } else {
        stop("'se' must be a RangedSummarizedExperiment object")
      }
    }
  } else {
    if(is.null(de_results))
      stop("You have to provide at least an se object or a de_results object")

    message("TODO: needs to create a mock SE that does not break the creation later?!")
    ## mock - todo ## # mock up the se from the de_results
    ## mock - todo ## first_de <- de_results[[1]]
    ## mock - todo ##
    ## mock - todo ## # independently of the class, the feature names are in the rownames slot, TODO: check
    ## mock - todo ## ids <- rownames(first_de)
    ## mock - todo ##
    ## mock - todo ## rd_mock <- DataFrame(
    ## mock - todo ##   gene_id = ids,
    ## mock - todo ##   row.names = ids)
    ## mock - todo ##
    ## mock - todo ## # way1
    ## mock - todo ## se_mock <- SummarizedExperiment(
    ## mock - todo ##   assays = SimpleList(),
    ## mock - todo ##   rowData = rd_mock
    ## mock - todo ## )
    ## mock - todo ## # se_mock@NAMES <- NULL
    ## mock - todo ## # rownames(se_mock) <- ids
    ## mock - todo ##
    ## mock - todo ## se <- se_mock
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
  se_out <- se

  # here is where I add the names in the rowData to make all info matched
  # checks on the names
  names(de_results)
  # if not there, "force add"
  # TODO

  dde_ids <- rownames(se_out)

  dea_contrasts <- list()

  for (i in names(de_results)){
    this_de <- de_results[[i]]

    # do different things according to what these objects are
    if(is(this_de, "DESeqResults")) {


      input_deseq2 <- .importDE_DESeq2(se_out, this_de, i)
      se_out <- input_deseq2$se
      dea_contrasts[[i]] <- input_deseq2$dea_contrast

      # matched_ids <- match(rownames(se_out), rownames(this_de))
      #
      # # if not tested, add NA - everywhere? -> pre-fill?
      # rowData(se_out)[[paste0(i,"_log2FoldChange")]] <- NA
      # rowData(se_out)[[paste0(i,"_pvalue")]] <- NA
      # rowData(se_out)[[paste0(i,"_padj")]] <- NA
      #
      # rowData(se_out)[[paste0(i,"_log2FoldChange")]][matched_ids] <- this_de$log2FoldChange
      # rowData(se_out)[[paste0(i,"_pvalue")]][matched_ids] <- this_de$pvalue
      # rowData(se_out)[[paste0(i,"_padj")]][matched_ids] <- this_de$padj
      #
      # dea_contrasts[[i]] <- list(
      #   alpha = metadata(this_de)$alpha,
      #   lfcThreshold = metadata(this_de)$lfcThreshold,
      #   metainfo_logFC = mcols(this_de)$description[colnames(this_de) == "log2FoldChange"],
      #   metainfo_pvalue = mcols(this_de)$description[colnames(this_de) == "pvalue"],
      #   original_object = this_de,
      #   package = "DESeq2"
      # )
    }
  }

  # rowData(dde)[["new_rd"]] <- de_name

  object <- new("DeeDeeExperiment",
                se_out,
                dea = dea_contrasts
  )

  # stash the package version
  metadata(object)[["version"]] <- packageVersion("DeeDee")

  return(object)
}





#' extends the rowData slot of the provided SE and returns also metadata
.importDE_DESeq2 <- function(se, res_de, de_name) {
  # checks TODO:
  # correct object format
  # contain the right columns
  # contain the feature ids
  # p value respect the 0-1 interval


  matched_ids <- match(rownames(se), rownames(res_de))

  # if not tested, add NA - everywhere? -> pre-fill?
  rowData(se)[[paste0(de_name,"_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name,"_pvalue")]] <- NA
  rowData(se)[[paste0(de_name,"_padj")]] <- NA

  rowData(se)[[paste0(de_name,"_log2FoldChange")]][matched_ids] <- res_de$log2FoldChange
  rowData(se)[[paste0(de_name,"_pvalue")]][matched_ids] <- res_de$pvalue
  rowData(se)[[paste0(de_name,"_padj")]][matched_ids] <- res_de$padj

  dea_contrast <- list(
    alpha = metadata(res_de)$alpha,
    lfcThreshold = metadata(res_de)$lfcThreshold,
    metainfo_logFC = mcols(res_de)$description[colnames(res_de) == "log2FoldChange"],
    metainfo_pvalue = mcols(res_de)$description[colnames(res_de) == "pvalue"],
    original_object = res_de,
    # object_name = deparse(substitute(res_de)),
    package = "DESeq2"
  )

  return(
    list(
      se = se,
      dea_contrast = dea_contrast
    )
  )
}


.importDE_edgeR <- function(x) {

}

.importDE_limma <- function(x) {

}

.importDE_custom <- function(x) {

}



.check_de_results <- function(x) {
  # checks that:
  # it is a list
  # its component are either of the expected/accepted elements
  # checks their columns?

  stopifnot(is.list(x))
  stopifnot(length(x) > 0)

  ok_types <- unlist(lapply(x, function(arg) {
    is(arg, "DESeqResults") | is(arg, "DGEExact")
  }))

  stopifnot(all(ok_types))

  return(TRUE)
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


