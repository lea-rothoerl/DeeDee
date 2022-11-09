
#' se: one SE object
#' de_results: a list of DE results, whatever format these are
#' will be handled by other function to split/convert/make uniform
#' has to have some names, otherwise will override with some defaults - ideally preferring names
#'
#' dde <- DeeDeeExperiment(se_macrophage, de_results = del)
#'
#' @param se ssss
#' @param de_results ssss
#'
#' @return tODO
#' @export
#'
#' @examples
#' # todo
DeeDeeExperiment <- function(se, de_results) {
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
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


.importDE_edgeR <- function(x) {

}

.importDE_limma <- function(x) {

}

.importDE_custom <- function(x) {

}



#' dde : the deedee exp object
#' dea: the named list, as for the main constructor
#' Title
#'
#' @param dde todo
#' @param dea todo
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
#'
#' new_del <- list(
#'   ifng2 = del$ifng_vs_naive,
#'   ifngsalmo2 = del$ifngsalmo_vs_naive
#' )
#' # add a new (set of) DE result(s)
#' dde_new <- add_dea(dde, new_del)
add_dea <- function(dde, dea) {
  # dde must be a DeeDeeExp
  # dea must be named list

  # check that names are all unique, and do not overlap with the existing ones
  names(dea)
  names(dea(dde))

  dea_contrasts <- dea(dde)
  dde_ids <- rownames(dde)

  # update rowData, naming them correctly
  for (i in names(dea)){
    this_de <- dea[[i]]

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
  # update the deslot
  ## TODO should be a replacement operator?
  dea(dde) <- dea_contrasts

  # here check some validity?
  # TODO

  # return the object
  return(dde)
}

# TODO: might need one where I also simply add ONE single DE object, and that gets autoconverted to a named list (of length 1)


#' Title
#'
#' @param dde  todo
#' @param dea_name todo
#'
#' @return todo
#' @export
#'
#' @examples
#' # todo
remove_dea <- function(dde, dea_name) {

}





#' Title
#'
#' @param dde TODO
#' @param dea_name TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
get_dea <- function(dde,
                    dea_name) {
  deas <- dea(dde)
  dea_names <- names(deas)

  if (!(dea_name %in% dea_names)) {
    stop("dea not found")
  }

  rd_info <- paste0(dea_name,
                    c("_log2FoldChange", "_pvalue", "_padj"))

  if (! all(rd_info %in% colnames(rowData(dde)))) {
    stop("Columns not found")
  }

  out <- rowData(dde)[, rd_info]

  return(out)
}


#' this reminds the format of the DeeDee list object
#'
#' @param dde todo
#'
#' @return todo
#' @export
#'
#' @examples
#' # todo
get_deas_list<- function(dde) {
  deas <- dea(dde)
  dea_names <- names(deas)

  dea_list <- list()

  for (i in dea_names) {
    dea_list[[i]] <- as.data.frame(get_dea(dde, i))
    colnames(dea_list[[i]]) <- c("log2FoldChange", "pvalue", "padj")
  }

  return(dea_list)
}

