#' DeeDee Prepare
#'
#' @description create DeeDee object from DEA result
#'
#' @param data result from DESeq2, limma or edgeR
#' @param input_type the program the data comes from ("DESeq2", "limma" or
#'                  "edgeR")
#'
#' @return table to be used as input for other DeeDee functions
#'
#' @examples
#'   library("macrophage")
#'   library("DESeq2")
#'   source("~/Development/DeeDee_wip/deedee_prepare.R")
#'
#'   data(gse, "macrophage")
#'
#'   dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
#'   rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#'
#'   keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
#'   dds_macrophage <- dds_macrophage[keep, ]
#'   dds_macrophage <- DESeq(dds_macrophage)
#'
#'   IFNg_naive <- results(dds_macrophage,
#'   contrast = c("condition", "IFNg", "naive"),
#'   lfcThreshold = 1, alpha = 0.05)
#'
#'   dd_table <- deedee_prepare(IFNg_naive, "DESeq2")
#'
#' @export
#'

deedee_prepare <- function(data, input_type) {

  # ----------------------------- argument check ------------------------------
  choices <- c("DESeq2", "edgeR", "limma")
  checkmate::assertChoice(input_type, choices)

  if (input_type == "DESeq2") {
    checkmate::assertClass(data, "DESeqResults")
    logFC <- data$log2FoldChange
    pval <- data$padj
    input <- data.frame(logFC, pval)
    rownames(input) <- data@rownames
  }

    else if (input_type == "edgeR") {
      checkmate::assertClass(data, "DGEExact")
      logFC <- data[["table"]][["logFC"]]
      pval <- data[["table"]][["PValue"]]
      input <- data.frame(logFC, pval)
      rownames(input) <- data[["genes"]][["genes"]]
  }

    else if (input_type == "limma") {
      checkmate::assertDataFrame(data, types = "numeric")
      logFC <- data$logFC
      pval <- data$adj.P.Val
      input <- data.frame(logFC, pval)
      rownames(input) <- rownames(data)
  }
  return(input)
}
