suppressPackageStartupMessages(
  library("SummarizedExperiment")
)

data("de_named_list", package = "DeeDee")

rd_macrophage <- DataFrame(
  gene_id = rownames(del$ifng_vs_naive))
rownames(rd_macrophage) <- rownames(del$ifng_vs_naive)
se_macrophage_noassays <- SummarizedExperiment(
  assays = SimpleList(),
  rowData = rd_macrophage
)

names(del)

