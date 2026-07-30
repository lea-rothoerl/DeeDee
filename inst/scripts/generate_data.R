library("macrophage")
library("DESeq2")
library("DeeDeeExperiment")

data("gse", package = "macrophage")
dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)

keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
dds_macrophage <- DESeq(dds_macrophage)

IFNg_naive <- results(dds_macrophage,
  contrast = c("condition", "IFNg", "naive"),
  lfcThreshold = 1, alpha = 0.05
)

IFNg_both <- results(dds_macrophage,
  contrast = c("condition", "IFNg_SL1344", "IFNg"),
  lfcThreshold = 1, alpha = 0.05
)

Salm_naive <- results(dds_macrophage,
  contrast = c("condition", "SL1344", "naive"),
  lfcThreshold = 1, alpha = 0.05
)

Salm_both <- results(dds_macrophage,
  contrast = c("condition", "IFNg_SL1344", "SL1344"),
  lfcThreshold = 1, alpha = 0.05
)

dde_macrophage <- DeeDeeExperiment(
  sce = dds_macrophage,
  de_results = list(
    IFNg_naive = IFNg_naive,
    IFNg_both  = IFNg_both,
    Salm_naive = Salm_naive,
    Salm_both  = Salm_both
  )
)

save(dde_macrophage, file = "data/dde_macrophage.RData", compress = "xz")
