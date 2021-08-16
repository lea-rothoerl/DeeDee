library("macrophage")
library("DESeq2")
source("~/Development/DeeDee_wip/deedee_prepare.R")

data(gse, "macrophage")

dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)

keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
dds_macrophage <- DESeq(dds_macrophage)

IFNg_naive <- results(dds_macrophage,
                      contrast = c("condition", "IFNg", "naive"),
                      lfcThreshold = 1, alpha = 0.05)

save(IFNg_naive, file = "data/DE_results_IFNg_naive.RData", compress = "xz")


IFNg_both <- results(dds_macrophage,
                     contrast = c("condition", "IFNg_SL1344", "IFNg"),
                     lfcThreshold = 1, alpha = 0.05)

save(IFNg_both, file = "data/DE_results_IFNg_both.RData", compress = "xz")


Salm_naive <- results(dds_macrophage,
                      contrast = c("condition", "SL1344", "naive"),
                      lfcThreshold = 1, alpha = 0.05)

save(Salm_naive, file = "data/DE_results_Salm_naive.RData", compress = "xz")


Salm_both <- results(dds_macrophage,
                     contrast = c("condition", "IFNg_SL1344", "SL1344"),
                     lfcThreshold = 1, alpha = 0.05)

save(Salm_both, file = "data/DE_results_Salm_both.RData", compress = "xz")
