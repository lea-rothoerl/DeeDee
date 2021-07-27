
gsea <- function(geneList,
                 universe,
                 orgDB,
                 select = 1) {

  # ---------------------------- data preparation -----------------------------
  if (select == 1) {
    genes <- geneList[1]
  }
  else {
    genes <- geneList[3]
  }
  genes <- row.names(genes)

  universe[[1]] <- tibble::rownames_to_column(universe[[1]])
  universe[[2]] <- tibble::rownames_to_column(universe[[2]])
  univ <- dplyr::inner_join(universe[[1]], universe[[2]], by = "rowname")
  univ <- univ[["rowname"]]

  # ---------------------------------- GSEA -----------------------------------
  res <- clusterProfiler::enrichGO(gene = genes,
                                   universe = univ,
                                   OrgDb = get(orgDB),
                                   ont = "BP",
                                   keyType = "ENSEMBL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.01,
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)

  # --------------------------------- return ----------------------------------
  return(res)
}
