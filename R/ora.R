
ora <- function(geneList,
                universe,
                orgDB,
                key_type) {

  # ---------------------------- data preparation -----------------------------
  genes <- geneList[1]
  genes <- row.names(genes)
  genes <- genes[!is.na(genes)]

  universe[[1]] <- tibble::rownames_to_column(universe[[1]])
  universe[[2]] <- tibble::rownames_to_column(universe[[2]])
  univ <- dplyr::inner_join(universe[[1]], universe[[2]], by = "rowname")
  univ <- univ[["rowname"]]

  # ---------------------------------- GSEA -----------------------------------
  res <- clusterProfiler::enrichGO(
    gene = genes,
    universe = univ,
    OrgDb = get(orgDB),
    ont = "BP",
    keyType = key_type,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  # --------------------------------- return ----------------------------------
  return(res)
}
