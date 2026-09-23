#' DeeDee Functional Enrichment (on Scatter)
#'
#' @description `deedee_fea` runs `topGO` analysis (via `mosdef::run_topGO()`)
#' on a set of genes, e.g. genes selected interactively from the
#' `deedee_scatter()` plot. The universe for the test is every gene tested in
#' *both* selected contrasts.
#'
#' @param data instance of the DeeDeeExperiment class
#' @param select1 index of first data-list element to be used (default = 1)
#' @param select2 index of second data-list element to be used (default = 2)
#' @param genes character vector of gene symbols to test for enrichment
#'              (e.g. the genes selected from a `deedee_scatter()` plot)
#' @param species species the data comes from (possible values = `Homo_sapiens`
#'              (default), `Mus_musculus`)
#' @param ontology which GO domain to test (possible values = `BP` (default),
#'                 `MF`, `CC`)
#'
#' @return a data.frame of GO enrichment results
#'

deedee_fea <- function(data,
                       select1 = 1,
                       select2 = 2,
                       genes,
                       species = "Homo_sapiens",
                       ontology = "BP") {

  # ----------------------------- argument check ------------------------------
  checkmate::assert(length(names(data@dea)) >= 2)
  data <- deedee_from_dde(data)
  checkmate::assert_number(select1, lower = 1, upper = length(data))
  checkmate::assert_number(select2, lower = 1, upper = length(data))
  checkmate::assert_character(genes, min.len = 1)
  checkmate::assert_choice(species, c("Homo_sapiens", "Mus_musculus"))
  checkmate::assert_choice(ontology, c("BP", "MF", "CC"))

  mapping <- switch(species,
                    "Homo_sapiens" = "org.Hs.eg.db",
                    "Mus_musculus" = "org.Mm.eg.db"
  )

  for (pkg in c("AnnotationDbi", "topGO")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", mapping, "' is required for functional enrichment analysis ",
        "but not installed."
      )
    }
    if (!paste0("package:", pkg) %in% search()) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }

  # --------------------------------- universe --------------------------------
  bg <- dplyr::inner_join(
    tibble::rownames_to_column(data[[select1]]),
    tibble::rownames_to_column(data[[select2]]),
    by = "rowname"
  )
  bg_genes <- unique(stats::na.omit(bg$symbol.x))

  de_genes <- unique(stats::na.omit(genes[genes %in% bg_genes]))

  if (length(de_genes) < 2) {
    return(NULL)
  }

  # --------------------------------- run FEA ----------------------------------
  res <- mosdef::run_topGO(
    de_genes = de_genes,
    bg_genes = bg_genes,
    ontology = ontology,
    mapping = mapping,
    gene_id = "symbol",
    verbose = FALSE
  )

  # --------------------------------- return ----------------------------------
  return(res)
}

# --- SMALL HELPER ---
#' @keywords internal
beautify_GO_table <- function(enrich_tbl, GO_id_column = "GO.ID") {
  rownames(enrich_tbl) <- enrich_tbl[[GO_id_column]]
  enrich_tbl[[GO_id_column]] <- mosdef::create_link_GO(enrich_tbl[[GO_id_column]])
  return(enrich_tbl)
}
