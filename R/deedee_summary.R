#' DeeDee Summary
#'
#' @param deedee_list named list of results from deedee_prepare()
#' @param output_path the path to save the resulting report in. must end with
#'                    a filename.html (default = "DeeDee_Summary.html" in the
#'                    working directory)
#' @param overwrite logical value specifying if the output is supposed to
#'                  overwrite a potential existing file at the location of
#'                  output_path (default = FALSE)
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#' @param scatter_select1 index of first data-list element to be used
#'                        (default = 1)
#' @param scatter_select2 index of second data-list element to be used
#'                        (default = 2)
#' @param scatter_color_by indicates which set of values the output should be
#'                         colored by (possible values = `pval1` (default),
#'                         `pval2`)
#' @param heatmap_show_first indicating the number of genes depicted
#'                           (default = 25)
#' @param heatmap_show_gene_names boolean, show row names next to heatmap
#'                                (default = FALSE)
#' @param heatmap_dist select the distance measure (`euclidean`, `manhattan`,
#'                     `pearson`, `spearman`)
#' @param heatmap_clust select the clustering method (`single`, `complete`,
#'                      `average`, `centroid`)
#' @param heatmap_show_na boolean, include genes with NA values in heatmap
#'                        (default = FALSE)
#' @param venn_mode show all overlapping DE genes (`both`, default),
#'                  only conjointly up-regulated (`up`)
#'                  or only conjointly down-regulated (`down`) genes
#' @param upset_mode show all overlapping DE genes (`both`),
#'                   all overlapping genes colored by DE direction
#'                   (`both_colored`, default),
#'                   only conjointly up-regulated (`up`)
#'                   or only conjointly down-regulated (`down`) genes
#' @param upset_min_setsize the minimum size of intersections to be displayed in
#'                          the UpSet plot (default = 10)
#' @param qqmult_ref index of the contrast in data to be used as reference
#'                   contrast (default = 1)
#' @param cat_ref index of the contrast in data to be used as reference contrast
#'                (default = 1)
#' @param cat_maxrank highest rank that should be displayed (default = 1000)
#' @param cat_mode sort by highest logFC (`up`, default), lowest logFC (`down`)
#'                 or greatest deviation from zero (`both`)
#'
#' @return Creates a html document containing the results of running the DeeDee
#'         functions on your input data and params. TODO where
#'
#' @export
#'
#' @examples
#'
#' data(DE_results_IFNg_naive, package = "DeeDee")
#' IFNg_naive <- deedee_prepare(IFNg_naive, "DESeq2")
#'
#' data(DE_results_IFNg_both, package = "DeeDee")
#' IFNg_both <- deedee_prepare(IFNg_both, "DESeq2")
#'
#' data(DE_results_Salm_naive, package = "DeeDee")
#' Salm_naive <- deedee_prepare(Salm_naive, "DESeq2")
#'
#' data(DE_results_Salm_both, package = "DeeDee")
#' Salm_both <- deedee_prepare(Salm_both, "DESeq2")
#'
#' dd_list <- list(
#'   IFNg_naive = IFNg_naive, IFNg_both = IFNg_both,
#'   Salm_naive = Salm_naive, Salm_both = Salm_both
#' )
#' \dontrun{
#' deedee_summary(dd_list)
#' }
#'
deedee_summary <- function(deedee_list,
                           output_path = "DeeDee_Summary.html",
                           overwrite = FALSE,
                           pthresh = 0.05,
                           scatter_select1 = 1,
                           scatter_select2 = 2,
                           scatter_color_by = "pval1",
                           heatmap_show_first = 25,
                           heatmap_show_gene_names = FALSE,
                           heatmap_dist = "euclidean",
                           heatmap_clust = "average",
                           heatmap_show_na = FALSE,
                           venn_mode = "both",
                           upset_mode = "both_colored",
                           upset_min_setsize = 10,
                           qqmult_ref = 1,
                           cat_ref = 1,
                           cat_maxrank = 1000,
                           cat_mode = "up") {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(deedee_list, type = "data.frame", min.len = 2)
  for (i in 1:length(deedee_list)) {
    checkmate::assert_data_frame(deedee_list[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0, upper = 1)

  checkmate::assert_logical(overwrite)

  checkmate::assert_number(scatter_select1, lower = 1, upper = length(deedee_list))
  checkmate::assert_number(scatter_select2, lower = 1, upper = length(deedee_list))
  choices <- c("pval1", "pval2")
  checkmate::assert_choice(scatter_color_by, choices)

  checkmate::assert_number(heatmap_show_first, lower = 1)
  checkmate::assert_logical(heatmap_show_gene_names)
  choices1 <- c("euclidean", "manhattan", "pearson", "spearman")
  checkmate::assert_choice(heatmap_dist, choices1)
  choices2 <- c("single", "complete", "average", "centroid")
  checkmate::assert_choice(heatmap_clust, choices2)
  checkmate::assert_logical(heatmap_show_na)

  choices <- c("up", "down", "both")
  checkmate::assert_choice(venn_mode, choices)

  choices <- c("up", "down", "both", "both_colored")
  checkmate::assert_choice(upset_mode, choices)
  checkmate::assert_number(upset_min_setsize, lower = 0)

  checkmate::assert_number(qqmult_ref, lower = 1, upper = length(deedee_list))

  checkmate::assert_number(cat_ref, lower = 1, upper = length(deedee_list))
  checkmate::assert_number(cat_maxrank, lower = 1)
  choices <- c("up", "down", "both")
  checkmate::assert_choice(cat_mode, choices)

  # ---------------------------- R Markdown output -----------------------------
  template <- system.file("extdata",
    "summary_template.Rmd",
    package = "DeeDee"
  )

  if (file.exists(output_path)) {
    if (!overwrite) {
      stop(" Your declared output file already exists. ",
        "Set overwrite = TRUE to proceed anyway.",
        call. = FALSE
      )
    }
  }

  output_rmd <- unlist(strsplit(output_path,
    split = ".",
    fixed = TRUE
  ))[1]

  file.copy(from = template, to = output_rmd, overwrite = TRUE)

  args <- list()
  args$input <- output_rmd
  args$output_format <- "html_document"
  args$output_path <- output_path

  output_path <- rmarkdown::render(output_rmd,
    params = args
  )
  utils::browseURL(output_path)

  file.remove(output_rmd)

  print("Your summary has been generated!")
}
