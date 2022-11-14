#' Title
#'
#' @param dde todo
#' @param template todo
#' @param output_path todo
#' @param overwrite todo
#' @param pthresh todo
#' @param scatter_select1 todo
#' @param scatter_select2 todo
#' @param scatter_color_by todo
#' @param heatmap_show_first todo
#' @param heatmap_show_gene_names todo
#' @param heatmap_dist todo
#' @param heatmap_clust todo
#' @param heatmap_show_na todo
#' @param venn_mode todo
#' @param upset_mode todo
#' @param upset_min_setsize todo
#' @param qqmult_ref todo
#' @param cat_ref todo
#' @param cat_maxrank todo
#' @param cat_mode todo
#' @param render_quiet todo
#' @param silent todo
#' @param open_file todo
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
#'   de_results = del
#' )
#' dde
#'
#' # deedee_summary(dde)
deedee_summary <- function(dde,
                           template = system.file("extdata",
                                                  "summary_template.Rmd",
                                                  package = "DeeDee"),
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
                           cat_mode = "up",
                           render_quiet = TRUE,
                           silent = FALSE,
                           open_file = TRUE) {

  # checkmate::assert_list(deedee_list, type = "data.frame", min.len = 2)
  # for (i in 1:length(deedee_list)) {
  #   checkmate::assert_data_frame(deedee_list[[i]], type = "numeric")
  # }
  # checkmate::assert_number(pthresh, lower = 0, upper = 1)
  #
  # checkmate::assert_logical(overwrite)
  #
  # checkmate::assert_number(scatter_select1, lower = 1, upper = length(deedee_list))
  # checkmate::assert_number(scatter_select2, lower = 1, upper = length(deedee_list))
  # choices <- c("pval1", "pval2")
  # checkmate::assert_choice(scatter_color_by, choices)
  #
  # checkmate::assert_number(heatmap_show_first, lower = 1)
  # checkmate::assert_logical(heatmap_show_gene_names)
  # choices1 <- c("euclidean", "manhattan", "pearson", "spearman")
  # checkmate::assert_choice(heatmap_dist, choices1)
  # choices2 <- c("single", "complete", "average", "centroid")
  # checkmate::assert_choice(heatmap_clust, choices2)
  # checkmate::assert_logical(heatmap_show_na)
  #
  # choices <- c("up", "down", "both")
  # checkmate::assert_choice(venn_mode, choices)
  #
  # choices <- c("up", "down", "both", "both_colored")
  # checkmate::assert_choice(upset_mode, choices)
  # checkmate::assert_number(upset_min_setsize, lower = 0)
  #
  # checkmate::assert_number(qqmult_ref, lower = 1, upper = length(deedee_list))
  #
  # checkmate::assert_number(cat_ref, lower = 1, upper = length(deedee_list))
  # checkmate::assert_number(cat_maxrank, lower = 1)
  # choices <- c("up", "down", "both")
  # checkmate::assert_choice(cat_mode, choices)

  # template <- system.file("extdata",
  #                         "summary_template.Rmd",
  #                         package = "DeeDee"
  # )

  file.exists(template)

  if (file.exists(output_path)) {
    if (!overwrite) {
      stop("Your declared output file already exists. ",
           "Set overwrite = TRUE to proceed anyway.",
           call. = FALSE
      )
    }
  }

  output_rmd <- paste(unlist(strsplit(output_path,
                                      split = ".",
                                      fixed = TRUE
  ))[1], ".Rmd", sep = "")

  file.copy(from = template, to = output_rmd, overwrite = TRUE)

  args <- list()
  args$input <- output_rmd
  args$output_format <- "html_document"
  args$output_path <- output_path

  output_path <- rmarkdown::render(output_rmd,
                                   params = args,
                                   quiet = render_quiet
  )
  if (open_file == TRUE) {
    utils::browseURL(output_path)
  }

  file.remove(output_rmd)

  if (silent == FALSE) {
    message("Your summary has been generated!")
  }
}

# deedee_summary(dde,
#                template = "inst/extdata/summary_template.Rmd", output_path = "v2.html")


