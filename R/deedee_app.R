#' DeeDee App
#'
#' @description `deedee_app` opens the DeeDee Shiny web application, combining
#' the functionalities of all other DeeDee functions with a user-friendly
#' graphical user interface.
#'
#' @param deedee_obj An object of the class DeeDeeObject to be analyzed.
#'
#' @return A shiny app
#' @export
#'
#' @examples
#'
#' data(dde_macrophage, package = "DeeDee")
#' deedee_app(dde_macrophage)
#'
deedee_app <- function(deedee_obj = NULL) {

  ui <- bslib::page_navbar(
    title = "DeeDee",
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    bslib::nav_panel("Input", mod_input_ui("input")),
    bslib::nav_menu("Overlap",
                    bslib::nav_panel("Venn", mod_venn_ui("venn")),
                    bslib::nav_panel("UpSet", mod_upset_ui("upset")),
                    bslib::nav_panel("Overlap List", mod_overlap_ui("overlap")),
    ),
    bslib::nav_menu("Concordance",
                    bslib::nav_panel("CAT", mod_cat_ui("cat")),
                    bslib::nav_panel("RRHO", mod_rrho_ui("rrho")),
                    bslib::nav_panel("QQ", mod_qq_ui("qq")),
    ),
    bslib::nav_panel("Scatter", mod_scatter_ui("scatter")),

    bslib::nav_panel("Heatmap",
                     bslib::layout_columns(
                       col_widths = c(4, 8),
                       bslib::card(
                         shiny::numericInput("heatmap_show_first", "Show first", value = 25, min = 1),
                         shiny::checkboxInput("heatmap_show_gene_names", "Show gene names", FALSE),
                         shiny::checkboxInput("heatmap_show_na", "Show NA", FALSE),
                         shiny::selectInput("heatmap_dist", "Distance measure",
                                            choices = c("Euclidean" = "euclidean",
                                                        "Manhattan" = "manhattan",
                                                        "Pearson" = "pearson",
                                                        "Spearman" = "spearman"),
                                            selected = "euclidean"
                         ),
                         shiny::selectInput("heatmap_clust", "Clustering method",
                                            choices = c("Single" = "single",
                                                        "Complete" = "complete",
                                                        "Average" = "average",
                                                        "Centroid" = "centroid"),
                                            selected = "average"
                         ),
                         shiny::numericInput("heatmap_pthresh", "P-value threshold",
                                             value = 0.05, min = 0.01, max = 1, step = 0.01
                         )
                       ),
                       bslib::card(
                         shinycssloaders::withSpinner(
                           InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput("heatmap_ht")
                         )
                       ),
                       bslib::accordion(
                         open = FALSE,
                         bslib::accordion_panel(
                           "INFO",
                           shiny::includeMarkdown(system.file("extdata",
                                                              "heatmap.md",
                                                              package = "DeeDee"))
                         )
                       )
                     )
    )
  )

  server <- function(input, output, session) {
    input_mod <- mod_input_server("input", dde_arg = deedee_obj)

    dde <- input_mod$dde
    show_symbols <- input_mod$show_symbols
    species <- input_mod$species

    mod_venn_server("venn", dde)
    mod_upset_server("upset", dde)
    mod_overlap_server("overlap", dde, show_symbols, species)
    mod_cat_server("cat", dde)
    mod_rrho_server("rrho", dde)
    mod_scatter_server("scatter", dde, show_symbols, species)
    mod_qq_server("qq", dde)

    heatmap_output <- shiny::reactive({
      shiny::req(dde(), input$heatmap_show_first)
      shiny::validate(shiny::need(
        length(DeeDeeExperiment::getDEANames(dde())) >= 2,
        "Please select at least two contrasts."
      ))
      res <- deedee_heatmap(dde(),
                            show_first = input$heatmap_show_first,
                            show_gene_names = input$heatmap_show_gene_names,
                            dist = input$heatmap_dist,
                            clust = input$heatmap_clust,
                            pthresh = input$heatmap_pthresh,
                            show_na = input$heatmap_show_na
      )
      shiny::validate(shiny::need(!is.null(res),
                                  "No common genes in input datasets."))
      ComplexHeatmap::draw(res)
    })

    shiny::observe({
      shiny::req(heatmap_output())
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input, output, session, heatmap_output()
      )
    })
  }

  shiny::shinyApp(ui, server)
}
