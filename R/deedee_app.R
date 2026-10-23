#' DeeDee App
#'
#' @description `deedee_app` opens the DeeDee Shiny web application, combining
#' the functionalities of all other DeeDee functions with a user-friendly
#' graphical user interface.
#'
#' @param deedee_obj An object of the class DeeDeeExperiment to be analyzed. If
#'                   None, uploading a DeeDeeExperiment object as RDS via the
#'                   GUI is still possible.
#'
#' @return A shiny app
#' @export
#'
#' @examples
#'
#' if (interactive()) {
#'   data(dde_macrophage, package = "DeeDee")
#'   deedee_app(dde_macrophage)
#' }
#'
deedee_app <- function(deedee_obj = NULL) {

  arg_label <- paste(deparse(substitute(deedee_obj)), collapse = " ")

  ui <- bslib::page_navbar(
    title = "DeeDee",
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    header = shiny::tags$head(shiny::tags$script(shiny::HTML(
      "function copyDeeDeeCode(elId, btn) {
         var text = document.getElementById(elId).innerText;
         var done = function() {
           var old = btn.innerText;
           btn.innerText = 'Copied!';
           setTimeout(function() { btn.innerText = old; }, 1500);
         };
         if (navigator.clipboard && navigator.clipboard.writeText) {
           navigator.clipboard.writeText(text).then(done);
         } else {
           var ta = document.createElement('textarea');
           ta.value = text;
           document.body.appendChild(ta);
           ta.select();
           document.execCommand('copy');
           document.body.removeChild(ta);
           done();
         }
       }"
    ))),
    bslib::nav_panel("Input",
                     mod_input_ui("input",
                                  collapse_upload = !is.null(deedee_obj))),
    bslib::nav_menu("Overlap",
                    bslib::nav_panel("Venn", mod_venn_ui("venn")),
                    bslib::nav_panel("UpSet", mod_upset_ui("upset")),
                    bslib::nav_panel("Overlap List", mod_overlap_ui("overlap")),
    ),
    bslib::nav_menu("Concordance",
                    bslib::nav_panel("CAT", mod_cat_ui("cat")),
   #                 bslib::nav_panel("RRHO", mod_rrho_ui("rrho")),
                    bslib::nav_panel("QQ", mod_qq_ui("qq")),
    ),
    bslib::nav_panel("Scatter", mod_scatter_ui("scatter")),

    bslib::nav_panel("Heatmap",
                     bslib::layout_columns(
                       col_widths = c(4, 8),
                       bslib::card(
                         shiny::numericInput("heatmap_show_first",
                                             "Show first",
                                             value = 25,
                                             min = 1),
                         shiny::checkboxInput("heatmap_show_gene_names",
                                              "Show gene identifiers",
                                              FALSE),
                         shiny::checkboxInput("heatmap_show_na",
                                              "Show NA",
                                              FALSE),
                         shiny::selectInput("heatmap_dist",
                                            "Distance measure",
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
                         shiny::numericInput("heatmap_pthresh",
                                             "P-value threshold",
                                             value = 0.05,
                                             min = 0.01,
                                             max = 1,
                                             step = 0.01
                         )
                       ),
                       bslib::card(
                         shiny::tags$style(
                           shiny::HTML("#heatmap_ht_show_code { display: none !important; }")
                         ),
                         shinycssloaders::withSpinner(
                           InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput("heatmap_ht")
                         ),
                         .deedee_download_ui(identity, id_prefix = "heatmap_")
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
    input_mod <- mod_input_server("input",
                                  dde_arg = deedee_obj,
                                  arg_label = arg_label)

    dde <- input_mod$dde
    show_symbols <- input_mod$show_symbols
    species <- input_mod$species
    source_label <- input_mod$source_label
    code_data <- input_mod$code_data

    mod_venn_server("venn", dde, source_label, code_data)
    mod_upset_server("upset", dde, source_label, code_data)
    mod_overlap_server("overlap", dde, show_symbols, species)
    mod_cat_server("cat", dde, source_label, code_data)
    mod_scatter_server("scatter",
                       dde,
                       show_symbols,
                       species,
                       source_label,
                       code_data)
    mod_qq_server("qq", dde, source_label, code_data)

    heatmap_obj <- shiny::reactive({
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
                            show_na = input$heatmap_show_na,
                            show_symbols = show_symbols()
      )
      shiny::validate(shiny::need(!is.null(res),
                                  "No common genes in input datasets."))
      res
    })

    heatmap_output <- shiny::reactive({
      shiny::req(heatmap_obj())
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
      ComplexHeatmap::draw(heatmap_obj())
    })

    shiny::observe({
      shiny::req(heatmap_output())
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input, output, session, heatmap_output()
      )
    })

    .deedee_download_server(input, output,
                            id_prefix = "heatmap_",
                            filename_prefix = "deedee_heatmap",
                            draw_fn = function() ComplexHeatmap::draw(heatmap_obj()),
                            source_label = source_label,
                            width = 8, height = 8
    )

    .deedee_code_server(input, session, id_prefix = "heatmap_",
                        code_fn = function() {
                          cd <- code_data()
                          .deedee_code_text("deedee_heatmap",
                                            cd$data,
                                            args = list(
                                              show_first = input$heatmap_show_first,
                                              show_gene_names = input$heatmap_show_gene_names,
                                              dist = input$heatmap_dist,
                                              clust = input$heatmap_clust,
                                              pthresh = input$heatmap_pthresh,
                                              show_na = input$heatmap_show_na,
                                              show_symbols = show_symbols()
                                            ),
                                            preamble = cd$preamble,
                                            show = "ComplexHeatmap::draw(p)"
                          )
                        }
    )
  }

  shiny::shinyApp(ui, server, options = list(launch.browser = TRUE))
}
