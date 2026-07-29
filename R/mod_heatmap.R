#' @keywords internal
mod_heatmap_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(4, 8),
    bslib::card(
      shiny::numericInput(ns("show_first"), "Show first", value = 25, min = 1),
      shiny::checkboxInput(ns("show_gene_names"), "Show gene names", FALSE),
      shiny::checkboxInput(ns("show_na"), "Show NA", FALSE),
      shiny::selectInput(ns("dist"), "Distance measure",
                         choices = c("Euclidean" = "euclidean", "Manhattan" = "manhattan",
                                     "Pearson" = "pearson", "Spearman" = "spearman"),
                         selected = "euclidean"
      ),
      shiny::selectInput(ns("clust"), "Clustering method",
                         choices = c("Single" = "single", "Complete" = "complete",
                                     "Average" = "average", "Centroid" = "centroid"),
                         selected = "average"
      ),
      shiny::numericInput(ns("pthresh"), "P-value threshold",
                          value = 0.05, min = 0.01, max = 1, step = 0.01
      ),
      shiny::actionButton(ns("create"), "Create heatmap")
    ),
    bslib::card(
      shinycssloaders::withSpinner(
        InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(ns("ht"))
      )
    )
  )
}

#' @keywords internal
#' @param dde reactive DeeDeeExperiment from mod_input_server()
mod_heatmap_server <- function(id, dde) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    heatmap_output <- shiny::eventReactive(input$create, {
      shiny::req(dde(), input$show_first)
      shiny::validate(shiny::need(
        length(DeeDeeExperiment::getDEANames(dde())) >= 2,
        "Please select at least two contrasts."
      ))

      res <- deedee_heatmap(dde(),
                            show_first       = input$show_first,
                            show_gene_names  = input$show_gene_names,
                            dist             = input$dist,
                            clust            = input$clust,
                            pthresh          = input$pthresh,
                            show_na          = input$show_na
      )

      shiny::validate(shiny::need(!is.null(res), "No common genes in input datasets."))
      ComplexHeatmap::draw(res)
    })

    shiny::observeEvent(input$create, {
      shiny::req(length(DeeDeeExperiment::getDEANames(dde())) >= 2)
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input, output, session,
        heatmap_output(),
        heatmap_id = ns("ht")
      )
    })
  })
}
