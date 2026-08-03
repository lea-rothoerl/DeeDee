#' @keywords internal
mod_cat_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(4, 8),
    bslib::card(
      shiny::selectInput(ns("mode"), "Mode",
                         choices = c("Up" = "up", "Down" = "down", "Both" = "both"),
                         selected = "up"
      ),
      shiny::numericInput(ns("maxrank"), "Max rank", value = 1000, min = 1),
      shiny::uiOutput(ns("ref_choice")),
      shiny::numericInput(ns("pthresh"), "P-value threshold",
                          value = 0.05, min = 0.01, max = 1, step = 0.01
      )
    ),
    bslib::card(
      shinycssloaders::withSpinner(
        shiny::uiOutput(ns("cat_out"))
      )
    ),
    bslib::accordion(
      open = FALSE,
      bslib::accordion_panel(
        "INFO",
        shiny::includeMarkdown(system.file("extdata",
                                           "cat.md",
                                           package = "DeeDee"))
      )
    )
  )
}

#' @keywords internal
#' @param dde reactive DeeDeeExperiment from mod_input_server()
mod_cat_server <- function(id, dde) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$ref_choice <- shiny::renderUI({
      shiny::req(dde())
      contrasts <- DeeDeeExperiment::getDEANames(dde())
      shiny::selectInput(ns("ref"), "Reference contrast",
                         choices = contrasts, selected = contrasts[1]
      )
    })

    plot_obj <- shiny::reactive({
      shiny::req(dde(), input$ref, input$maxrank, input$mode)
      shiny::validate(shiny::need(
        length(DeeDeeExperiment::getDEANames(dde())) >= 2,
        "Please select at least two contrasts."
      ))

      ref_idx <- match(input$ref, DeeDeeExperiment::getDEANames(dde()))
      shiny::req(ref_idx)

      res <- deedee_cat(dde(),
                        ref = ref_idx,
                        maxrank = input$maxrank,
                        mode = input$mode,
                        pthresh = input$pthresh
      )

      shiny::validate(shiny::need(
        !is.null(res),
        "No genes in your datasets. Maybe your specified p-value threshold is too low?"
      ))
      res
    })

    output$cat_plotly <- plotly::renderPlotly({
      shiny::req(plot_obj())
      plotly::ggplotly(plot_obj())
    })

    output$cat_out <- shiny::renderUI({
      plotly::plotlyOutput(ns("cat_plotly"))
    })
  })
}
