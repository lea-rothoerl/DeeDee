#' @keywords internal
mod_scatter_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(4, 8),
    bslib::card(
      shiny::uiOutput(ns("select1")),
      shiny::uiOutput(ns("select2")),
      shiny::selectInput(ns("color_by"), "Color by",
                         choices = c("1st p-value" = "pval1",
                                     "2nd p-value" = "pval2")
      ),
      shiny::numericInput(ns("pthresh"), "P-value threshold",
                          value = 0.05, min = 0.01, max = 1, step = 0.01
      ),
    ),
    bslib::card(
      shinycssloaders::withSpinner(
        shiny::uiOutput(ns("scatter_out"))
      )
    )
  )
}

#' @keywords internal
#' @param dde reactive DeeDeeExperiment from mod_input_server()
mod_scatter_server <- function(id, dde) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$select1 <- shiny::renderUI({
      shiny::selectInput(ns("s1"), "1st contrast", choices = names(dde()@dea))
    })
    output$select2 <- shiny::renderUI({
      shiny::selectInput(ns("s2"), "2nd contrast",
                         choices = names(dde()@dea), selected = names(dde()@dea)[2]
      )
    })

    plot_obj <- shiny::reactive({
      shiny::req(input$s1, input$s2)
      shiny::validate(shiny::need(
        length(names(dde()@dea)) >= 2, "Please provide at least two contrasts."
      ))

      sel1 <- match(input$s1, names(dde()@dea))
      sel2 <- match(input$s2, names(dde()@dea))
      shiny::req(sel1, sel2)

      deedee_scatter(dde(),
                     select1 = sel1,
                     select2 = sel2,
                     color_by = input$color_by,
                     pthresh = input$pthresh
      )
    })

    output$scatter_plotly <- plotly::renderPlotly({
      shiny::req(plot_obj())
      plotly::ggplotly(plot_obj(), tooltip = c("text", "x", "y"))
    })

    output$scatter_out <- shiny::renderUI({
      plotly::plotlyOutput(ns("scatter_plotly"))
    })
  })
}
