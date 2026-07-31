#' @keywords internal
mod_upset_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(4, 8),
    bslib::card(
      shiny::selectInput(ns("mode"), "Mode",
                         choices = c("Up" = "up", "Down" = "down", "Both" = "both"),
                         selected = "both"
      ),
      shiny::conditionalPanel(
        condition = "input.mode == 'both'", ns = ns,
        shiny::checkboxInput(ns("colored"), "Coloring", TRUE)
      ),
      shiny::numericInput(ns("min_setsize"), "Minimum set size",
                          value = 10, min = 0, step = 1
      ),
      shiny::numericInput(ns("pthresh"), "P-value threshold",
                          value = 0.05, min = 0.01, max = 1, step = 0.01
      )
    ),
    bslib::card(
      shinycssloaders::withSpinner(
        shiny::uiOutput(ns("upset_out"))
      )
    )
  )
}

#' @keywords internal
#' @param dde reactive DeeDeeExperiment from mod_input_server()
mod_upset_server <- function(id, dde) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_obj <- shiny::reactive({
      shiny::req(dde(), input$mode)
      shiny::validate(shiny::need(
        length(DeeDeeExperiment::getDEANames(dde())) >= 2,
        "Please select at least two contrasts."
      ))

      mode <- if (input$mode == "both" && isTRUE(input$colored)) {
        "both_colored"
      } else {
        input$mode
      }

      res <- deedee_upset(dde(),
                          mode = mode,
                          pthresh = input$pthresh,
                          min_setsize = input$min_setsize
      )

      shiny::validate(shiny::need(
        !is.null(res),
        "No genes in your datasets. Maybe your specified p-value threshold is too low?"
      ))
      res
    })

    output$upset_static <- shiny::renderPlot({
      shiny::req(plot_obj())
      plot_obj()
    })

    output$upset_out <- shiny::renderUI({
      shiny::plotOutput(ns("upset_static"))
    })
  })
}
