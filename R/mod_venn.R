#' @keywords internal
mod_venn_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(4, 8),
    bslib::card(
      shiny::uiOutput(ns("select_contrasts")),
      shiny::selectInput(ns("mode"), "Mode",
                         choices = c("Up" = "up", "Down" = "down", "Both" = "both"),
                         selected = "both"
      ),
      shiny::numericInput(ns("pthresh"), "P-value threshold",
                          value = 0.05, min = 0.01, max = 1, step = 0.01
      )
    ),
    bslib::card(
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("venn"))
      )
    )
  )
}

#' @keywords internal
#' @param dde reactive DeeDeeExperiment from mod_input_server()
mod_venn_server <- function(id, dde) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$select_contrasts <- shiny::renderUI({
      shiny::req(dde())
      contrasts <- DeeDeeExperiment::getDEANames(dde())
      shiny::selectizeInput(ns("contrasts"), "Contrasts (2-4)",
                            choices = contrasts,
                            selected = utils::head(contrasts, 4),
                            multiple = TRUE,
                            options = list(maxItems = 4)
      )
    })

    output$venn <- shiny::renderPlot({
      shiny::req(dde(), input$contrasts, input$mode)
      shiny::validate(shiny::need(
        length(input$contrasts) >= 2,
        "Please select at least two contrasts."
      ))

      contrasts <- DeeDeeExperiment::getDEANames(dde())
      sel <- match(input$contrasts, contrasts)
      shiny::req(sel)

      res <- deedee_venn(dde(),
                         mode    = input$mode,
                         pthresh = input$pthresh,
                         select  = sel
      )

      shiny::validate(shiny::need(
        !is.null(res),
        "No genes in your datasets. Maybe your specified p-value threshold is too low?"
      ))
      res
    })
  })
}
