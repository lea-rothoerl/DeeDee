#' @keywords internal
mod_rrho_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(4, 8),
    bslib::card(
      shiny::uiOutput(ns("select1")),
      shiny::uiOutput(ns("select2")),
      shiny::selectInput(ns("alternative"), "Test type",
                         choices = c("Two-sided" = "two.sided", "Enrichment" = "enrichment"),
                         selected = "two.sided"
      ),
      shiny::checkboxInput(ns("by_correction"), "Benjamini-Yekutieli FDR correction", FALSE)
    ),
    bslib::card(
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("rrho"))
      )
    ),
    bslib::accordion(
      open = FALSE,
      bslib::accordion_panel(
        "INFO",
        shiny::includeMarkdown(system.file("extdata",
                                           "rrho.md",
                                           package = "DeeDee"))
      )
    )
  )
}

#' @keywords internal
#' @param dde reactive DeeDeeExperiment from mod_input_server()
mod_rrho_server <- function(id, dde) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$select1 <- shiny::renderUI({
      shiny::req(dde())
      contrasts <- DeeDeeExperiment::getDEANames(dde())
      shiny::selectInput(ns("s1"), "1st contrast", choices = contrasts)
    })
    output$select2 <- shiny::renderUI({
      shiny::req(dde())
      contrasts <- DeeDeeExperiment::getDEANames(dde())
      shiny::selectInput(ns("s2"), "2nd contrast",
                         choices = contrasts, selected = contrasts[2]
      )
    })

    rrho_obj <- shiny::reactive({
      shiny::req(dde(), input$s1, input$s2, input$alternative)
      shiny::validate(shiny::need(
        length(DeeDeeExperiment::getDEANames(dde())) >= 2,
        "Please select at least two contrasts."
      ))

      contrasts <- DeeDeeExperiment::getDEANames(dde())
      sel1 <- match(input$s1, contrasts)
      sel2 <- match(input$s2, contrasts)
      shiny::req(sel1, sel2)
      shiny::validate(shiny::need(sel1 != sel2, "Please select two different contrasts."))

      deedee_rrho(dde(),
                  select1 = sel1,
                  select2 = sel2,
                  alternative = input$alternative,
                  corr = input$by_correction
      )
    })

    output$rrho <- shiny::renderPlot({
      shiny::req(rrho_obj())
      mat <- if (isTRUE(input$by_correction)) {
        rrho_obj()$hypermat.by
      } else {
        rrho_obj()$hypermat
      }
      lattice::levelplot(mat,
                         col.regions = viridis::viridis(50, option = "magma"),
                         xlab = "", ylab = ""
      )
    })
  })
}
