#' @keywords internal
mod_overlap_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(4, 8),
    bslib::card(
      shiny::uiOutput(ns("select_contrasts")),
      shiny::selectInput(ns("mode"), "Mode",
                         choices = c("Up" = "up",
                                     "Down" = "down",
                                     "Both" = "both"),
                         selected = "both"
      ),
      shiny::numericInput(ns("pthresh"), "P-value threshold",
                          value = 0.05, min = 0.01, max = 1, step = 0.01
      ),
      shiny::selectInput(ns("species"), "Species (for gene links)",
                         choices = c("Human" = "Homo_sapiens",
                                     "Mouse" = "Mus_musculus",
                                     "Rat" = "Rattus_norvegicus"),
                         selected = "Homo_sapiens"
      )
    ),
    bslib::card(
      bslib::card_header("Overlapping genes"),
      shinycssloaders::withSpinner(
        DT::DTOutput(ns("table"))
      )
    ),
    bslib::accordion(
      open = FALSE,
      bslib::accordion_panel(
        "INFO",
        shiny::includeMarkdown(system.file("extdata",
                                           "overlap.md",
                                           package = "DeeDee"))
      )
    )
  )
}

#' @keywords internal
#' @param dde reactive DeeDeeExperiment from mod_input_server()
mod_overlap_server <- function(id, dde) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$select_contrasts <- shiny::renderUI({
      shiny::req(dde())
      contrasts <- DeeDeeExperiment::getDEANames(dde())
      shiny::selectizeInput(ns("contrasts"), "Contrasts",
                            choices = contrasts,
                            selected = utils::head(contrasts),
                            multiple = TRUE
      )
    })

    sel <- shiny::reactive({
      shiny::req(dde(), input$contrasts)
      contrasts <- DeeDeeExperiment::getDEANames(dde())
      match(input$contrasts, contrasts)
    })

    overlap_tbl <- shiny::reactive({
      shiny::req(sel(), input$mode, input$pthresh)
      shiny::validate(shiny::need(length(sel()) >= 2, "Please select at least two contrasts."))

      res <- deedee_overlap(dde(), mode = input$mode, pthresh = input$pthresh, select = sel())
      shiny::validate(shiny::need(
        !is.null(res),
        "No overlapping genes found. Maybe your specified p-value threshold is too low?"
      ))
      res
    })

    output$table <- DT::renderDT({
      tbl <- overlap_tbl()
      shiny::req(input$species)

      tbl$Ensembl <- mosdef::create_link_ENSEMBL(tbl$gene, species = input$species)
      tbl <- tbl[, c("Ensembl", setdiff(names(tbl), c("gene", "Ensembl"))), drop = FALSE]

      DT::datatable(tbl,
                    escape = FALSE,
                    rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 25)
      ) |>
        DT::formatRound(columns = grep("_logFC$", names(tbl), value = TRUE), digits = 3) |>
        DT::formatSignif(columns = grep("_pval$", names(tbl), value = TRUE), digits = 3)
    })
  })
}
