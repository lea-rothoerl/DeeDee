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
    ),
    bslib::accordion(
      open = FALSE,
      bslib::accordion_panel(
        "INFO",
        shiny::includeMarkdown(system.file("extdata",
                                           "scatter.md",
                                           package = "DeeDee"))
      )
    ),
    bslib::accordion(
      open = FALSE,
      bslib::accordion_panel(
        "Functional enrichment analysis on selected genes",

        shiny::p(
          "Draw a lasso or box selection on the scatter plot above to pick a ",
          "set of genes, then run a topGO enrichment analysis against the ",
          "background of all genes tested in both contrasts."
        ),
        bslib::layout_columns(
          col_widths = c(4, 4, 4),
          shiny::selectInput(ns("fea_ontology"), "GO ontology",
                             choices = c("Biological Process" = "BP",
                                         "Molecular Function" = "MF",
                                         "Cellular Component" = "CC"),
                             selected = "BP"
          ),
          shiny::div(
            style = "padding-top: 1.9em;",
            shiny::textOutput(ns("fea_selected_count"))
          ),
          shiny::div(
            style = "padding-top: 1.5em;",
            shiny::actionButton(ns("run_fea"), "Run enrichment analysis")
          )
        ),
        shinycssloaders::withSpinner(
          DT::DTOutput(ns("fea_table"))
        )
      )
    ),
  )
}

#' @keywords internal
#' @param dde reactive DeeDeeExperiment from mod_input_server()
mod_scatter_server <- function(id, dde, show_symbols, species) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    scatter_src <- ns("scatter_plotly")

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
                     pthresh = input$pthresh,
                     show_symbols = show_symbols(),
                     species = species(),
                     source = scatter_src
      )
    })

    scatter_data <- shiny::reactive({
      shiny::req(input$s1, input$s2)
      sel1 <- match(input$s1, names(dde()@dea))
      sel2 <- match(input$s2, names(dde()@dea))
      shiny::req(sel1, sel2)

      .deedee_scatter_data(dde(),
                           select1 = sel1,
                           select2 = sel2,
                           pthresh = input$pthresh,
                           show_symbols = show_symbols(),
                           species = species()
      )
    })

    selected_genes <- shiny::reactive({
      ed <- plotly::event_data("plotly_selected", source = scatter_src)
      prep <- scatter_data()
      if (is.null(ed) || nrow(ed) == 0 || is.null(prep)) {
        return(character(0))
      }
      comp <- prep$comp
      unique(stats::na.omit(comp$symbol1[match(ed$key, comp$rowname)]))
    })

    output$fea_selected_count <- shiny::renderText({
      paste(length(selected_genes()), "gene currently selected")
    })

    fea_result <- shiny::reactiveVal(NULL)
    shiny::observeEvent(list(input$s1, input$s2, input$pthresh), {
      fea_result(NULL)
    })

    shiny::observeEvent(input$run_fea, {
      genes <- selected_genes()
      if (length(genes) < 2) {
        shiny::showNotification(
          "Please lasso- or box-select at least two genes on the plot first.",
          type = "warning"
        )
        return(invisible(NULL))
      }

      sel1 <- match(input$s1, names(dde()@dea))
      sel2 <- match(input$s2, names(dde()@dea))

      res <- tryCatch(
        deedee_fea(dde(),
                   select1 = sel1,
                   select2 = sel2,
                   genes = genes,
                   species = species(),
                   ontology = input$fea_ontology
        ),
        error = function(e) {
          shiny::showNotification(
            paste("Enrichment analysis failed:", conditionMessage(e)),
            type = "error", duration = NULL
          )
          NULL
        }
      )
      fea_result(res)
    })

    output$fea_table <- DT::renderDT({
      res <- fea_result()
      shiny::validate(shiny::need(
        !is.null(res),
        "No enrichment results yet. Select genes on the plot and click 'Run enrichment analysis'."
      ))
      res <- beautify_GO_table(res)
      DT::datatable(res,
                    escape = FALSE,
                    rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 10)
      )
    })

    output$scatter_plotly <- plotly::renderPlotly({
      shiny::req(plot_obj())
      plot_obj()
    })

    output$scatter_out <- shiny::renderUI({
      plotly::plotlyOutput(ns("scatter_plotly"))
    })
  })
}
