#' @keywords internal
mod_input_ui <- function(id, collapse_upload = FALSE) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(8, 4),

    # LEFT
    bslib::card(
      bslib::accordion(
        open = !collapse_upload,
        bslib::accordion_panel(
          "",
          shiny::fileInput(
            ns("dde_file"),
            "Upload a DeeDeeExperiment (.rds)",
            accept = ".rds",
            placeholder = "No file selected"
          )
        )
      ),

      bslib::card_header("Overview of your DDE"),
      shiny::tableOutput(ns("inp_overview")),
      shiny::uiOutput(ns("inp_barplot"))
    ),

    # RIGHT
    bslib::card(
      bslib::card_header("Contrasts to Analyze"),
      shiny::uiOutput(ns("datasets")),

      bslib::input_switch(
        ns("show_symbols"),
        "Show gene symbols instead of Ensembl IDs",
        value = FALSE
      ),

      shiny::textOutput(ns("species_detected")),

      shiny::numericInput(
        ns("overview_pthresh"),
        "Adjusted p-value threshold (for the DE count overview)",
        value = 0.05, min = 0.01, max = 1, step = 0.01
      )
    )
  )
}

#' @keywords internal
#' @param dde_arg a DeeDeeExperiment passed in the call of deedee_app(),
#'   or NULL if the user will upload one instead
mod_input_server <- function(id, dde_arg = NULL, arg_label = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    dde_raw <- shiny::reactive({
      shiny::req(shiny::isTruthy(input$dde_file) || !is.null(dde_arg))

      if (!is.null(dde_arg)) {
        checkmate::assertClass(dde_arg, "DeeDeeExperiment")
        return(dde_arg)
      }
      obj <- readRDS(input$dde_file$datapath)
      checkmate::assertClass(obj, "DeeDeeExperiment")
      obj
    })

    source_label <- shiny::reactive({
      if (!is.null(dde_arg)) {
        arg_label
      } else if (shiny::isTruthy(input$dde_file)) {
        input$dde_file$name
      } else {
        NULL
      }
    })

    output$datasets <- shiny::renderUI({
      shiny::req(dde_raw())
      all_names <- DeeDeeExperiment::getDEANames(dde_raw())
      shiny::checkboxGroupInput(session$ns("select_contrasts"),
                                "Select contrasts to use",
                                choices = all_names,
                                selected = all_names
      )
    })

    dde <- shiny::reactive({
      shiny::req(dde_raw(), input$select_contrasts)

      all_names <- DeeDeeExperiment::getDEANames(dde_raw())
      to_drop <- setdiff(all_names, input$select_contrasts)

      if (length(to_drop) == 0) {
        dde_raw()
      } else {
        DeeDeeExperiment::removeDEA(dde_raw(), to_drop)
      }
    })

    code_data <- shiny::reactive({
      shiny::req(dde_raw(), input$select_contrasts)
      fmt <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")

      if (!is.null(dde_arg)) {
        data <- arg_label
        preamble <- character(0)
      } else {
        shiny::req(input$dde_file)
        data <- "dde"
        preamble <- paste0("dde <- readRDS(", fmt(input$dde_file$name), ")")
      }

      to_drop <- setdiff(DeeDeeExperiment::getDEANames(dde_raw()),
                         input$select_contrasts)
      if (length(to_drop) > 0) {
        preamble <- c(preamble, paste0(
          "dde <- DeeDeeExperiment::removeDEA(", data, ", ", fmt(to_drop), ")"
        ))
        data <- "dde"
      }

      list(data = data, preamble = preamble)
    })

    output$inp_overview <- shiny::renderTable({
      shiny::req(dde())
      deedee_overview(
        data = dde(),
        pthresh = input$overview_pthresh
      )
    })

    output$inp_barplot <- shiny::renderUI({
      shiny::req(dde())
      shiny::plotOutput(session$ns("overview_plot"), height = "350px")
    })

    output$overview_plot <- shiny::renderPlot({
      shiny::req(dde())
      deedee_overview_plot(
        data = dde(),
        pthresh = input$overview_pthresh
      )
    })

    species <- shiny::reactive({
      shiny::req(dde_raw())
      ids <- rownames(deedee_from_dde(dde_raw())[[1]])
      if (any(grepl("^ENSMUSG", ids))) {
        "Mus_musculus"
      } else if (any(grepl("^ENSG", ids))) {
        "Homo_sapiens"
      } else {
        NA_character_
      }
    })

    output$species_detected <- shiny::renderText({
      sp <- species()
      shiny::validate(shiny::need(!is.na(sp),
          "Could not detect species from gene IDs
          (expected ENSG... or ENSMUSG...)."
      ))
      paste("Detected species:", sp)
    })

    list(
      dde = dde,
      show_symbols = shiny::reactive(input$show_symbols),
      species = species,
      source_label = source_label,
      code_data = code_data
    )
  })
}

