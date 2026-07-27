#' @keywords internal
mod_input_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(8, 4),
    bslib::card(
      bslib::card_header("Load data"),
      shiny::fileInput(ns("dde_file"), "Upload a DeeDeeExperiment (.rds)",
                       accept = ".rds", placeholder = "No file selected"
      ),
      shiny::tableOutput(ns("inp_infobox"))
    ),
    bslib::card(
      bslib::card_header("Contrasts"),
      shiny::uiOutput(ns("datasets"))
    )
  )
}

#' @keywords internal
#' @param dde_arg a DeeDeeExperiment passed in the call of deedee_app(),
#'   or NULL if the user will upload one instead
mod_input_server <- function(id, dde_arg = NULL) {
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

    output$datasets <- shiny::renderUI({
      shiny::req(dde_raw())
      all_names <- DeeDeeExperiment::getDEANames(dde_raw())
      shiny::checkboxGroupInput(session$ns("select_contrasts"),
                                "Select contrasts to use",
                                choices  = all_names,
                                selected = all_names
      )
    })

    shiny::reactive({
      shiny::req(dde_raw(), input$select_contrasts)

      all_names <- DeeDeeExperiment::getDEANames(dde_raw())
      to_drop   <- setdiff(all_names, input$select_contrasts)

      if (length(to_drop) == 0) {
        dde_raw()
      } else {
        DeeDeeExperiment::removeDEA(dde_raw(), to_drop)
      }
    })
  })
}
