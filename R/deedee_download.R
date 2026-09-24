#' @keywords internal
.deedee_download_ui <- function(ns, label = "Download plot (PNG)", id_prefix = "") {
  shiny::tagList(
    shiny::downloadButton(ns(paste0(id_prefix, "download_plot")), label),
    bslib::input_switch(
      ns(paste0(id_prefix, "show_footer")),
      "Include data source in downloaded image",
      value = TRUE
    )
  )
}

#' @keywords internal
.deedee_download_server <- function(input, output, id_prefix = "",
                                    filename_prefix, draw_fn, source_label,
                                    width = 8, height = 6) {
  output[[paste0(id_prefix, "download_plot")]] <- shiny::downloadHandler(
    filename = function() {
      paste0(filename_prefix, "_", format(Sys.Date(), "%Y%m%d"), ".png")
    },
    content = function(file) {
      shiny::isolate({
        grDevices::png(file, width = width, height = height, units = "in", res = 300)
        on.exit(grDevices::dev.off(), add = TRUE)

        draw_fn()

        show_footer <- isTRUE(input[[paste0(id_prefix, "show_footer")]])
        lbl <- source_label()
        if (show_footer && checkmate::test_string(lbl, min.chars = 1)) {
          grid::upViewport(0)
          grid::grid.text(
            label = paste("Source:", lbl),
            x = grid::unit(1, "npc") - grid::unit(0.1, "in"),
            y = grid::unit(0.1, "in"),
            just = c("right", "bottom"),
            gp = grid::gpar(fontsize = 8, col = "grey40")
          )
        }
      })
    }
  )
}

#' @keywords internal
.deedee_download_table_server <- function(output, filename_prefix, table_fn) {
  output$download_table <- shiny::downloadHandler(
    filename = function() {
      paste0(filename_prefix, "_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
    },
    content = function(file) {
      tbl <- shiny::isolate(table_fn())
      shiny::req(tbl)
      writexl::write_xlsx(tbl, path = file)
    }
  )
}
