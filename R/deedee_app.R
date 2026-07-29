deedee_app <- function(deedee_obj = NULL) {

  ui <- bslib::page_navbar(
    title = "DeeDee",
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    bslib::nav_panel("Input", mod_input_ui("input")),
    bslib::nav_menu("TODO groups",
#                    bslib::nav_panel("Venn", mod_venn_ui("venn")),
                    bslib::nav_panel("UpSet", mod_upset_ui("upset")),
                    bslib::nav_panel("CAT", mod_cat_ui("cat")),
                    bslib::nav_panel("RRHO", mod_rrho_ui("rrho")),
                    bslib::nav_panel("Scatter", mod_scatter_ui("scatter")),
                    bslib::nav_panel("QQ", mod_qq_ui("qq")),
                    bslib::nav_panel("Heatmap", mod_heatmap_ui("heatmap")),
    )
  )

  server <- function(input, output, session) {
    dde <- mod_input_server("input", dde_arg = deedee_obj)

#    mod_venn_server("venn", dde)
    mod_upset_server("upset", dde)
    mod_cat_server("cat", dde)
    mod_rrho_server("rrho", dde)
    mod_scatter_server("scatter", dde)
    mod_qq_server("qq", dde)
    mod_heatmap_server("heatmap", dde)
  }

  shiny::shinyApp(ui, server)
}
