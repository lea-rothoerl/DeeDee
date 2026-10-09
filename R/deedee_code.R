#' @keywords internal
.deedee_code_text <- function(fn, data, args, preamble = character(0),
                              comments = list(), show = "p") {
  fmt <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")

  arg_lines <- vapply(seq_along(args), function(i) {
    nm <- names(args)[i]
    cm <- comments[[nm]]
    paste0(
      "  ", nm, " = ", fmt(args[[i]]),
      if (i < length(args)) "," else "",
      if (!is.null(cm)) paste0("  # ", cm) else ""
    )
  }, character(1))

  paste(
    c(
      "library(DeeDee)",
      "",
      if (length(preamble) > 0) c(preamble, ""),
      paste0("p <- ", fn, "("),
      paste0("  ", data, ","),
      arg_lines,
      ")",
      show
    ),
    collapse = "\n"
  )
}

#' @keywords internal
.deedee_code_server <- function(input, session, code_fn, id_prefix = "") {
  shiny::observeEvent(input[[paste0(id_prefix, "show_code")]], {
    el_id <- session$ns(paste0(id_prefix, "code_text"))

    shiny::showModal(shiny::modalDialog(
      title = "Code to reproduce this plot",
      shiny::tags$pre(
        id = el_id,
        style = "white-space: pre-wrap; word-break: break-word;",
        code_fn()
      ),
      shiny::tags$button(
        "Copy to clipboard",
        type = "button",
        class = "btn btn-outline-secondary",
        onclick = sprintf("copyDeeDeeCode('%s', this)", el_id)
      ),
      easyClose = TRUE,
      footer = shiny::modalButton("Close"),
      size = "l"
    ))
  })
}
