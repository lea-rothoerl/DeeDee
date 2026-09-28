#' @keywords internal
deedee_overview <- function(data, pthresh = 0.05) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_class(data, "DeeDeeExperiment")
  checkmate::assert_number(pthresh, lower = 0, upper = 1)

  # --------------------------------- summary ---------------------------------
  dea_info <- DeeDeeExperiment::getDEAInfo(data)
  dea_list <- deedee_from_dde(data)
  contrasts <- DeeDeeExperiment::getDEANames(data)

  dea_package <- vapply(contrasts, function(name) {
    pkg <- dea_info[[name]]$package
    if (is.null(pkg) || length(pkg) != 1 || is.na(pkg)) {
      NA_character_
    } else {
      as.character(pkg)
    }
  }, character(1))

  res <- data.frame(
    contrast = contrasts,
    package = unname(dea_package),
    genes = unname(vapply(dea_list, nrow, integer(1))),
    de_genes = unname(vapply(dea_list, function(d) {
      sum(d$pval < pthresh, na.rm = TRUE)
    }, integer(1))),
    up = unname(vapply(dea_list, function(d) {
      sum(d$pval < pthresh & d$logFC > 0, na.rm = TRUE)
    }, integer(1))),
    down = unname(vapply(dea_list, function(d) {
      sum(d$pval < pthresh & d$logFC < 0, na.rm = TRUE)
    }, integer(1))),
    stringsAsFactors = FALSE
  )

  # --------------------------------- return ----------------------------------
  return(res)
}

# --- PLOTTER ---
#' @keywords internal
deedee_overview_plot <- function(data, pthresh = 0.05) {
  overview <- deedee_overview(data, pthresh)

  long <- data.frame(
    contrast = rep(overview$contrast, times = 2),
    direction = rep(c("up", "down"), each = nrow(overview)),
    n = c(overview$up, overview$down),
    stringsAsFactors = FALSE)

  long$direction <- factor(long$direction, levels = c("down", "up"))
  long$contrast <- factor(long$contrast, levels = rev(overview$contrast))

  col_up <- viridis::viridis(n = 1, begin = 0.4, option = "magma")
  col_down <- viridis::viridis(n = 1, begin = 0.9, option = "magma")

  ggplot2::ggplot(long, ggplot2::aes(x = contrast, y = n, fill = direction )) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::scale_fill_manual(values = c(up = col_up, down = col_down)) +
    ggplot2::labs(x = NULL, y = paste0("DE genes (p < ", pthresh, ")"), fill = NULL) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
