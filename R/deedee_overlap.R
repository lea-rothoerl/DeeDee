#' DeeDee Overlap Gene List
#'
#' @description `deedee_overlap` computes the differentially expressed genes
#' shared across all selected contrasts, filtered by p-value threshold and
#' regulation direction.
#'
#' @param data instance of the DeeDeeExperiment class
#' @param mode show all overlapping DE genes (`both`, default),
#'             only conjointly up-regulated (`up`)
#'             or only conjointly down-regulated (`down`) genes
#' @param pthresh threshold for p-values to be in-/excluded (default = 0.05)
#' @param select vector of DEA slot indexes to use for plotting, must be at least 2
#'
#' @return data.frame with one row per overlapping gene
#'
#' @examples
#'
#' data(dde_macrophage, package = "DeeDee")
#' deedee_overlap(dde_macrophage, mode = "both", pthresh = 0.05)
#'
#' @export
#'

deedee_overlap <- function(data,
                        mode = "both",
                        pthresh = 0.05,
                        select = NULL) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert(length(names(data@dea)) >= 2)
  data <- deedee_from_dde(data)
  checkmate::assert_number(pthresh, lower = 0, upper = 1)
  choices <- c("up", "down", "both")
  checkmate::assert_choice(mode, choices)

  n_contrasts <- length(data)

  # --------------------------- selection handling ----------------------------
  if (is.null(select)) {
    select <- seq_len(n_contrasts)
  } else {
    checkmate::assert_integerish(select, lower = 1)
    checkmate::assert_true(length(select) >= 2)
    checkmate::assert_true(length(unique(select)) == length(select))
    checkmate::assert_true(all(select <= n_contrasts))
  }

  data <- data[select]

  # ---------------------------- data preparation -----------------------------
  for (i in 1:length(data)) {
    data[i][[1]] <- subset(
      data[i][[1]],
      data[i][[1]]$pval < pthresh
    )

    if (length(data[i][[1]][[1]]) == 0) {
      return(NULL)
    }

    if (mode == "up") {
      data[i][[1]] <- subset(
        data[i][[1]],
        data[i][[1]]$logFC > 0
      )
    }
    if (mode == "down") {
      data[i][[1]] <- subset(
        data[i][[1]],
        data[i][[1]]$logFC < 0
      )
    }

    common_genes <- Reduce(intersect, lapply(data, rownames))

    if (length(common_genes) == 0) {
      return(NULL)
    }
  }

  # --------------------- creation of the resulting list ----------------------
  res <- data.frame(gene = common_genes, stringsAsFactors = FALSE)
  for (i in seq_along(data)) {
    nm <- names(data)[i]
    res[[paste0(nm, "_logFC")]] <- data[[i]][common_genes, "logFC"]
    res[[paste0(nm, "_pval")]]  <- data[[i]][common_genes, "pval"]
  }

  # --------------------------------- return ----------------------------------
  return(res)
}
