#' DeeDee RRHO Plot
#'
#' @description
#'
#' @param data instance of the DeeDeeExperiment class
#' @param select1 index of the 1st contrast to be compared (default = 1)
#' @param select2 index of the 2nd contrast to be compared (default = 2)
#' @param alternative `"two.sided"` (default) for concordant and discordant
#'                     regions, or `"enrichment"` for a one-sided test
#' @param corr boolean, apply Benjamini-Yekutieli FDR correction to the map
#'           (default = FALSE)
#'
#' @return plottable RRHO object
#'
#' @examples
#'
#' data(dde_macrophage, package = "DeeDee")
#' deedee_rrho(dde_macrophage, select1 = 1, select2 = 2)
#'
#' @export
#'
deedee_rrho <- function(data,
                        select1 = 1,
                        select2 = 2,
                        alternative = "two.sided",
                        corr = FALSE) {

  # ----------------------------- argument check ------------------------------
  checkmate::assert(length(names(data@dea)) >= 2)
  data <- deedee_from_dde(data)
  checkmate::assert_number(select1, lower = 1, upper = length(data))
  checkmate::assert_number(select2, lower = 1, upper = length(data))
  checkmate::assert_choice(alternative, c("two.sided", "enrichment"))
  checkmate::assert_logical(corr)

  # ---------------------------- data preparation -----------------------------
  d1 <- data[[select1]]
  d2 <- data[[select2]]

  score <- function(d) {
    p <- pmax(d$pval, .Machine$double.xmin)
    -log10(p) * sign(d$logFC)
  }

  list1 <- data.frame(gene = rownames(d1), score = score(d1))
  list2 <- data.frame(gene = rownames(d2), score = score(d2))

  # ------------------- creation of the resulting RRHO object -----------------
  res <- RRHO::RRHO(list1, list2,
                    labels = c(names(data)[select1], names(data)[select2]),
                    alternative = alternative,
                    plots = FALSE,
                    BY = corr,
                    log10.ind = TRUE
  )

  # --------------------------------- return ----------------------------------
  return(res)
}
