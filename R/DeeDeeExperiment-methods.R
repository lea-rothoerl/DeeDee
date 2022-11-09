#' @name DeeDeeExperiment-methods
#'
#' @title Methods for DeeDeeExperiment objects
#'
#' @aliases
#' dea
#' dea<-
#'
#'
#' @description
#' Todo
#'
#' @param x A \code{\link{DeeDeeExperiment}} object
#' @param value Replacement value for replacement methods.
#'
#' @return Return value varies depending on method, as described below.
NULL


#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("dea",
          signature = "DeeDeeExperiment",
          definition = function(x) {
            x@dea
          })

#' @rdname DeeDeeExperiment-methods
#' @export
setReplaceMethod("dea",
                 c("DeeDeeExperiment", "list"),
                 definition = function(x, value) {
                   x@dea <- value
                   validObject(x)
                   x
                 })



#' @name DeeDeeExperiment-misc
#'
#' @title Miscellaneous DeeDeeExperiment methods
#'
#' @description
#' Miscellaneous methods for the \code{\link{DeeDeeExperiment}} class and its
#' descendants that do not fit into any other documentation category such as,
#' for example, show methods.
#'
#' @param object a \code{\link{DeeDeeExperiment}} object
#'
#' @return Returns NULL
#'
NULL


#' @rdname DeeDeeExperiment-misc
#' @export
#'
setMethod("show",
          signature = signature(object = "DeeDeeExperiment"),
          definition = function(object) {

            callNextMethod()
            cat(
              "Including ", length(object@dea), " DE analyses:\n",
              paste(names(object@dea), collapse = ", "),
              sep=""
            )
          })
