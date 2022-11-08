#' @name DeeDeeExperiment-dea
#'
#' @title Methods for deedee exps
#'
#' @aliases
#' dea
#'
#' @description
#' Todo
#'
#' @param x deex
#'
#' @return Return value varies depending on method, as described below.
NULL


#' @rdname DeeDeeExperiment-dea
#' @export
setMethod("dea",
          signature = "DeeDeeExperiment",
          definition = function(x) {
            x@dea
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
