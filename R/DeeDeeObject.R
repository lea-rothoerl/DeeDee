DeeDeeObject <- setClass(

  "DeeDeeObject",

  slots = c(
    DEA_results = "list"
  ),

  prototype = list(
    DEA_results = NULL
  ),

  validity = function(object) {
    res <- object@DEA_results
    if (class(res) != "list") {
      return("DEA_results must be a list of DeeDee data.frames.")
    }
    for (i in 1:length(res)) {
      if (class(res[[i]]) != "data.frame") {
        return("DEA_results must be a list of DeeDee data.frames.")
      }
      if (checkmate::test_subset(names(res[[i]]), c("logFC", "pval")) == FALSE) {
        return("DEA_results must be a list of DeeDee data.frames.")
      }
      if (class(res[[i]]$logFC) != "numeric") {
        return("LogFC values must be of type 'numeric'.")
      }
      if (class(res[[i]]$pval) != "numeric") {
        return("P-values must be of type 'numeric'.")
      }
      if (max(res[[i]]$pval) > 1) {
        return("P-value cannot be greater than 1.")
      }
      if (min(res[[i]]$pval) < 0) {
        return("P-value cannot be smaller than 0.")
      }
    }
    return(TRUE)
  }
)
