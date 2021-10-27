DeeDeeObject <- setClass(

  "DeeDeeObject",
  slots = c(
    DeeDeeList = "list"
  ),
  prototype = list(
    DeeDeeList = NULL
  ),
  validity = function(object) {
    res <- object@DeeDeeList
    if (class(res) != "list") {
      return("DeeDeeList must be a list of DeeDee data.frames.")
    }
    for (i in 1:length(res)) {
      if (class(res[[i]]) != "data.frame") {
        return("DeeDeeList must be a list of DeeDee data.frames.")
      }
      if (checkmate::test_subset(names(res[[i]]), c("logFC", "pval")) == FALSE) {
        return("DeeDeeList must be a list of DeeDee data.frames.")
      }
      if (class(res[[i]]$logFC) != "numeric") {
        return("LogFC values must be a vector of type 'numeric'.")
      }
      if (class(res[[i]]$pval) != "numeric") {
        return("P-values must be a vector of type 'numeric'.")
      }
      if (max(res[[i]]$pval) > 1) {
        return("A p-value cannot be greater than 1.")
      }
      if (min(res[[i]]$pval) < 0) {
        return("A p-value cannot be lower than 0.")
      }
    }
    return(TRUE)
  }
)
