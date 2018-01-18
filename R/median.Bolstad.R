#' Median generic
#' 
#' @param x an object.
#' @param na.rm Ideally if \code{TRUE} then missing values will be removed, but not currently used.
#' @param ... [>=R3.4.0 only] Not currently used.
#' @details If \code{x} is an object of class \code{Bolstad} then the posterior 
#'   median of the parameter of interest will be calculated.
#' @author James Curran
#' @method median Bolstad
if(is.na(match("...", names(formals(median))))) {
  median.Bolstad = function(x, na.rm = FALSE) {
    return(quantile(x, probs = 0.5))
  }
}else{
  median.Bolstad = function(x, na.rm = FALSE, ...) {
    return(quantile(x, probs = 0.5))
  }
}

