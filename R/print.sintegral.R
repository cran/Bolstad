#' Generic print method
#' 
#' Print the value of an \code{sintegral} object, specifically the 
#' (approximate) value of the integral.
#' 
#' @param x An object of type \code{sintegral}--see \code{\link{sintegral}}.
#' @param \dots other parameters passed to \code{paste0}.
#' @export
print.sintegral = function(x, ...){
  cat(paste0("Value: ", x$value, "\n", ...))
}