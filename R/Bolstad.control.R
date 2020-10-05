#' Control Bolstad functions
#' @param plot if \code{TRUE} then draw a plot (for functions that actually have plots)
#' @param quiet if \code{TRUE} then suppress the function output
#' @param ... additional parameters
#'
#' @return an invisible list of options and their values
#' @export
Bolstad.control = function(plot = TRUE, quiet = FALSE, ...){
  invisible(list(
            plot = plot,
            quiet = quiet))
}