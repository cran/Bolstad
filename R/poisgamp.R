#' Poisson sampling with a gamma prior
#' 
#' Evaluates and plots the posterior density for \eqn{\mu}{mu}, the mean rate
#' of occurance in a Poisson process and a \eqn{gamma} prior on \eqn{\mu}{mu}
#' 
#' 
#' @param y a random sample from a Poisson distribution.
#' @param shape the shape parameter of the \eqn{gamma} prior.
#' @param rate the rate parameter of the \eqn{gamma} prior. Note that the scale
#' is \eqn{1 / rate}
#' @param scale the scale parameter of the \eqn{gamma} prior
#' @param alpha the width of the credible interval is controlled by the
#' parameter alpha.
#' @param \dots additional arguments that are passed to \code{Bolstad.control}
#' @return An object of class 'Bolstad' is returned. This is a list with the
#' following components:
#' 
#' \item{prior}{the prior density assigned to \eqn{\mu}{mu}}
#' \item{likelihood}{the scaled likelihood function for \eqn{\mu}{mu} given
#' \eqn{y}} \item{posterior}{the posterior probability of \eqn{\mu}{mu} given
#' \eqn{y}} \item{shape}{the shape parameter for the \eqn{gamma} posterior}
#' \item{rate}{the rate parameter for the \eqn{gamma} posterior}
#' @seealso \code{\link{poisdp}} \code{\link{poisgcp}}
#' @keywords misc
#' @examples
#' 
#' ## simplest call with an observation of 4 and a gamma(1, 1), i.e. an exponential prior on the
#' ## mu
#' poisgamp(4, 1, 1)
#' 
#' ##  Same as the previous example but a gamma(10, ) prior
#' poisgamp(4, 10, 1)
#' 
#' ##  Same as the previous example but an improper gamma(1, ) prior
#' poisgamp(4, 1, 0)
#' 
#' ## A random sample of 50 observations from a Poisson distribution with
#' ## parameter mu = 3 and  gamma(6,3) prior
#' set.seed(123)
#' y = rpois(50,3)
#' poisgamp(y,6,3)
#' 
#' ## In this example we have a random sample from a Poisson distribution
#' ## with an unknown mean. We will use a gamma(6,3) prior to obtain the
#' ## posterior gamma distribution, and use the R function qgamma to get a
#' ## 95% credible interval for mu
#' y = c(3,4,4,3,3,4,2,3,1,7)
#' results = poisgamp(y,6,3)
#' ci = qgamma(c(0.025,0.975),results$shape, results$rate)
#' cat(paste("95% credible interval for mu: [",round(ci[1],3), ",", round(ci[2],3)),"]\n")
#' 
#' ## In this example we have a random sample from a Poisson distribution
#' ## with an unknown mean. We will use a gamma(6,3) prior to obtain the
#' ## posterior gamma distribution, and use the R function qgamma to get a
#' ## 95% credible interval for mu
#' y = c(3,4,4,3,3,4,2,3,1,7)
#' results = poisgamp(y, 6, 3)
#' ci = quantile(results, c(0.025, 0.975))
#' cat(paste("95% credible interval for mu: [",round(ci[1],3), ",", round(ci[2],3)),"]\n")
#' 
#' 
#' @export poisgamp
poisgamp = function(y, shape, rate = 1, scale = 1 / rate, 
                    alpha = 0.05,  ...){
  n = length(y)
  y.sum = sum(y)
  
  if(is.null(y) || length(y)==0)
    stop("Error: y has no data")
  
  if(any(y < 0))
    stop("Error: y contains negative values")
  
  if(scale !=1 & rate == 1){
    rate = 1 / scale
  }
  
  if(shape < 0 || rate < 0)
    stop("Shape parameter and rate parameter must be greater than or equal to zero")
  
  quiet = Bolstad.control(...)$quiet
  if(!quiet){
    cat("Summary statistics for data\n")
    cat("---------------------------\n")
    cat(paste("Number of observations:\t", n, "\n"))
    cat(paste("Sum of observations:\t", y.sum, "\n\n"))
  }
  
  if(rate > 0){                              ##proper gamma prior
    upperBnd = qgamma(0.9999, shape, rate)
    stepSize = upperBnd / 1000
    mu = seq(0, upperBnd, by = stepSize)
    shapePost = shape + y.sum
    ratePost = rate + n
    
    prior = dgamma(mu, shape, rate)
    posterior = dgamma(mu, shapePost, ratePost)
  }else if(rate == 0){
    shapePost = shape + y.sum
    ratePost = rate + n
    upperBnd = qgamma(0.9999, shapePost, ratePost)
    stepSize = upperBnd / 1000
    mu = seq(0, upperBnd, by = stepSize)
    
    ##mu[1] = mu[2] ## fixes infinite upper bound problem
    
    prior = mu^(shape - 1)
    priorInt = integrate(function(x)x^(shape - 1), mu[1], mu[1001])$value
    prior = prior / priorInt
    
  }else{
    stop("Error: rate must be greater or equal to zero")
  }

  likelihood = sapply(mu, dpois, x = y, simplify = "array")
  if(is.matrix(likelihood)){
    likelihood = apply(likelihood, 2, prod)
  }
  
  posterior = dgamma(mu, shapePost, ratePost)
  credInt = qgamma(c(alpha * 0.5 , 1 - alpha * 0.5), shapePost, ratePost)
  
  if(!quiet){
    cat("Summary statistics for posterior\n")
    cat("--------------------------------\n")
    cat(paste("Shape parameter (r):\t", shapePost, "\n"))
    cat(paste("Rate parameter (v):\t", ratePost, "\n"))
    cat(sprintf("%d%% credible interval for mu:\t[%.2f, %.2f]\n",
                round(100 * (1 - alpha)),
                credInt[1], 
                credInt[2]))
  }
  
  if(Bolstad.control(...)$plot){
    y.max = max(prior[is.finite(prior)], posterior)
    plot(mu[is.finite(prior)], prior[is.finite(prior)], 
         ylim = c(0, 1.1 * y.max), xlab = expression(mu),          
         ylab = "Density",          main = "Shape of gamma prior and posterior\n for Poisson mean",          
         type = "l", 
         lty = 2, 
         col = "red")
    lines(mu, posterior, lty = 3, col = "blue")
    legend("topleft", bty = "n", lty = 2:3, 
           col = c("red", "blue"),
           legend = c("Prior", "Posterior"), 
           cex = 0.7)
  }
 
  
  results = list(name = 'mu', 
                 param.x = mu, 
                 prior = prior, 
                 likelihood = likelihood, 
                 posterior = posterior,
                 mean = shapePost / ratePost, 
                 var = shapePost / ratePost^2,
                 cdf = function(m, ...){pgamma(m, shape = shapePost, rate = ratePost, ...)},
                 quantileFun = function(probs, ...){qgamma(probs, shape = shapePost, rate = ratePost, ...)},
                 mu = mu, # for backwards compatibility only
                 shape = shapePost, 
                 rate = ratePost)
  
  class(results) = 'Bolstad'
  invisible(results)
}



