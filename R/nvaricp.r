#' Bayesian inference for a normal standard deviation with a scaled inverse
#' chi-squared distribution
#' 
#' Evaluates and plots the posterior density for \eqn{\sigma}{sigma}, the
#' standard deviation of a Normal distribution where the mean \eqn{\mu}{mu} is
#' known
#' 
#' 
#' @param y a random sample from a
#' \eqn{normal(\mu,\sigma^2)}{normal(mu,sigma^2)} distribution.
#' @param mu the known population mean of the random sample.
#' @param S0 the prior scaling factor.
#' @param kappa the degrees of freedom of the prior.
#' @param cred.int if TRUE then a 100(1-alpha) percent credible interval will
#' be calculated for \eqn{\sigma}{sigma}
#' @param alpha controls the width of the credible interval. Ignored if
#' cred.int is FALSE
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @return A list will be returned with the following components:
#' 
#' \item{sigma}{the vaules of \eqn{\sigma}{sigma} for which the prior,
#' likelihood and posterior have been calculated} \item{prior}{the prior
#' density for \eqn{\sigma}{sigma}} \item{likelihood}{the likelihood function
#' for \eqn{\sigma}{sigma} given \eqn{y}{y}} \item{posterior}{the posterior
#' density of \eqn{\mu}{sigma} given \eqn{y}{y}} \item{S1}{the posterior
#' scaling constant} \item{kappa1}{the posterior degrees of freedom}
#' @keywords misc
#' @examples
#' 
#' ## Suppose we have five observations from a normal(mu, sigma^2)
#' ## distribution mu = 200 which are 206.4, 197.4, 212.7, 208.5.
#' y = c(206.4, 197.4, 212.7, 208.5, 203.4)
#' 
#' ## We wish to choose a prior that has a median of 8. This happens when
#' ## S0 = 29.11 and kappa = 1
#' nvaricp(y,200,29.11,1)
#' 
#' ##  Same as the previous example but a calculate a 95% credible
#' ## interval for sigma
#' nvaricp(y,200,29.11,1,cred.int=TRUE)
#' 
#' ##  Same as the previous example but a calculate a 95% credible
#' ## interval for sigma by hand. Note that the syntax of sintegral has
#' ## changed
#' results = nvaricp(y,200,29.11,1,cred.int=TRUE)
#' attach(results)
#' cdf = sintegral(sigma,posterior)$cdf
#' Finv = approxfun(cdf$y,cdf$x)
#' lb = Finv(0.025)
#' ub = Finv(0.975)
#' cat(paste("95% credible interval for sigma: [",
#'               signif(lb,4),", ", signif(ub,4),"]\n",sep=""))
#' 
#' ##  Same as the previous example but a calculate a 95% credible
#' ## interval for sigma by hand using the new quantile function
#' results = nvaricp(y, 200, 29.11, 1)
#' ci = quantile(results, c(0.025, 0.975))
#' lb = ci[1]
#' ub = ci[2]
#' cat(paste("95% credible interval for sigma: [",
#'               signif(lb,4),", ", signif(ub,4),"]\n",sep=""))
#' 
#' 
#' @export nvaricp
nvaricp = function(y, mu, S0, kappa, cred.int = FALSE, alpha = 0.05, plot = TRUE){
  n = length(y)
  SST = sum((y-mu)^2)

  if(kappa > 0){
    S1 = S0 + SST
    kappa1 = kappa + n

    k1 = qchisq(0.01, kappa)
    k2 = S0/k1
    k3 = sqrt(k2)


    sigma = seq(0,k3,length = 1002)[-1]
    k4 = diff(sigma)[1]
    sigma.sq = sigma^2
    log.prior =  -((kappa-1)/2+1)*log(sigma.sq)-S0/(2*sigma.sq)
    prior = exp(log.prior)

    log.like =  -(n/2)*log(sigma.sq)-SST/(2*sigma.sq)
    likelihood =  exp(log.like)


    kint = ((2*sum(prior))-prior[1]-prior[1001])*k4/2*0.99
    prior = prior/kint

    posterior = prior*likelihood
    kint = ((2*sum(posterior))-posterior[1]-posterior[1001])*k4/2
    posterior = posterior/kint

    y.max = max(c(prior,posterior))
    k1 = qchisq(0.01,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)

    plot(sigma, prior, type = "l", col = "blue", ylim = c(0, 1.1 * y.max),
         xlim = c(0,k3),
         main = expression(
             paste("Shape of Inverse ", chi^2," and posterior for ",
                   sigma, sep = "")),
         xlab = expression(sigma),
         ylab = "Density")
    lines(sigma, posterior, lty = 2, col = "red")
  }else if(kappa == 0){                   ## Jeffrey's prior
    S = 0
    S1 = S+SST
    kappa1 = kappa+n

    k1 = qchisq(0.001,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)
    k4 = k3/1000

    sigma = seq(0,k3,length = 1002)[-1]
    sigma.sq = sigma^2
    likelihood = NULL

    log.posterior =  -((kappa1-1)/2+1)*log(sigma.sq)-S1/(2*sigma.sq)
    posterior = exp(log.posterior)
    kint = ((2*sum(posterior))-posterior[1]-posterior[1001])*k4/(2*.999)
    posterior = posterior/kint

    log.prior =  -((kappa-1)/2+1)*log(sigma.sq)-S0/(2*sigma.sq)
    prior = exp(log.prior)
    kint = ((2*sum(prior))-prior[1]-prior[1001])*k4/2
    prior = prior/kint

    k1 = qchisq(0.01,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)
    k4 = 1.2*max(posterior)

    if(plot){
      plot(sigma,prior,type = "l",col = "blue", ylim = c(0,k4),
           main = expression(paste("Shape of prior and posterior for ", sigma,
               sep = "")),
           xlab = expression(sigma),ylab = "Density")
      lines(sigma,posterior,lty = 2,col = "red")
      }
  }else if(kappa<0){
    S0 = 0
    S1 = S0+SST
    kappa1 = kappa+n

    k1 = qchisq(0.001,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)
    k4 = k3/1000

    sigma = seq(0,k3,length = 1002)[-1]
    sigma.sq = sigma^2

    log.posterior =  -((kappa1-1)/2+1)*log(sigma.sq)-S1/(2*sigma.sq)
    posterior = exp(log.posterior)
    kint = ((2*sum(posterior))-posterior[1]-posterior[1001])*k4/(2*.999)
    posterior = posterior/kint

    log.prior =  -((kappa-1)/2+1)*log(sigma.sq)-S0/(2*sigma.sq)
    prior = exp(log.prior)
    kint = ((2*sum(prior))-prior[1]-prior[1001])*k4/2
    prior = prior/kint

    likelihood = NULL

    k1 = qchisq(0.01,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)
    k4 = 1.2*max(posterior)

    if(plot){
      plot(sigma,prior,type = "l",col = "blue", xlim = c(0,k3),ylim = c(0,k4),
           main = expression(paste("Shape of prior and posterior for ", sigma,
               sep = "")),
           xlab = expression(sigma),ylab = "Density")
      lines(sigma,posterior,lty = 2,col = "red")
      legend(sigma[1],0.9*k4,lty = 1:2,col = c("blue","red")
             ,legend = c("Prior","Posterior"))
      }
  }

  cat(paste("S1: ",signif(S1,4)," kappa1 :", signif(kappa1,3),"\n",sep = ""))

  if(cred.int){
    if(kappa1<2)
      cat("Unable to calculate credible interval for sigma if kappa1<= 2\n")
    else{
      sigmahat.post.mean = sqrt(S1/(kappa1-2))
      cat(paste("Estimate of sigma using posterior mean: ",
                signif(sigmahat.post.mean,4),"\n",sep = ""))
    }

    q50 = qchisq(0.5,kappa1)
    sigmahat.post.median = sqrt(S1/q50)
    cat(paste("Estimate of sigma using posterior median: ",
              signif(sigmahat.post.median,4),"\n",sep = ""))

    cdf = sintegral(sigma,posterior)$cdf
    Finv = approxfun(cdf$y,cdf$x)
    lb = Finv(alpha/2)
    ub = Finv(1-alpha/2)
    cat(paste(round(100*(1-alpha)),"% credible interval for sigma: [",
              signif(lb,4),", ", signif(ub,4),"]\n",sep = ""))
    if(plot)
      abline(v = c(lb,ub),col = "blue",lty = 3)

  }


  results = list(param.x = sigma, prior = prior, likelihood = likelihood,
                 posterior = posterior, 
                 sigma = sigma, # for backwards compat. only
                 S1 = S1, kappa1 = kappa1)
  class(results) = 'Bolstad'
  invisible(results)
}

