poisgamp = function(y, r, v, plot = TRUE, suppressOutput = FALSE){
  n = length(y)
  y.sum = sum(y)

  if(is.null(y) || length(y)==0)
    stop("Error: y has no data")

  if(any(y < 0))
    stop("Error: y contains negative values")

  if(r<0 || v<0)
    stop("Shape parameter r and rate parameter v must be greater than or equal to zero")

  if(!suppressOutput){
    cat("Summary statistics for data\n")
    cat("---------------------------\n")
    cat(paste("Number of observations:\t", n,"\n"))
    cat(paste("Sum of observations:\t", y.sum,"\n\n"))
  }
  
  if(v>0){                              ##proper gamma prior
    v.inv = 1/v
    k1 = qgamma(0.9999,r,v)
    k2 = k1/1000
    mu = seq(0,k1, by = k2)
    r1 = r+y.sum
    v1 = v+n
    v1.inv = 1/v1

    prior = dgamma(mu,r, v)
    likelihood = matrix(0, ncol = length(mu), nrow = length(y))
    for(i in 1:length(mu)){
        likelihood[,i] = dpois(y,mu[i])
    }
    likelihood = apply(likelihood, 2, prod)
    posterior = dgamma(mu, r1, v1)

    k3 = qgamma(c(0.005,0.995),r1,v1)

    if(!suppressOutput){
      cat("Summary statistics for posterior\n")
      cat("--------------------------------\n")
      cat(paste("Shape parameter r:\t", r1,"\n"))
      cat(paste("Rate parameter v:\t",v1,"\n"))
      cat(paste("99% credible interval for mu:\t[",round(k3[1],2), ",",round(k3[2],2), "]\n"))
    }
    
    if(plot){
      y.max = max(prior,posterior)
      plot(mu,prior,ylim = c(0,1.1*y.max),xlab = expression(mu)
           ,ylab = "Density",
           ,main = "Shape of gamma prior and posterior for Poisson mean"
           ,type = "l",lty = 2,col = "red")
      lines(mu,posterior,lty = 3,col = "blue")
      legend(mu[1],y.max,lty = 2:3,col = c("red","blue"),legend = c("Prior","Posterior"))
    }
  }else if(v == 0){
    r1 = r+y.sum
    v1 = v+n
    v1.inv = 1/v1
    k3 = qgamma(c(0.005,0.995),r1,v1)
    k4 = k3[2]/1000
    mu = seq(0,k3[2],by = k4)

    mu[1] = mu[2] ## fixes infinite upper bound problem

    if(!suppressOutput){
      cat("Summary statistics for posterior\n")
      cat("--------------------------------\n")
      cat(paste("Shape parameter r:\t", r1,"\n"))
      cat(paste("Rate parameter v:\t",v1,"\n"))
      cat(paste("99% credible interval :\t[",round(k3[1],2),", ",round(k3[2],2), "]\n"))
    }
    
    prior = mu^(r-1)
    kint = (2*sum(prior)-prior[1001])*k4/2
    prior = prior/kint

    likelihood = matrix(0, ncol = length(mu), nrow = length(y))
    for(i in 1:length(mu)){
        likelihood[,i] = dpois(y,mu[i])
    }
    likelihood = apply(likelihood, 2, prod)

    posterior = dgamma(mu, r1, v1)
    
    if(plot){
      y.max = max(prior,posterior)
      plot(mu,prior,ylim = c(0,1.1*y.max),xlab = expression(mu)
           ,ylab = "Density",
           ,main = "Shape of gamma prior and posterior for Poisson mean"
           ,type = "l",lty = 2,col = "red")
      lines(mu,posterior,lty = 3,col = "blue")
      legend(mu[1],y.max,lty = 2:3,col = c("red","blue"),legend = c("Prior","Posterior"))
    }
  }else{
    stop("Error: v must be greater or equal to zero")
  }

  results = list(name = 'mu', param.x = mu, prior = prior, likelihood = likelihood, posterior = posterior,
                 mu = mu, # for backwards compatibility only
                 r = r1,v = v1)
  
  class(results) = 'Bolstad'
  invisible(results)
}



