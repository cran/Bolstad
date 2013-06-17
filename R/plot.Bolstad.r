plot.Bolstad = function(x, overlay = TRUE, which = c(1,3), densCols = c("red","green","blue")[which],  ...){
  
  which = sort(which)
  
  if(is.null(which) || length(which) <= 0 || length(which) > 3 || any(!grepl('^[1-3]+$', which))){
    stop("parameter which can only take vectors of length 3 containing the values 1, 2 and 3")
  }
  
  if(overlay){
    with(x,{
      Y = cbind(prior, likelihood, posterior)[,which]
      ylim = c(0, 1.1 * max(Y))
      
      plot(param.x, Y[,1], ylim = ylim, type="l",
           lty = (3:1)[which[1]], col = densCols[1],
           xlab = eval(expression(name)), ylab = "",
           main = "Shape of prior and posterior", ...)
      
      for(i in 2:ncol(Y)){
        lines(param.x, Y[,i], lty = (3:1)[which[i]], col = densCols[i])
      }
      
      legend("topleft", lty = (3:1)[which], col = densCols, legend = c("Prior", "Likelihood", "Posterior")[which], bty = 'n')
    })
  }else{
    oldpar = par(mfrow = c(1, length(which)), mai = c(0.7, 0.1, 0.2, 0.1), yaxs = 'i', xaxs = 'i')
    
    with(x,{
      Y = cbind(prior, likelihood, posterior)[,which]
      ylim = c(0, 1.1 * max(Y))
      legend = c("Prior", "Likelihood", "Posterior")[which]
      
      plot(param.x, Y[,1], ylim = ylim, type="l",
           col = densCols[1],
           xlab = eval(expression(name)), ylab = "",
           main = legend[1], axes = FALSE, ...)
      axis(1)
      box()
      
      for(i in 2:ncol(Y)){
        plot(param.x, Y[,i], ylim = ylim, col = densCols[i], type = 'l', xlab = "", 
             main = legend[i], axes = FALSE, ...)
        box()
      }
      
      par(oldpar)
    })
  }
}