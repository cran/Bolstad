normnp<-function(x,sigma.x,m.x=0,s.x=1,n.mu = 100,ret=FALSE){

	# x - the vector of observations
	# sigma.x - the population standard deviation
	# m.x - the mean of the normal prior
	# s.x - the standard deviation of the normal prior
	# ret - if true then the prior, likelihood, posterior, mean, variance, and
	# quantiles are returned as a list

	mean.x<-mean(x)
	
	if(n.mu<100)
	{
		warning("Number of prior values of mu must be greater than 100")		
		n.mu<-100
	}

	if(s.x<=0)
		stop("The std. deviation of the prior must be greater than zero")

	lb<-m.x-3.5*s.x
	ub<-m.x+3.5*s.x

	mu<-seq(lb,ub,length=n.mu)
	mu.prior<-dnorm(mu,m.x,s.x)
	
	n.x<-length(x)
	
	likelihood<-exp(-n.x/(2*sigma.x^2)*(mean.x-mu)^2)

	precision<-1/s.x^2
	post.precision<-precision+(n.x/sigma.x^2)
	post.sd<-sqrt(1/post.precision)
	post.mean<-(precision/post.precision*m.x)+((n.x/sigma.x^2)/post.precision*mean.x)

	cat(paste("Posterior mean           : ",round(post.mean,7),"\n",sep=""))
	cat(paste("Posterior std. deviation : ",round(post.sd,7),"\n",sep=""))

	posterior<-dnorm(mu,post.mean,post.sd)
	
	plot(mu,posterior,ylim=c(0,1.1*max(posterior,mu.prior)),type="l",
		lty=1,col="blue",
		xlab=expression(mu),ylab=expression(Probabilty(mu)))
	lines(mu,mu.prior,lty=2,col="red")

	left<-min(mu)+diff(range(mu))*0.05
	legend(left,max(posterior,mu.prior),lty=1:2,col=c("blue","red"),legend=c("Posterior","Prior"))

	probs<-c(0.005,0.01,0.025,0.05,0.5,0.95,0.975,0.99,0.995)
	qtls<-qnorm(probs,post.mean,post.sd)
	names(qtls)<-probs

	cat("\nProb.\tQuantile \n")
	cat("\------\t---------\n")
	for(i in 1:length(probs))
		cat(paste(round(probs[i],3),"\t",round(qtls[i],7),"\n",sep=""))

	if(ret)
		return(list(prior=prior,likelihood=likelihood,posterior=posterior,mean=post.mean,sd=post.sd,quantiles=qtls))
	
}
