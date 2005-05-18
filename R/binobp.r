binobp<-function(x,n,a = 1, b = 1 , ret = FALSE){

	# n - the number of trials in the binomial
	# x - the number of observed successes
	# a,b  - the parameters of the Beta prior density (must be > 0)
	# ret - if true then the prior, likelihood, posterior, mean, variance and
  # std. deviation are returned as a list

	if(x>n)
		stop("The number of observed successes (x) must be smaller than the number of trials (n)")
	if(a<=0||b<=0)
		stop("The parameters of the prior must be greater than zero")

	theta<-seq(0.01,0.999,by=0.001)
	prior<-dbeta(theta,a,b)
	likelihood<-dbinom(x,n,prob=theta)
	posterior<-dbeta(theta,a+x,b+n-x)

	plot(theta,posterior,ylim=c(0,1.1*max(posterior,prior)),type="l"
		,lty=1
		,xlab=expression(theta)
		,ylab="Density")
	lines(theta,prior,lty=2)
	
	left<-min(theta)+diff(range(theta))*0.05
	legend(left,max(posterior,prior),lty=1:2,legend=c("Posterior","Prior"))
	
	m1<-(a+x)/(a+b+n)
	v1<-m1*(1-m1)/(a+b+n+1)
	s1<-sqrt(v1)

	cat(paste("Posterior Mean           : ",round(m1,7),"\n"))
	cat(paste("Posterior Variance       : ",round(v1,7),"\n"))
	cat(paste("Posterior Std. Deviation : ",round(s1,7),"\n"))

	probs<-c(0.005,0.01,0.025,0.05,0.5,0.95,0.975,0.99,0.995)
	qtls<-qbeta(probs,a+x,b+n-x)
	names(qtls)<-probs

	cat("\nProb.\tQuantile \n")
	cat("\------\t---------\n")
	for(i in 1:length(probs))
		cat(paste(round(probs[i],3),"\t",round(qtls[i],7),"\n",sep=""))


	if(ret)
		return(list(posterior=posterior,likelihood=likelihood,prior=prior,theta=theta,mean=m1,var=v1,sd=s1,quantiles=qtls))
	
}