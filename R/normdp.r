normdp<-function(x,sigma.x,uniform = TRUE, n.mu = 50, mu = NULL, mu.prior = NULL, ret = FALSE){

	# x - the vector of observations
	# sigma.x - the population standard deviation
	# uniform - default to a discrete uniform prior
	# mu - vector of possible values of the population mean 
	# mu.prior - the associated prior probability mass
	# ret - if true then the likelihood and posterior are returned as a 
	# list
	
	if(n.mu<3)
		stop("Number of prior values of theta must be greater than 2")

	if(uniform){
		mu<-seq(min(x)-sigma.x,max(x)+sigma.x,length = n.mu)
		mu.prior<-rep(1/n.mu,n.mu)
	}

	if(any(mu.prior<0) | any(mu.prior>1))
		stop("Prior probabilities must be between 0 and 1 inclusive")

	if(round(sum(mu.prior),7)!=1){
		warning("The prior probabilities did not sum to 1, therefore the prior has been normalized")
		mu.prior<-mu.prior/sum(mu.prior)
	}

	if(!uniform&(is.null(mu)|is.null(mu.prior)))
		stop("If you wish to use a non-uniform discrete prior then you must supply a mean vector, mu, and an associated probability vector, mu.prior")


	n.mu<-length(mu)
	mx<-mean(x)
	nx<-length(x)
	snx<-sigma.x^2/nx
	likelihood<-exp(-0.5*(mx-mu)^2/snx)

	posterior<-likelihood*mu.prior/sum(likelihood*mu.prior)

	plot(mu,posterior,ylim=c(0,1.1*max(posterior,mu.prior)),pch="o",
		xlab=expression(mu),ylab=expression(Probabilty(mu)))
	points(mu,mu.prior,pch="+")
	
	left<-min(mu)+diff(range(mu))*0.05
	legend(left,max(posterior,mu.prior),pch=c("o","+"),legend=c("Posterior","Prior"))

	if(ret)
		return(list(likelihood=likelihood,posterior=posterior,mu=mu,mu.prior=mu.prior))
	
}
