binogcp<-function(x,n,density="uniform", params = c(0,1), n.theta = 1000, theta = NULL, theta.prior = NULL, ret = FALSE){

	# n - the number of trials in the binomial
	# x - the number of observed successes
	# density - may be one of "exp","normal","uniform" or "user"
	# params - if the density is not "user" then a vector of parameters
	# must be supplied. 
	#	exp:		rate
	#	normal: 	mean,sd
	#	uniform: 	min,max
	# n.theta - the number of points to divide the [0,1] interval into  

	# theta and theta.prior are only specified if density == "user"
	# theta - the probability of success 
	# theta.prior - the associated prior probability mass
	# ret - if true then the likelihood and posterior are returned as a 
	# list
	
	if(x>n)
		stop("The number of observed successes (x) must be smaller than the number of trials") 
	if(n.theta<100)
		stop("Number of prior values of theta must be greater than 100")
	
	if(is.null(theta)||is.null(theta.prior))
		theta<-seq(0+1/n.theta,1-1/n.theta,length=n.theta)
	else{
		if(length(theta)!=length(theta.prior))
			stop("theta and theta.prior must have same length")
		
		if(sum(theta<0|theta>1)>0) # check that probabilities lie on [0,1]
			stop("Values of theta must be between 0 and 1 inclusive")
	}

	if(density=="beta"){
		if(length(params)<2){
			warning("Beta prior requires two shape parameters. Default value Beta(1,1) = Uniform is being used")
			a<-1
			b<-1
		}else{
			if(params[1]<=0|params[2]<0)
				stop("Beta prior shape parameters must be greater than zero")
			a<-params[1]
			b<-params[2]
		}
		theta.prior<-dbeta(theta,a,b)
	}else	if(density=="exp"){
		if(params[1]<=0){
			stop("Parameter for exponential density must be greater than zero")
		}else{
			rate<-params[1]
			theta.prior<-dexp(theta,rate)
		}
	}else if(density=="normal"){
		if(length(params)<2)
			stop("Normal prior requires a mean and std. deviation")
		else{
			mx<-params[1]
			sx<-params[2]
			if(sx<=0)
				stop("Std. deviation for normal prior must be greater than zero")
			theta.prior<-dnorm(theta,mx,sx)
		}
	}else if(density=="uniform"){
		if(length(params)<2)
			stop("Uniform prior requires a minimum and a maximum")
		else{
			minx<-params[1]
			maxx<-params[2]
			
			if(maxx<=minx)
				stop("Maximum must be greater than minimum for a uniform prior")
			theta.prior<-dunif(theta,minx,maxx)
		}
	}else if (density!="user"){
		stop(paste("Unrecognized density :",density))
	}

	likelihood<-(theta^x)*((1-theta)^(n-x))

	# Numerically integrate the denominator
	# First calculate the height of the function to be integrated

	f.x.theta<-likelihood*theta.prior

	# Now get a linear approximation so that we don't have to worry about
	# the number of points specified by the user

	ap<-approx(theta,f.x.theta,n=513)
	integral<-sum(ap$y[2*(1:256)-1]+4*ap$y[2*(1:256)]+ap$y[2*(1:256)+1])
	integral<-(ap$x[2]-ap$x[1])*integral/3

	posterior<-likelihood*theta.prior/integral

	plot(theta,posterior,ylim=c(0,1.1*max(posterior,theta.prior)),lty=1,type="l",col="blue",
		xlab=expression(theta),ylab="Density")
	lines(theta,theta.prior,lty=2,col="red")
	
	left<-min(theta)+diff(range(theta))*0.05
	legend(left,max(posterior,theta.prior),lty=1:2,col=c("blue","red"),legend=c("Posterior","Prior"))
	if(ret)
		return(list(likelihood=likelihood,posterior=posterior,theta=theta,theta.prior=theta.prior))
	
}
