binodp<-function(x,n,uniform = TRUE, n.theta = 10, theta = NULL, theta.prior = NULL, ret = FALSE){

	# n - the number of trials in the binomial
	# x - the number of observed successes
	# theta - the probability of success 
	# theta.prior - the associated prior probability mass
	# ret - if true then the likelihood and posterior are returned as a 
	# list
	
	if(x>n)
		stop("The number of observed successes (x) must be smaller than the number of trials (n)") 
	if(n.theta<3)
		stop("Number of prior values of theta must be greater than 2")

	if(is.null(theta)&is.null(theta.prior)&uniform){
		theta<-seq(0,1, length = n.theta)
		theta.prior<-rep(1/n.theta,n.theta)
	}

	if(sum(theta<0) > 0 | sum(theta > 1) > 0) # check that probabilities lie on [0,1]
		stop("Values of theta must be between 0 and 1 inclusive")

	if(sum(theta.prior<0)>0|sum(theta.prior>1)>1)
		stop("Prior probabilities must be between 0 and 1 inclusive")

	if(round(sum(theta.prior),7)!=1){
		warning("The prior probabilities did not sum to 1, therefore the prior has been normalized")
		theta.prior<-theta.prior/sum(theta.prior)
	}

	if(!uniform&(is.null(theta)|is.null(theta.prior)))
		stop("If you wish to use a non-uniform discrete prior then you must supply a theta vector and an associated probability vector theta.prior")


	n.theta<-length(theta)
	likelihood<-dbinom(x,n,theta)*theta.prior


	posterior<-likelihood/sum(likelihood)

	plot(theta,posterior,ylim=c(0,1.1*max(posterior,theta.prior)),pch="o",
		xlab=expression(theta),ylab=expression(Probabilty(theta)))
	points(theta,theta.prior,pch="+")

	legend(max(c(0.05,min(theta))),max(posterior,theta.prior),pch=c("o","+"),legend=c("Posterior","Prior"))

	# calculate the Conditional distribution

	f.cond<-matrix(0,nrow=n.theta,ncol=n+1)
	rownames(f.cond)<-as.character(round(theta,3))
	colnames(f.cond)<-as.character(0:n)

	for(i in 1:n.theta)
		f.cond[i,]<-dbinom(0:n,n,theta[i])

	cat("Conditional distribution of x given theta and  n:\n\n")
	print(round(f.cond,4))	

	# caculate the joint distribution of theta and x given n

	f.joint<-diag(theta.prior)%*%f.cond
	cat("\nJoint distribution:\n\n")
	print(round(f.joint,4))	

	# calculate the marginal distribtion

	f.marg<-matrix(1,nrow=1,ncol=n.theta)%*%f.joint
	cat("\nMarginal distribution of x:\n\n")
	print(round(f.marg,4))
	cat("\n\n")	

	# finally display the prior, likelihood, and posterior

	results<-cbind(theta.prior,likelihood,posterior)
	rownames(results)<-as.character(round(theta,3))
	colnames(results)<-c("Prior","Likelihood","Posterior")

	print(results)






	if(ret)
		return(list(f.cond=f.cond,f.joint=f.joint,f.marg=f.marg,likelihood=likelihood,posterior=posterior,theta=theta,theta.prior=theta.prior))
	
}
