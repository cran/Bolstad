\name{binogcp}
\alias{binogcp}
\title{Binomial sampling with a general continuous prior}
\description{
	Evaluates and plots the posterior density for theta, the probability of a success in a Bernoulli trial, with binomial sampling and a general continuous prior on theta
}
\usage{
binogcp(x, n, density = "uniform", params = c(0,1), n.theta = 1000,
	theta = NULL, theta.prior = NULL, ret = FALSE)
}
\arguments{
	\item{x}{the number of observed successes in the binomial experiment.}
	\item{n}{the number of trials in the binomial experiment.}
	\item{density}{may be one of "beta", "exp", "normal", "student", "uniform" or "user"}
	\item{params}{if density is one of the parameteric forms then 
		then a vector of parameters  must be supplied. 
		beta:		a,b
		exp:		rate
		normal: 	mean,sd
		uniform: 	min,max}
	\item{n.theta}{the number of possible theta values in the prior}
	\item{theta}{a vector of possibilities for the probability of success in a single trial. This must be set if density="user".}
	\item{theta.prior}{the associated prior probability mass. This must be set if density="user".}
	\item{ret}{if true then the likelihood and posterior are returned as a list.}
}
\seealso{
	\code{\link{binobp}}
	\code{\link{binodp}}
}
\value{
	If ret is true, then a list will be returned with the following components:
	\item{likelihood}{the scaled likelihood function of x given theta and n}
	\item{posterior}{the posterior probability of theta given x and n}
	\item{theta}{the vector of possible theta values used in the prior}
	\item{theta.prior}{the associated probability mass for the values in theta}
}
\examples{
## simplest call with 6 successes observed in 8 trials and a continuous 
## uniform prior
binogcp(6,8)

## 6 successes, 8 trials and a Beta(2,2) prior
binogcp(6,8,density="beta",params=c(2,2))

## 5 successes, 10 trials and a N(0.5,0.25) prior
binogcp(5,10,density="normal",params=c(0.5,0.25))

## 4 successes, 12 trials with a user specified triangular continuous prior
theta<-seq(0,1,by=0.001)
theta.prior<-rep(0,length(theta))
theta.prior[theta<=0.5]<-4*theta[theta<=0.5]
theta.prior[theta>0.5]<-4-4*theta[theta>0.5]
results<-binogcp(4,12,"user",theta=theta,theta.prior=theta.prior,ret=TRUE)

## find the posterior CDF using the previous example and Simpson's rule
cdf<-sintegral(theta,results$posterior,n.pts=length(theta),ret=TRUE)
plot(cdf,type="l",xlab=expression(theta[0])
	,ylab=expression(Pr(theta<=theta[0])))

## use the cdf to find the 95\% credible region. Thanks to John Wilkinson for this simplified code.
lcb<-cdf$x[with(cdf,which.max(x[y<=0.025]))]
ucb<-cdf$x[with(cdf,which.max(x[y<=0.975]))]
cat(paste("Approximate 95\% credible interval : ["
	,round(lcb,4)," ",round(ucb,4),"]\n",sep=""))

## find the posterior mean, variance and std. deviation
## using Simpson's rule and the output from the previous example
dens<-theta*results$posterior # calculate theta*f(theta | x, n)
post.mean<-sintegral(theta,dens)

dens<-(theta-post.mean)^2*results$posterior
post.var<-sintegral(theta,dens)
post.sd<-sqrt(post.var)

# calculate an approximate 95\% credible region using the posterior mean and 
# std. deviation
lb<-post.mean-qnorm(0.975)*post.sd
ub<-post.mean+qnorm(0.975)*post.sd

cat(paste("Approximate 95\% credible interval : ["
	,round(lb,4)," ",round(ub,4),"]\n",sep=""))
}
\keyword{misc}


