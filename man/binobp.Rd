\name{binobp}
\alias{binobp}
\title{Binomial sampling with a beta prior}
\description{
	Evaluates and plots the posterior density for theta, the probability of a success in a Bernoulli trial, with binomial sampling and a continous Beta(a,b) prior.
}
\usage{
binobp(x, n, a = 1, b = 1, ret = FALSE)
}
\arguments{
	\item{x}{the number of observed successes in the binomial experiment.}
	\item{n}{the number of trials in the binomial experiment.}
	\item{a}{parameter for the beta prior - must be greater than zero}
	\item{b}{parameter for the beta prior - must be greater than zero}
	\item{ret}{if true then the prior and posterior are returned as a list.}
}
\seealso{
	\code{\link{binodp}}
	\code{\link{binogcp}}
}
\value{
	If ret is true, then a list will be returned with the following components:
	\item{posterior}{the posterior density of theta given x and n = the Beta(a+x,b+n-x) density}
	\item{likelihood}{the likelihood of x given theta and n = the Binomial(n,theta) density}
	\item{prior}{the prior density of theta = the Beta(a,b) density}
	\item{theta}{the values of theta for which the posterior density was evaluated}
	\item{mean}{the posterior mean}
	\item{var}{the posterior variance}
	\item{sd}{the posterior std. deviation}
	\item{quantiles}{a set of quantiles from the posterior}
}
\examples{
## simplest call with 6 successes observed in 8 trials and a Beta(1,1) uniform 
## prior
binobp(6,8)

## 6 successes observed in 8 trials and a non-uniform Beta(0.5,6) prior
binobp(6,8,0.5,6)

## 4 successes observed in 12 trials with a non uniform Beta(3,3) prior
## plot the stored prior, likelihood and posterior
results<-binobp(4,12,3,3,ret=TRUE)

par(mfrow=c(3,1))
y.lims<-c(0,1.1*max(results$posterior,results$prior))

plot(results$theta,results$prior,ylim=y.lims,type="l"
	,xlab=expression(theta),ylab="Density",main="Prior")
polygon(results$theta,results$prior,col="red")

plot(results$theta,results$likelihood,ylim=c(0,0.25),type="l"
	,xlab=expression(theta),ylab="Density",main="Likelihood")
polygon(results$theta,results$likelihood,col="green")

plot(results$theta,results$posterior,ylim=y.lims,type="l"
	,xlab=expression(theta),ylab="Density",main="Posterior")
polygon(results$theta,results$posterior,col="blue")



}
\keyword{misc}

