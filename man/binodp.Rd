\name{binodp}
\alias{binodp}
\title{Binomial sampling with a discrete prior}
\description{
	Evaluates and plots the posterior density for theta, the probability of a success in a Bernoulli trial, with binomial sampling and a discrete prior on theta
}
\usage{
binodp(x,n,uniform = TRUE, n.theta = 10, theta = NULL, theta.prior = NULL, ret = FALSE)
}
\arguments{
	\item{x}{the number of observed successes in the binomial experiment.}
	\item{n}{the number of trials in the binomial experiment.}
	\item{uniform}{if true then a discrete uniform prior for theta will be used.}
	\item{n.theta}{the number of possible theta values in the prior}
	\item{theta}{a vector of possibilities for the probability of success in a single trial.}
	\item{theta.prior}{the associated prior probability mass.}
	\item{ret}{if true then the likelihood and posterior are returned as a list.}
}
\seealso{
	\code{\link{binobp}}
	\code{\link{binogcp}}
}
\value{
	If ret is true, then a list will be returned with the following components:
	\item{f.cond}{the conditional distribution of x given theta and n}
	\item{f.joint}{the joint distribution of x and theta given n}
	\item{f.marg}{the marginal distribution of x}
	\item{likelihood}{the scaled likelihood function of x given theta and n}
	\item{posterior}{the posterior probability of theta given x and n}
	\item{theta}{the vector of possible theta values used in the prior}
	\item{theta.prior}{the associated probability mass for the values in theta}
}
\examples{
## simplest call with 6 successes observed in 8 trials and a uniform prior
binodp(6,8)

## same as previous example but with more possibilities for theta
binodp(6,8,n.theta=100)

## 6 successes, 8 trials and a non-uniform discrete prior
theta<-seq(0,1,by=0.01)
theta.prior<-runif(101)
theta.prior<-sort(theta.prior/sum(theta.prior))
binodp(6,8,uniform=FALSE,theta=theta,theta.prior=theta.prior)

## 5 successes, 6 trials, non-uniform prior
theta<-c(0.3,0.4,0.5)
theta.prior<-c(0.2,0.3,0.5)
results<-binodp(5,6,uniform=FALSE,theta=theta,theta.prior=theta.prior,ret=TRUE)

## plot the results from the previous example using a side-by-side barplot
results.matrix<-rbind(results$theta.prior,results$posterior)
colnames(results.matrix)<-theta
barplot(results.matrix,col=c("red","blue"),beside=TRUE
	,xlab=expression(theta),ylab=expression(Probability(theta)))
box()
legend(1,0.65,legend=c("Prior","Posterior"),fill=c("red","blue"))
}
\keyword{misc}

