\name{normdp}
\alias{normdp}
\title{Bayesian inference on a normal mean with a discrete prior}
\description{
	Evaluates and plots the posterior density for mu, the mean of a normal distribution, with a discrete prior on mu
}
\usage{
normdp(x, sigma.x, uniform = TRUE, n.mu = 50, mu = NULL, mu.prior = NULL,
	ret = FALSE)
}
\arguments{
	\item{x}{a vector of observations from a normal distribution with unknown mean and known std. deviation.}
	\item{sigma.x}{the population std. deviation of the normal distribution}
	\item{uniform}{if true then a discrete uniform prior for mu will be used.}
	\item{n.mu}{the number of possible mu values in the prior}
	\item{mu}{a vector of possibilities for the probability of success in a single trial.}
	\item{mu.prior}{the associated prior probability mass.}
	\item{ret}{if true then the likelihood and posterior are returned as a list.}
}
\seealso{
	\code{\link{normnp}}
	\code{\link{normgcp}}
}
\value{
	If ret is true, then a list will be returned with the following components:
	\item{likelihood}{the scaled likelihood function of x given mu and sigma.x}
	\item{posterior}{the posterior probability of mu given x and sigma.x}
	\item{mu}{the vector of possible mu values used in the prior}
	\item{mu.prior}{the associated probability mass for the values in mu}
}
\examples{
## generate a sample of 20 observations from a N(-0.5,1) population
x<-rnorm(20,-0.5,1)

## find the posterior density with a uniform prior on mu
normdp(x,1)

## find the posterior density with a non-uniform prior on mu
mu<-seq(-3,3,by=0.1)
mu.prior<-runif(length(mu))
mu.prior<-sort(mu.prior/sum(mu.prior))
normdp(x,1,uniform=FALSE,mu=mu,mu.prior=mu.prior)

## Let mu have the discrete distribution with 5 possible
## values, 2, 2.5, 3, 3.5 and 4, and associated prior probability of
## 0.1, 0.2, 0.4, 0.2, 0.1 respectively. Find the posterior 
## distribution after a drawing random sample of n = 5 observations 
## from a N(mu,1) distribution y = [1.52, 0.02, 3.35, 3.49, 1.82]
mu<-seq(2,4,by=0.5)
mu.prior<-c(0.1,0.2,0.4,0.2,0.1)
y<-c(1.52,0.02,3.35,3.49,1.82)
normdp(y,1,uniform=FALSE,n.mu=5,mu,mu.prior)
}
\keyword{misc}

