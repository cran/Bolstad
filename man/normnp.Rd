\name{normnp}
\alias{normnp}
\title{Bayesian inference on a normal mean with a normal prior}
\description{
	Evaluates and plots the posterior density for mu, the mean of a normal distribution, with a normal prior on mu
}
\usage{
normnp(x, sigma.x, m.x = 0, s.x = 1, n.mu = 100, ret = FALSE)
}
\arguments{
	\item{x}{a vector of observations from a normal distribution with unknown mean and known std. deviation.}
	\item{sigma.x}{the population std. deviation of the normal distribution}	
	\item{m.x}{the mean of the normal prior}
	\item{s.x}{the standard deviation of the normal prior}
	\item{n.mu}{the number of possible mu values in the prior}
	\item{ret}{if true then the likelihood and posterior are returned as a list.}
}
\seealso{
	\code{\link{normdp}}
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

## find the posterior density with a N(0,1) prior on mu
normnp(x,1)

## find the posterior density with N(0.5,3) prior on mu
normnp(x,1,0.5,3)

## Find the posterior density for mu, given a random sample of 4 
## observations from N(mu,sigma^2=1), y = [2.99, 5.56, 2.83, 3.47], 
## and a N(3,sd=2)$ prior for mu
y<-c(2.99,5.56,2.83,3.47)
normnp(y,1,3,2)
}
\keyword{misc}

