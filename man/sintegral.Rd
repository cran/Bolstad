\name{sintegral}
\alias{sintegral}
\title{Numerical integration using Simpson's Rule}
\description{
  Takes a vector of \eqn{x} values and a corresponding set of postive
  \eqn{f(x)=y} values and evaluates the area under the curve:
  \deqn{
    \int{f(x)dx}
    }.
}
\usage{
	sintegral(x, fx, n.pts = 256, ret = FALSE)
}
\arguments{
	\item{x}{a sequence of \eqn{x} values.}
	\item{fx}{the value of the function to be integrated at \eqn{x}.}
	\item{n.pts}{the number of points to be used in the integration.}
	\item{ret}{if true returns the partial sums of the integration.}
}
\examples{
## integrate the normal density from -3 to 3
x<-seq(-3,3,length=100)
fx<-dnorm(x)
estimate<-sintegral(x,fx)
true.val<-diff(pnorm(c(-3,3)))
cat(paste("Absolute error :",round(abs(estimate-true.val),7),"\n"))
cat(paste("Relative percentage error :", 100*round((abs(estimate-true.val)/true.val),6),"\%\n"))
}
\keyword{misc}


