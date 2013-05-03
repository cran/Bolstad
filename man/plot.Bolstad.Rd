\name{plot.Bolstad}
\alias{plot.Bolstad}
\title{Plot method for objects of type Bolstad}
\description{
A unified plotting method for plotting the prior, likelihood and posterior from any of the analyses in the book
}
\usage{
\method{plot}{Bolstad}(x, overlay = TRUE, which = c(1,3), densCols = c("red","green","blue")[which],  \dots)
}
\arguments{
  \item{x}{A S3 object of class Bolstad}
  \item{overlay}{if \code{FALSE} then up to three plots will be drawn side-by-side}
  \item{which}{Control which of the prior = 1, likelihood = 2, and posterior = 3, are plots. This is set to prior and posterior by     default to retain compatibility with the book}
  \item{densCols}{The colors of the lines for each of the prior, likelihood and posterior}
  \item{\dots}{Any remaining arguments are fed to the \code{plot} command}
  
}
\details{
The function provides a unified way of plotting the prior, likelihood and posterior from any of the functions in the library that return these quantities. It will produce an overlay of the lines by default, or separate panels if \code{overlay = FALSE}.
}
\author{James Curran}

\examples{
x = rnorm(20,-0.5,1)
## find the posterior density with a N(0,1) prior on mu
b = normnp(x,sigma=1)
plot(b)
plot(b, which = 1:3)
plot(b, overlay = FALSE, which = 1:3)
}
\keyword{plot}
