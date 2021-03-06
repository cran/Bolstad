% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bolstad-package.R
\docType{data}
\name{slug}
\alias{slug}
\title{Slug data}
\format{
A data frame with 100 observations on the following 4 variables.
\itemize{ \item{length. length (mm) of the slug}
\item{weight. weight (g) of the slug} \item{log.len. natural
logarithm of the \code{length}} \item{log.wt. natural logarithm of
the \code{weight}} }
}
\description{
Lengths and weights of 100 slugs from the species Limax maximus collected
around Hamilton, New Zealand.
}
\examples{

data(slug)
plot(weight~length, data = slug)
plot(log.wt~log.len, data = slug)


}
\references{
Barker, G. and McGhie, R. (1984). The Biology of Introduced
Slugs (Pulmonata) in New Zealand: Introduction and Notes on Limax Maximus,
NZ Entomologist 8, 106--111.
}
\keyword{datasets}
