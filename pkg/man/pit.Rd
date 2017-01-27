\name{pit}
\alias{pit}

\alias{pit.default}
\alias{pit.glm}
\alias{pit.hurdle}
\alias{pit.zeroinfl}
\alias{pit.zerotrunc}

\title{Probability Integral Transform (PIT)}

\description{
  Generic function and methods for computing probability integral transform (PIT)
  from fitted model objects.
}

\usage{
pit(object, \dots)
}

\arguments{
  \item{object}{a fitted model object.}
  \item{\dots}{further parameters passed to methods.}
}

\details{
The probability integral transform (PIT) is defined as \eqn{\hat F(y)}, i.e., a fitted
or estimated cumulative distribution function \eqn{\hat F(\cdot)} evaluated at
observations \eqn{y}. It is a basic quantity employed for model diagnostics, e.g.,
with quantile residuals (\code{\link{qresiduals}}), corresponding Q-Q plots
(\code{\link{qqrplot}}), or the PIT histogram (\code{\link{pithist}}).

In case of discrete distributions \code{pit} should compute a matrix with two-columns:
the left-handed limit of the PIT and the PIT at the observed response, respectively.
}

\value{
  A vector or matrix of probabilities.
}

\seealso{\code{\link{qresiduals}}, \code{\link{qqrplot}}}

\examples{
## count data regression models: crab satellites
data("CrabSatellites", package = "countreg")
cs_p <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
pit(cs_p)[1:5,]
}

\keyword{regression}
