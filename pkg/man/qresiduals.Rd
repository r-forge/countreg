\name{qresiduals}
\alias{qresiduals}
\alias{qresiduals.default}

\title{(Randomized) Quantile Residuals}

\description{
  Generic function and methods for computing (randomized) quantile residuals.
}

\usage{
qresiduals(object, \dots)

\method{qresiduals}{default}(object, type = c("random", "quantile"), nsim = 1L, prob = 0.5, \dots)
}

\arguments{
  \item{object}{an object. For the \code{default} method this needs to
    be either a specification of probabilities (vector or 2-dimensional matrix of
    probabilities) or an object from which the these can be obtained with
    \code{\link{pit}}.}
  \item{type}{character specifying whether - in the case of discrete response distributions - 
    randomized quantile residuals or their corresponding quantiles should be computed.}
  \item{nsim}{numeric. The number of simulated randomized quantile residuals
    per observation (for \code{type = "numeric"}).}
  \item{prob}{numeric. The probabilities at which quantile residuals should be
    computed (for \code{type = "quantile"}), defaulting to the median.}
  \item{\dots}{further parameters passed to methods.}
}

\details{
(Randomized) quantile residuals have been suggested by Dunn and Smyth (1996).
For regression models with a continuous response distribution this simply computes
theoretical standard normal quantile corresponding to the probability integral transform of the
fitted distribution. For discrete distributions, a random theoretical normal
quantile is drawn from the range of probabilities corresponding to each observation.

The default \code{qresiduals} method can compute randomized quantile residuals
from a vector (which essentially just calls \code{\link[stats]{qnorm}}) or
a 2-column matrix of probabilities. The latter offers to either draw \code{"random"}
samples from the distribution or compute corresponding \code{"quantile"}s such as
the median etc.
}

\note{
Note that there is also a \code{\link[statmod]{qresiduals}} function in the \pkg{statmod}
package that is not generic and always returns a single random quantile residual.
}

\value{
  A vector or matrix of quantile residuals.
}

\references{
Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile Residuals.}
  \emph{Journal of Computational and Graphical Statistics}, \bold{5}, 1--10.
}

\seealso{\code{\link[stats]{qnorm}}, \code{\link{qqrplot}}}

\examples{
## count data regression models: crab satellites
data("CrabSatellites", package = "countreg")
cs_p <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)

qres <- cbind(
  sample = qresiduals(cs_p, nsim = 3),
  median = qresiduals(cs_p, type = "quantile"),
  mean100 = rowMeans(qresiduals(cs_p, nsim = 100)),
  range = qresiduals(cs_p, type = "quantile", prob = c(0, 1))
)
qres[1:5, ]
}

\keyword{regression}
