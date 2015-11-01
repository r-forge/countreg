\name{qresiduals}
\alias{qresiduals}

\alias{qresiduals.default}
\alias{qresiduals.glm}
\alias{qresiduals.hurdle}
\alias{qresiduals.zeroinfl}
\alias{qresiduals.zerotrunc}

\title{(Randomized) Quantile Residuals}

\description{
  Generic function and methods for computing (randomized) quantile residuals.
}

\usage{
qresiduals(object, \dots)

\method{qresiduals}{default}(object, aggregate = c("median", "mean", "range", "none"),
  nsim = 10L, \dots)
}

\arguments{
  \item{object}{an object. For the \code{default} method this needs to
    be a vector or 2-dimensional matrix of probabilities. Currently, there
    are also methods for \code{lm}, \code{glm}, \code{hurdle}, \code{zeroinfl}, and
    \code{zerotrunc} objects.}
  \item{aggregate}{character specifying which type of aggregation should be used
    for the distribution of the randomized quantile residuals in case of discrete
    response distributions.}
  \item{nsim}{numeric. The number of simulated randomized quantile residuals
    per observation (for \code{aggregate = "mean"} or \code{"none"}).}
  \item{\dots}{further parameters passed to methods.}
}

\details{
(Randomized) quantile residuals have been suggested by Dunn and Smyth (1996).
For regression models with a continuous response distribution this simply computes
theoretical standard normal quantile corresponding to the probabilities of the
fitted distribution. For discrete distributions, a random theoretical normal
quantile is drawn from the range of probabilities corresponding to each observation.

The default \code{qresiduals} method can compute randomized quantile residuals
from a vector (which essentially just call \code{\link[stats]{qnorm}}) or
a 2-column matrix of probabilities. The latter offers various options for aggregation:
The default is to return only the \code{"median"} of the quantile residual distribution.
Alternatively, the \code{"range"} can be returned or the \code{"mean"} from \code{nsim}
simulations. Finally, one or more random quantile residuals per observation can
be simulated if \code{aggregate = "none"}.

Currently, there are methods for \code{\link[stats]{lm}},  \code{\link[stats]{glm}}, 
\code{\link{hurdle}}, \code{\link{zeroinfl}}, and \code{\link{zerotrunc}} models. All
of these compute the empirical probabilities first and then call the default method.
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
  median = qresiduals(cs_p),
  mean10 = qresiduals(cs_p, aggregate = "mean"),
  mean100 = qresiduals(cs_p, aggregate = "mean", nsim = 100),
  range = qresiduals(cs_p, aggregate = "range"),
  sample = qresiduals(cs_p, aggregate = "none", nsim = 3)
)
qres[1:5, ]

}

\keyword{regression}
