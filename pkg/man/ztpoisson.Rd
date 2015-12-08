\name{ztpoisson}
\alias{ztpoisson}

\title{Family Object for the Zero-Truncated Poisson Distribution}

\description{
  Family object for specification of zero-truncated Poisson models
  as a \code{\link[stats]{glm}}.
}

\usage{
ztpoisson()
}

\details{
The \code{ztpoisson} family allows to estimate zero-truncated
Poisson regression models as generalized linear models. As in
the \code{\link{zerotrunc}} function, the link function is a
log-link between the mean \eqn{\lambda} of the untruncated Poisson
distribution and the linear predictor. This corresponds to a
non-canonical link between for the mean of the zero-truncated
Poisson distribution which does not have a closed-form representation.
}

\value{
  An object of class \code{"family"}.
}


\seealso{\code{\link{dztpois}}, \code{\link[stats]{poisson}},
  \code{\link{zerotrunc}}
}

\examples{
## data
data("CrabSatellites", package = "countreg")
cs <- subset(CrabSatellites, subset = satellites > 0)
cs$color <- as.numeric(cs$color)

## model
ztp1 <- glm(satellites ~ width + color, data = cs, family = ztpoisson)
ztp2 <- zerotrunc(satellites ~ width + color, data = cs)
summary(ztp1)
summary(ztp2)
}

\keyword{models}
