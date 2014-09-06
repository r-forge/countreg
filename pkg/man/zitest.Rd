\name{zitest}
\alias{zitest}

\title{Zero Inflation Tests}

\description{
 Tests the null hypothesis of a Poisson GLM against
 the alternative of a zero-inflated version.
}

\usage{
zitest(object, type = c("scoreZIP"))
}

\arguments{
 \item{object}{a fitted Poisson GLM of class \code{"glm"} as fitted
    by \code{\link{glm}} with family \code{\link{poisson}}.}
 \item{type}{type of test, currently only \code{scoreZIP}. See details.}     
}

\details{
Currently alternative contains only intercept term in binary part, as in van den Broek (1995).

Note that under the null hypothesis the parameter is on the boundary of the parameter space, hence the p-value is non-standard.

}

\value{An object of class \code{"htest"}.}

\references{
van den Broek J (1995). \dQuote{A Score Test for Zero Inflation in a Poisson Distribution}. 
\emph{Biometrics}, \bold{51}, 738--743.

}

\seealso{\code{\link{glm}}, \code{\link{poisson}}, \code{\link[MASS]{glm.nb}}}

\examples{
data("CrabSatellites", package = "countreg")
CrabSatellites <- transform(CrabSatellites,
  color = as.numeric(color),
  spine = as.numeric(spine),
  cwidth = cut(width, c(-Inf, seq(23.25, 29.25), Inf))
)

cs_p <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
zitest(cs_p)

}

\keyword{htest}
