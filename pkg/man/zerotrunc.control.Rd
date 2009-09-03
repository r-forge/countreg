\name{zerotrunc.control}
\alias{zerotrunc.control}

\title{Control Parameters for Zero-Truncated Count Data Regression}
\description{
  Various parameters that control fitting of zero-truncated count regression models
  using \code{\link{zerotrunc}}.
}

\usage{
zerotrunc.control(method = "BFGS", maxit = 10000, start = NULL, \dots)
}

\arguments{
  \item{method}{characters string specifying the \code{method} argument
    passed to \code{\link[stats]{optim}}.}
  \item{maxit}{integer specifying the \code{maxit} argument (maximal number
    of iterations) passed to \code{\link[stats]{optim}}.}
  \item{start}{an optional vector of starting values, see details.}
  \item{\dots}{arguments passed to \code{\link[stats]{optim}}.}
}

\details{
  All parameters in \code{\link{zerotrunc}} are estimated by maximum likelihood
  using \code{\link[stats]{optim}} with control options set in \code{\link{zerotrunc.control}}.
  Most arguments are passed on directly to \code{optim}, only \code{start} is used to
  control how \code{optim} is called.

  Starting values can be supplied via \code{start} or estimated by
  \code{\link[stats]{glm.fit}} (default). Standard errors are derived numerically using
  the Hessian matrix returned by \code{\link[stats]{optim}}. To supply starting 
  values, \code{start} should be a vector with (at least) starting values for the
  regression coefficients. In case a negative binomial distribution with unknown theta
  is used, a starting value for theta may be supplied by adding an additional vector
  element (e.g., \code{start = c(coef, theta)}); by default \code{theta = 1} is used as
  the starting value otherwise.
}

\value{
  A list with the arguments specified.
}

\author{Achim Zeileis <Achim.Zeileis@R-project.org>}

\seealso{\code{\link{zerotrunc}}}

\examples{
data("bioChemists", package = "pscl")
bc <- subset(bioChemists, subset = art > 0)

## default start values
zt_nb <- zerotrunc(art ~ fem + ment, data = bc, dist = "negbin")

## user-supplied start values and other options
zt_nb2 <- zerotrunc(art ~ fem + ment, data = bc, dist = "negbin",
  start = c(0.3, -0.2, 0, 2))
}

\keyword{regression}
