\name{nbreg.control}
\alias{nbreg.control}

\title{Control Parameters for Negative Binomial Count Data Regression}
\description{
  Various parameters that control fitting of negative binomial regression models
  using \code{\link{nbreg}}.
}

\usage{
nbreg.control(method = "BFGS", maxit = 10000, start = NULL, hessian = TRUE,
    dot = "separate", \dots)
}

\arguments{
  \item{method}{characters string specifying the \code{method} argument
    passed to \code{\link[stats]{optim}}.}
  \item{maxit}{integer specifying the \code{maxit} argument (maximal number
    of iterations) passed to \code{\link[stats]{optim}}.}
  \item{start}{an optional list with elements \code{"mu"} and \code{"theta"}
    containing the coefficients for the corresponding component.}
  \item{hessian}{logical. Should the numerically approximated Hessian be computed
    to derive an estimate of the variance-covariance matrix? If \code{FALSE} and
    parameter \code{hessA = FALSE} in \code{nbreg()}, the variance-covariance
    matrix contains only \code{NA}s.}
  \item{dot}{character. Controls how two-part \code{Formula}'s are processed.
    See \link{model.frame.Formula}.}
  \item{\dots}{arguments passed to \code{\link[stats]{optim}}.}
}

\details{
  All parameters in \code{\link{nbreg}} are estimated by maximum likelihood
  using \code{\link[stats]{optim}} with control options set in \code{\link{nbreg.control}}.
  Most arguments are passed on directly to \code{optim} and \code{start} controls
  the choice of starting values for calling \code{optim}.

  Starting values can be supplied or are estimated by a Poisson regression in
  \code{\link[stats]{glm.fit}} (the default, starting values of coefficients in \eqn{\theta}
  are set to zero to ensure compatibility with NB1). Standard errors are derived
  using the analytical Hessian matrix or by numerical approximation of the Hessian.
}

\value{
  A list with the arguments specified.
}

\seealso{\code{\link{nbreg}}}

\examples{
data("CrabSatellites", package = "countreg")

## default start values
fm1 <- nbreg(satellites ~ width + as.numeric(color), data = CrabSatellites)

## user-supplied start values
fm2 <- nbreg(satellites ~ width + as.numeric(color), data = CrabSatellites,
                start = list(mu = c(0, 0, 0), theta = c(0.5)))
}

\keyword{regression}
