\name{zeroinfl.control}
\alias{zeroinfl.control}

\title{Control Parameters for Zero-inflated Count Data Regression}
\description{
  Various parameters that control fitting of zero-inflated regression models
  using \code{\link{zeroinfl}}.
}

\usage{
zeroinfl.control(method = "BFGS", maxit = 10000, trace = FALSE,
  EM = FALSE, start = NULL, hessian = TRUE, \dots)
}

\arguments{
  \item{method}{characters string specifying the \code{method} argument
    passed to \code{\link[stats]{optim}}.}
  \item{maxit}{integer specifying the \code{maxit} argument (maximal number
    of iterations) passed to \code{\link[stats]{optim}}.}
  \item{trace}{logical or integer controlling whether tracing information on 
    the progress of the optimization should be produced (passed to \code{\link[stats]{optim}}).}
  \item{EM}{logical. Should starting values be estimated by the EM (expectation
    maximization) algorithm? See details.}
  \item{start}{an optional list with elements \code{"count"} and \code{"zero"}
    (and potentially \code{"theta"}) containing the coefficients for the corresponding component.}
  \item{hessian}{logical. Should the Hessian be computed to derive an estimate of the
    variance-covariance matrix? If \code{FALSE}, the variance-covariance matrix contains only
    \code{NA}s.}
  \item{\dots}{arguments passed to \code{\link[stats]{optim}}.}
}

\details{
  All parameters in \code{\link{zeroinfl}} are estimated by maximum likelihood
  using \code{\link[stats]{optim}} with control options set in \code{\link{zeroinfl.control}}.
  Most arguments are passed on directly to \code{optim}, only \code{trace} is also
  used within \code{zeroinfl} and \code{EM}/\code{start} control the choice
  of starting values for calling \code{optim}.

  Starting values can be supplied, estimated by the EM (expectation maximization)
  algorithm, or by \code{\link[stats]{glm.fit}} (the default). Standard errors are
  derived numerically using
  the Hessian matrix returned by \code{\link[stats]{optim}}. To supply starting 
  values, \code{start} should be a list with elements \code{"count"} and \code{"zero"}
  and potentially \code{"theta"} (for negative binomial components only) containing
  the starting values for the coefficients of the corresponding component of the
  model.
}

\value{
  A list with the arguments specified.
}

\seealso{\code{\link{zeroinfl}}}

\examples{
data("CrabSatellites", package = "countreg")

## default start values
fm1 <- zeroinfl(satellites ~ width + as.numeric(color), data = CrabSatellites)

## use EM algorithm for start values
fm2 <- zeroinfl(satellites ~ width + as.numeric(color), data = CrabSatellites, EM = TRUE)

## user-supplied start values
fm3 <- zeroinfl(satellites ~ width + as.numeric(color), data = CrabSatellites,
  start = list(count = c(0.5, 0, 0), zero = c(10, -0.5, 0.5)))
}

\keyword{regression}
