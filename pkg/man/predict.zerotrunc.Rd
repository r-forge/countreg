\name{predict.zerotrunc}
\alias{predict.zerotrunc}
\alias{residuals.zerotrunc}
\alias{estfun.zerotrunc}
\alias{terms.zerotrunc}
\alias{model.frame.zerotrunc}
\alias{model.matrix.zerotrunc}
\alias{coef.zerotrunc}
\alias{vcov.zerotrunc}
\alias{summary.zerotrunc}
\alias{print.summary.zerotrunc}
\alias{logLik.zerotrunc}
\alias{nobs.zerotrunc}
\alias{fitted.zerotrunc}
\alias{predprob.zerotrunc}
\alias{extractAIC.zerotrunc}
\alias{getSummary.zerotrunc}

\title{Methods for zerotrunc Objects}
\description{
  Methods for extracting information from fitted zero-truncated
  count regression model objects of class \code{"zerotrunc"}.
}

\usage{
\method{predict}{zerotrunc}(object, newdata,
  type = c("response", "prob", "count", "zero"), na.action = na.pass, \dots)
\method{residuals}{zerotrunc}(object, type = c("deviance", "pearson", "response"), \dots)
}

\arguments{
  \item{object}{an object of class \code{"zerotrunc"} as returned by
    \code{\link{zerotrunc}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{type}{character specifying the type of predictions or residuals,
    respectively. For details see below.}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}. The default is to predict \code{NA}.}
  \item{\dots}{currently not used.}
}

\details{
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"zerotrunc"}, including methods to the generic functions
  \code{\link[base]{print}} and \code{\link[base]{summary}} which print the estimated
  coefficients along with some further information. The \code{summary} in particular
  supplies partial Wald tests based on the coefficients and the covariance matrix
  (estimated from the Hessian in the numerical optimization of the log-likelihood).
  As usual, the \code{summary} method returns an object of class \code{"summary.zerotrunc"}
  containing the relevant summary statistics which can subsequently be printed
  using the associated \code{print} method.
  
  Both the \code{\link[stats]{fitted}} and \code{\link[stats]{predict}} methods can
  compute fitted responses. The latter additionally provides the predicted density
  (i.e., probabilities for the observed counts), the predicted mean from the count
  component (without zero truncation) and the predicted probability for observing
  a non-zero count (in the un-truncated model). The \code{\link[stats]{residuals}}
  method can compute raw residuals (observed - fitted), Pearson residuals
  (raw residuals scaled by square root of variance function), and deviance
  residuals (contributions to the centered log-likelihood).
  
  A \code{\link[stats]{logLik}} method is provided, hence \code{\link[stats]{AIC}}
  can be called to compute information criteria.
}

\seealso{\code{\link{zerotrunc}}}

\examples{
data("CrabSatellites", package = "countreg")
fm <- zerotrunc(satellites ~ width + color, data = CrabSatellites, subset = satellites > 0)

plot(residuals(fm, type = "deviance") ~ fitted(fm))
plot(residuals(fm, type = "pearson") ~ fitted(fm))

coef(fm)
summary(fm)
logLik(fm)
AIC(fm)
}

\keyword{regression}
