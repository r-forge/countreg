\name{predict.nbreg}
\alias{predict.nbreg}
\alias{residuals.nbreg}
\alias{terms.nbreg}
\alias{model.frame.nbreg}
\alias{model.matrix.nbreg}
\alias{coef.nbreg}
\alias{vcov.nbreg}
\alias{summary.nbreg}
\alias{getSummary.nbreg}
\alias{print.summary.nbreg}
\alias{logLik.nbreg}
\alias{nobs.nbreg}
\alias{fitted.nbreg}
\alias{extractAIC.nbreg}

\title{Methods for nbreg Objects}
\description{
  Methods for extracting information from fitted negative binomial
  count regression model objects of class \code{"nbreg"}.
}

\usage{
\method{predict}{nbreg}(object, newdata,
  type = c("response", "prob", "theta", "parameters"), na.action = na.pass, \dots)
\method{residuals}{nbreg}(object, type = c("pearson", "deviance", "response"), \dots)

\method{coef}{nbreg}(object, model = c("full", "mu", "theta"), \dots)
\method{vcov}{nbreg}(object, model = c("full", "mu", "theta"), \dots)

\method{terms}{nbreg}(x, model = c("full", "mu", "theta"), \dots)
\method{model.matrix}{nbreg}(object, model = c("mu", "theta"), \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"nbreg"} as returned by
    \code{\link{nbreg}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{type}{character specifying the type of predictions or residuals,
    respectively. For details see below.}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}. The default is to predict \code{NA}.}
  \item{model}{character specifying for which component of the model the
    terms or model matrix should be extracted.}
  \item{\dots}{currently not used.}
}

\details{
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"nbreg"}, including methods to the generic functions
  \code{\link[base]{print}} and \code{\link[base]{summary}} which print the estimated
  coefficients along with some further information. The \code{summary} in particular
  supplies partial Wald tests based on the coefficients and the covariance matrix.
  As usual, the \code{summary} method returns an object of class \code{"summary.nbreg"}
  containing the relevant summary statistics which can subsequently be printed
  using the associated \code{print} method.
  
  The methods for \code{\link[stats]{coef}} and \code{\link[stats]{vcov}} by default
  return a single vector of coefficients and their associated covariance matrix,
  respectively, i.e., all coefficients are concatenated. By setting the \code{model}
  argument, the estimates for the corresponding model component can be extracted.
  
  Both the \code{\link[stats]{fitted}} and \code{\link[stats]{predict}} methods can
  compute fitted responses. The latter additionally provides the predicted density
  (i.e., probabilities for the observed counts) and the predicted dispersion
  parameter theta. The \code{\link[stats]{residuals}} method can compute raw residuals  
  (observed - fitted), Pearson residuals (raw residuals scaled by square root of
  variance function), and deviance residuals. The latter are only supported for
  negative binomial type 2 models \code{(dist = NB2)} (includes NBH).
  
  A \code{\link[stats]{logLik}} method is provided, hence \code{\link[stats]{AIC}}
  can be called to compute information criteria.
}

\seealso{\code{\link{nbreg}}}

\examples{
data("CrabSatellites", package = "countreg")
fm <- nbreg(satellites ~ width + color, data = CrabSatellites)

plot(residuals(fm, type = "pearson") ~ fitted(fm))

coef(fm)
summary(fm)
logLik(fm)
AIC(fm)
}

\keyword{regression}
