\name{nbreg}
\alias{nbreg}
\alias{print.nbreg}

\title{Negative Binomial Count Data Regression}

\description{
  Fit negative binomial regression models for count data via maximum likelihood
}

\usage{
nbreg(formula, data, subset, na.action, weights, offset, theta = NULL,
  dist = "NB2", link = "log", link.theta = "log", control = nbreg.control(\dots),
  model = TRUE, y = TRUE, x = FALSE, z = FALSE, hessA = TRUE, \dots)
}

\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor. See below for more information on offsets.}
  \item{theta}{numeric. Optional. If specified, then the dispersion parameter is 
    not estimated.}
  \item{dist}{character specification of the NB type. Either \code{"NB2"} or \code{"NB1"}}
  \item{link}{character specification of the link function for the mean.
    Currently, only \code{"log"} is supported.}
  \item{link.theta}{character specification of the link function for the
    dispersion parameter. Currently, only \code{"log"} is supported.}
  \item{control}{a list of control arguments specified via
    \code{\link{nbreg.control}}.}
  \item{model, y, x, z}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.}
  \item{hessA}{logical. If \code{TRUE}, then the analytical Hessian is used to
    compute the covariance matrix of the estimator.}
  \item{\dots}{currently not used.}
}

\details{
  The Negative Binomial Distribution is often used to model count data with overdispersion.
  Cameron and Trivedi (2013) offer two parametrization, negative binomial type 2 (NB2) and
  type 1 (NB1). NB2 is parametrized as follows
  
  \deqn{f(y | \mu, \theta) = \dfrac{\Gamma(\theta + y)}{\Gamma(\theta) y!}  \left(\dfrac{\theta}{\theta + \mu}\right)^\theta  \left(\dfrac{\mu}{\theta + \mu}\right)^y.}
  
  For NB1 replace \eqn{\theta} with \eqn{\mu\theta} on the RHS.
  
  This function further allows us to model the dispersion parameter with
  covariates via a two-part \code{formula}. If a \code{formula} of type
  \code{y ~ x1 + x2} is supplied, then the regressors are employed in the mean
  and the dispersion parameter is estimated as a constant, i.e. a standard
  NB2 or NB1 is estimated. This is equivalent to \code{y ~ x1 + x2 | 1}. If a
  \code{formula} of type \code{y ~ x1 + x2 + x3 | z1 + z2} is given, then the
  mean \eqn{mu} is modeled using \code{x1 + x2} and the dispersion parameter \eqn{\theta}
  with \code{z1 + z2}. If \code{dist = "NB2"}, then the function estimates the NBH model.
  
  Offsets can be specified in both the mean and dispersion parameter \eqn{\theta}:
  \code{y ~ x1 + x2 + offset(x3) | z1 + offset(z2)}, where
  \code{x3} is used as an offset (i.e., with coefficient fixed to 1) in the
  mean \eqn{mu} and \code{z2} analogously in \eqn{\theta}. By the rule
  stated above \code{y ~ x1 + offset(x2)} is equivalent to
  \code{y ~ x1 + offset(x2) | 1}. Instead of using the \code{offset()} wrapper
  within the \code{formula}, the \code{offset} argument can also be employed 
  which sets an offset only for \eqn{mu}. Thus, \code{formula = y ~ x1} and 
  \code{offset = x2} is equivalent to \code{formula = y ~ x1 + offset(x2) | 1}.  
  
  All parameters are estimated by maximum likelihood using \code{\link[stats]{optim}},
  with control options set in \code{\link{nbreg.control}}.
  Starting values can be supplied or are estimated by a Poisson regression in
  \code{\link[stats]{glm.fit}} (the default, starting values of coefficients in
  \eqn{\theta} are set to zero to ensure compatibility with NB1). Standard errors
  are derived analytically or numerically using the Hessian matrix returned by      
  \code{\link[stats]{optim}}.
  See \code{\link{nbreg.control}} for details.
  
  The returned fitted model object is of class \code{"nbreg"} and is similar
  to fitted \code{"glm"} objects.
  
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"nbreg"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{fitted}}, \code{\link[stats]{terms}},
  \code{\link[stats]{model.matrix}}. See \code{\link{predict.nbreg}} for more details
  on all methods.
  
}
\value{
  An object of class \code{"nbreg"}, i.e., a list with components including
  \item{coefficients}{a vector containing the coefficients from the mean,}
  \item{coefficients.theta}{a vector containing the coefficients from the
    dispersion parameter theta,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a vector of fitted means,}
  \item{optim}{a list with the output from the \code{optim} call for
    maximizing the log-likelihood,}
  \item{control}{the control arguments passed to the \code{optim} call,}
  \item{start}{the starting values for the parameters passed to the \code{optim} call,}
  \item{weights}{the case weights used,}
  \item{offset}{a list with elements \code{"mu"} and \code{"theta"}
    containing the offset vectors (if any) from the respective parameters,}
  \item{n}{number of observations (with weights > 0),}
  \item{df.null}{residual degrees of freedom for the null model,}
  \item{df.residual}{residual degrees of freedom for fitted model,}
  \item{terms}{a list with elements \code{"mu"}, \code{"theta"} and
    \code{"full"} containing the terms objects for the respective parameters,}
  \item{SE.logtheta}{standard error for \eqn{\log(\theta)}{log(theta)},}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{vcov}{covariance matrix of all coefficients in the model (derived from the
    analytical Hessian (\code{hessA = TRUE}) or from the Hessian of the
    \code{optim} output (\code{hessA = FALSE})),}
  \item{dist}{character string describing the type of NB distribution used,}
  \item{link}{character string describing the link of the mean,}
  \item{link.theta}{character string describing the link of the dispersion
    parameter theta,}
  \item{converged}{logical indicating successful convergence of \code{optim},}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}
  \item{levels}{levels of the categorical regressors,}
  \item{contrasts}{a list with elements \code{"mu"} and \code{"theta"}
    containing the contrasts corresponding to \code{levels} from the
    respective parts,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response count vector (if \code{y = TRUE}),}
  \item{x}{the model matrix for the mean (if \code{x = TRUE}),}
  \item{z}{the model matrix for the mean (if \code{z = TRUE}),}
}

\references{
  Cameron AC, Trivedi PK (2013). \emph{Regression Analysis of Count Data}, 2nd ed.
    New York: Cambridge University Press.
}

\seealso{\code{\link{nbreg.control}}, \code{\link[stats]{glm}},
  \code{\link[stats]{glm.fit}}, \code{\link[MASS]{glm.nb}},
}

\examples{
data("CrabSatellites", package = "countreg")

## NB2
fm_nb2 <- nbreg(satellites ~ width + color, data = CrabSatellites)

## NB1
fm_nb1 <- nbreg(satellites ~ width + color, data = CrabSatellites, dist = "NB1")

## NBH
fm_nbh <- nbreg(satellites ~ width + color | weight, data = CrabSatellites)

## NB1 with variable theta
fm_nb1h <- nbreg(satellites ~ width + color | weight, data = CrabSatellites,
                dist = "NB1")

## Example not run:
## data
# data("GSOEP", package = "countreg")
# gsoep <- subset(GSOEP, year == "1984")

## NB2
# fm_nb2 <- nbreg(docvis ~ educ + public + addon,
#                 data = gsoep)
                
## NB1
# fm_nb1 <- nbreg(docvis ~ educ + public + addon,
#                 data = gsoep, dist = "NB1")                

## NBH
# fm_nbh <- nbreg(docvis ~ educ + public + addon | married + public,
#                 data = gsoep)
                
## NB1 with variable theta
# fm_nb1h <- nbreg(docvis ~ educ + public + addon | married + public,
#                 data = gsoep, dist = "NB1")
}

\keyword{regression}
