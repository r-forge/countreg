\name{zerotrunc}
\alias{zerotrunc}
\alias{print.zerotrunc}

\title{Zero-Truncated Count Data Regression}

\description{
  Fit zero-truncated regression models for count data via maximum likelihood. 
}

\usage{
zerotrunc(formula, data, subset, na.action, weights, offset,
  dist = c("poisson", "negbin", "geometric"), theta = Inf,
  control = zerotrunc.control(\dots),
  model = TRUE, y = TRUE, x = FALSE, \dots)
}

\arguments{
  \item{formula}{symbolic description of the model.}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor.}
  \item{dist}{character specification of the count distribution family.}
  \item{theta}{numeric. Alternative (and more flexible) specification of the
    count distribution family. Some values correspond to \code{dist} values:
    \code{theta = Inf} (\code{"poisson"}), \code{theta = 1}
    (\code{"geometric"}), \code{theta = NULL} (\code{"negbin"}). But every
    non-negative value for \code{theta} is allowed. When \code{theta} is
    given, \code{dist} must not be specified and vice versa.}
  \item{control}{a list of control arguments specified via
    \code{\link{zerotrunc.control}}.}
  \item{model, y, x}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.}
  \item{\dots}{arguments passed to \code{\link{zerotrunc.control}} in the
    default setup.}
}

\details{
  Hurdle count models are two-component models with a truncated count
  component for positive counts and a zerotrunc component that models the
  zero counts. Thus, unlike zero-inflation models, there are \emph{not} two
  sources of zeros: the count model is only employed if the zerotrunc for
  modeling the occurence of zeros is exceeded. The count model is typically
  a truncated Poisson or negative binomial regression (with log link).
  The geometric distribution is a special case of the negative binomial with
  size parameter equal to 1. For modeling the zerotrunc (occurence of positive
  counts) either a binomial model can be employed or a censored count distribution.
  Binomial logit and censored geometric models as the zerotrunc part both lead to 
  the same likelihood function and thus to the same coefficient estimates.
  
  All parameters are estimated by maximum likelihood using \code{\link[stats]{optim}},
  with control options set in \code{\link{zerotrunc.control}}.
  Starting values can be supplied, otherwise they are estimated by \code{\link[stats]{glm.fit}}
  (the default). Standard errors are derived numerically using
  the Hessian matrix returned by \code{\link[stats]{optim}}. See
  \code{\link{zerotrunc.control}} for details.
  
  The returned fitted model object is of class \code{"zerotrunc"} and is similar
  to fitted \code{"glm"} objects.
  
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"zerotrunc"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{fitted}}, \code{\link[stats]{terms}},
  \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}}. See
  \code{\link{predict.zerotrunc}} for more details on all methods.
}

\value{
  An object of class \code{"zerotrunc"}, i.e., a list with components including
  \item{coefficients}{estimated coefficients,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a vector of fitted means,}
  \item{optim}{a list with the output from the \code{optim} call for
    minimizing the negative log-likelihood,}
  \item{control}{the control arguments passed to the \code{optim} call,}
  \item{start}{the starting values for the parameters passed to the \code{optim} call(s),}
  \item{weights}{the case weights used (if any),}
  \item{offset}{the offset vector used (if any),}
  \item{n}{number of observations,}
  \item{df.null}{residual degrees of freedom for the null model,}
  \item{df.residual}{residual degrees of freedom for fitted model,}
  \item{terms}{terms objects for the model,}
  \item{theta}{(estimated) \eqn{\theta}{theta} parameter of the
    negative binomial model,}
  \item{SE.logtheta}{standard error for \eqn{\log(\theta)}{log(theta)},}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{vcov}{covariance matrix of the coefficients in the model (derived from the
    Hessian of the \code{optim} output),}
  \item{dist}{character describing the distribution used,}
  \item{converged}{logical indicating successful convergence of \code{optim},}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}
  \item{levels}{levels of the categorical regressors,}
  \item{contrasts}{contrasts corresponding to \code{levels} from the model,}
  \item{model}{the model frame (if \code{model = TRUE}),}
  \item{y}{the response count vector (if \code{y = TRUE}),}
  \item{x}{model matrix (if \code{x = TRUE}).}
}

\references{
Cameron AC, Trivedi PK (1998). \emph{Regression Analysis of Count Data}.
  New York: Cambridge University Press.

Cameron AC, Trivedi PK (2005). \emph{Microeconometrics: Methods and Applications}.
  Cambridge: Cambridge University Press.

Zeileis A, Kleiber C, Jackman S (2008).
  \dQuote{Regression Models for Count Data in R.} 
  \emph{Journal of Statistical Software}, \bold{27}(8), 1--25.
  URL \url{http://www.jstatsoft.org/v27/i08/}.
}

\seealso{\code{\link{zerotrunc.control}}, \code{\link[stats]{glm}},
  \code{\link[stats]{glm.fit}}, \code{\link[MASS]{glm.nb}},
  \code{\link{zeroinfl}}, \code{\link{hurdle}}
}

\examples{
## data
data("CrabSatellites", package = "countreg")
cs <- CrabSatellites[, c("satellites", "width", "color")]
cs$color <- as.numeric(cs$color)
cs <- subset(cs, subset = satellites > 0)

## poisson
zt_p <- zerotrunc(satellites ~ ., data = cs)
## or equivalently
zt_p <- zerotrunc(satellites ~ ., data = cs, theta = Inf)
summary(zt_p)

## negbin
zt_nb <- zerotrunc(satellites ~ ., data = cs, dist = "negbin")
## or equivalently
zt_nb <- zerotrunc(satellites ~ ., data = cs, theta = NULL)
summary(zt_nb)
}

\keyword{regression}
