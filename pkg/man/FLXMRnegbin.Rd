\name{FLXMRnegbin}
\alias{FLXMRnegbin}

\title{FlexMix Interface to Negative Binomial Regression Models}

\description{
FlexMix driver for fitting of negative binomial regression models.
}

\usage{
FLXMRnegbin(formula = . ~ ., theta = NULL, offset = NULL,
  control = list(reltol = .Machine$double.eps^(1/1.5), maxit = 500))
}

\arguments{
  \item{formula}{formula. This is interpreted relative to the formula
    specified in the call to \code{\link[flexmix]{flexmix}} using
    \code{\link{update.formula}}. Default is to use the original
    \code{flexmix} model formula.}
  \item{theta}{numeric or \code{NULL}. Value of the theta parameter
    of the negative binomial model. If \code{NULL}, theta is estimated
    along with the regression coefficients.}
  \item{offset}{numeric. Optional offset vector for the linear
    predictor.}
  \item{control}{list with control parameters passed to
    \code{\link[stats]{optim}}.}
}

\details{
The driver function \code{FLXMRnegbin} enables estimation of finite mixtures
of negative binomial regression models via \code{\link[flexmix]{flexmix}} or
\code{\link[flexmix]{stepFlexmix}}. The driver is modeled after \code{\link[flexmix]{FLXMRglm}}
and supports both fixed and unknown theta. In the M-step for fixed theta,
\code{\link[stats]{glm.fit}} is employed along with the \code{\link[MASS]{negative.binomial}}
family. If the \code{theta} is unknown and has be estimated
along with the regression coefficients, direct optimization using
\code{\link[stats]{optim}} with analytical gradients is employed.
}

\value{
An object of class \code{FLXMRglm}.
}

\seealso{
  \code{\link[flexmix]{flexmix}},
  \code{\link[flexmix]{stepFlexmix}}, 
  \code{\link[flexmix]{FLXMRglm}}, 
  \code{\link[MASS]{negative.binomial}}
}

\examples{
## artificial data from a two-component mixture of geometric regressions
set.seed(1)
d <- data.frame(x = runif(500, -1, 1))
d$cluster <- rep(1:2, each = 250)
d$y <- rnbinom(500, mu = exp(c(1, -1)[d$cluster] + c(0, 3)[d$cluster] * d$x), size = 1)

if(require("flexmix")) {
## fit mixture models with known correct theta and unknown theta
fm1 <- flexmix(y ~ x, data = d, k = 2, model = FLXMRnegbin(theta = 1))
fm0 <- flexmix(y ~ x, data = d, k = 2, model = FLXMRnegbin())

## parameter recovery
parameters(fm1)
parameters(fm0)

## refit to obtain joint summary
summary(refit(fm1, gradient = NULL))
summary(refit(fm0, gradient = NULL))

## refitting both components manually for rootograms
rf1 <- lapply(1:2, function(i)
  glm(y ~ x, data = d, family = negative.binomial(1), weights = posterior(fm1)[,i]))
rf0 <- lapply(1:2, function(i)
  glm.nb(y ~ x, data = d, weights = posterior(fm0)[,i]))

## Rootograms
par(mfrow = c(1, 2))

r11 <- rootogram(rf1[[1]])
r12 <- rootogram(rf1[[2]])

r01 <- rootogram(rf0[[1]])
r02 <- rootogram(rf0[[2]])

rootogram(glm.nb(y ~ x, data = d))
plot(r01 + r02)
}

\dontrun{
## two-component mixture model fro NMES1988 physician office visits
## (fitting takes some time...)
if(require("flexmix") & require("AER")) {

## data from AER
data("NMES1988", package = "AER")
nmes <- NMES1988[, c(1, 7:8, 13, 15, 18:19)] 

## single-component model
nmes_nb <- glm.nb(visits ~ ., data = nmes)

## two-component model
set.seed(1090)
nmes_fnb <- stepFlexmix(visits ~ ., data = nmes, k = 2, model = FLXMRnegbin())

## refit to obtain summary with estimate of joint covariance matrix
summary(refit(nmes_fnb, gradient = NULL))

## refit individual models manually for rootograms
nmes_fnb_rf <- lapply(1:2, function(i)
  glm.nb(visits ~ ., data = nmes, weights = posterior(nmes_fnb)[,i]))

r1 <- rootogram(nmes_fnb_rf[[1]], max = 50, plot = FALSE)
r2 <- rootogram(nmes_fnb_rf[[2]], max = 50, plot = FALSE)

par(mfrow = c(2, 2))
rootogram(nmes_nb, max = 50, main = "Negative Binomial", ylim = c(-1, 25))
plot(r1 + r2, xlab = "visits", main = "Mixture Negative Binomial", ylim = c(-1, 25))
plot(r1, main = "Mixture Negative Binomial (Component 1)", ylim = c(-1, 25))
plot(r2, main = "Mixture Negative Binomial (Component 2)", ylim = c(-1, 25))
}
}
}

\keyword{models}
