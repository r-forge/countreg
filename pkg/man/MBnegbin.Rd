\name{MBnegbin}
\alias{MBbinomial}
\alias{MBnegbin}
\alias{MBztpoisson}
\alias{MBztnegbin}

\title{mboost Families for Binary, (Zero-Truncated) Negative Binomial and Zero-Truncated Poisson Regression}

\description{
Family generators for model-based boosting of count data regressions
using mboost.
}

\usage{
MBbinomial(link = "logit")

MBnegbin(theta = NULL, link = "log",
  control = list(reltol = .Machine$double.eps^(1/1.5), maxit = 500))

MBztpoisson <- function(link = "log",
  control = list(reltol = .Machine$double.eps^(1/1.5), maxit = 500))

MBztnegbin(theta = NULL, link = "log",
  control = list(reltol = .Machine$double.eps^(1/1.5), maxit = 500))
}

\arguments{
  \item{link}{character or object of class \code{"link-glm"} for
    the link function linking the expectation and the predictor.}
  \item{theta}{numeric or \code{NULL}. Value of the theta parameter
    of the negative binomial model. If \code{NULL}, theta is estimated
    along with the regression coefficients.}
  \item{control}{list with control parameters passed to
    \code{\link[stats]{optim}}.}
}

\details{
The family generators \code{MBbinomial}, \code{MBnegbin}, \code{MBztpoisson}, \code{MBztnegbin} enable
boosting of binary regressions, negative binomial count regressions, zero-truncated
Poisson count regressions, and zero-truncated negative binomial count regressions, respectively. 
Family \code{MBbinomial} is comparable
to \code{\link[mboost]{Binomial}} but supports any link function (not just logit and probit).
Family \code{MBnegbin} is comparable to \code{\link[mboost]{NBinomial}} but is typically
much faster because the nuisance parameter \code{theta} is estimated using analytical
gradients (via \code{\link[stats]{optim}}) and setting better starting values.
\code{MBztpoisson} and \code{MBztnegbin} enable zero-truncated Poisson and negative binomial regressions so that
also the count parts of hurdle models can be easily estimated.
}

\value{
An object of class \code{boost_family_glm}.
}

\seealso{
  \code{\link[mboost]{mboost}},
  \code{\link[mboost]{glmboost}}, 
  \code{\link[mboost]{gamboost}},
  \code{\link[mboost]{Binomial}},
  \code{\link[mboost]{Poisson}},
  \code{\link[mboost]{NBinomial}}
}

\examples{
### Negative binomial regression for CrabSatellites ----------------------------

if(require("mboost")) {
## crab satellite data using ordered factors as numeric
data("CrabSatellites", package = "countreg")
CrabSatellites <- transform(CrabSatellites,
  color = as.numeric(color),
  spine = as.numeric(spine)
)  

## comparison of ML and boosting with NBinomial() vs. MBnegbin()
system.time(m0 <- glm.nb(satellites ~ width + color, data = CrabSatellites))
system.time(m1 <- glmboost(satellites ~ width + color, data = CrabSatellites,
  family = NBinomial(), control = boost_control(mstop = 500)))
system.time(m2 <- glmboost(satellites ~ width + color, data = CrabSatellites,
  family = MBnegbin(), control = boost_control(mstop = 500)))
## note that mstop is _not_ tuned here to (ab)use mboost to get the ML estimator

## compare coefficients
cbind(c(coef(m0), "theta" = m0$theta),
  c(coef(m1, off2int = TRUE, which = ""), nuisance(m1)),
  c(coef(m1, off2int = TRUE, which = ""), nuisance(m1))
)
}

### Hurdle regression for CrabSatellites using spline terms --------------------

\donttest{
if(require("mboost")) {
## ML estimation
g <- hurdle(satellites ~ width + color, data = CrabSatellites, dist = "negbin")
summary(g)

## boosting of zero hurdle
g0 <- gamboost(factor(satellites > 0) ~ bbs(width) + bbs(color),
  data = CrabSatellites, family = MBbinomial())
set.seed(0)
g0cv <- cvrisk(g0)
g0[mstop(g0cv)]

## boosting of count regression
g1 <- gamboost(satellites ~ bbs(width) + bbs(color),
  data = subset(CrabSatellites, satellites > 0), family = MBztnegbin())
set.seed(1)
g1cv <- cvrisk(g1)
g1[mstop(g1cv)]

par(mfrow = c(1, 2))

## optimal mstop values
plot(g0cv)
plot(g1cv)
## -> no effects in covariates for count part

## partial effects in zero hurdle
plot(g0)
## -> large effect of width, moderate effect of color with
## both effects almost linear
}
}

### Hurdle regression for RecreationDemand using linear terms ------------------

\dontrun{
library("mboost")
data("RecreationDemand", package = "AER")

### Zero hurdle ##

## ML vs. boosting
z0 <- glm(factor(trips > 0) ~ ., data = RecreationDemand, family = binomial)
z1 <- glmboost(factor(trips > 0) ~ ., data = RecreationDemand, family = MBbinomial(),
  control = boost_control(mstop = 5000))
plot(z1)

## tune mstop
set.seed(0)
z1cv <- cvrisk(z1)
z1cv
plot(z1cv)
## very flat (presumably due to separation?)
## -> stop earlier manually
z1[3000]

## compare coefficients
cbind(coef(z0), coef(z1, off2int = TRUE, which = ""))
## -> some shrunken entirely to zero,
## coefficient of variable with separation (userfee) shrunken considerably


### Count (zero-truncated)

## ML and boosting count part
c0 <- zerotrunc(trips ~ ., data = subset(RecreationDemand, trips > 0), dist = "negbin")
c1 <- glmboost(trips ~ ., data = subset(RecreationDemand, trips > 0),
  family = MBztnegbin(), control = boost_control(mstop = 5000))
plot(c1)

## tune mstop
set.seed(0)
c1cv <- cvrisk(c1)
c1cv
plot(c1cv)

## use mstop from cvrisk
c1[mstop(c1cv)]

## compare coefficients
cbind(c(coef(c0), "theta" = c0$theta),
  c(coef(c1, off2int = TRUE, which = ""), nuisance(c1)))
## -> not much difference
}
}

\keyword{models}
