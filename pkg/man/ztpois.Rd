\name{Zero-Truncated-Poisson}
\alias{dztpois}
\alias{pztpois}
\alias{qztpois}
\alias{rztpois}
\alias{sztpois}
\alias{hztpois}
\alias{mean_ztpois}
\alias{var_ztpois}

\title{The Zero-Truncated Poisson Distribution}

\description{
  Density, distribution function, quantile function, random
  generation, score function, hessian, mean, and variance
  for the zero-truncated Poisson
  distribution with parameter \code{lambda} (= mean of the
  untruncated distribution) or \code{mean} (= of the truncated
  distribution).
}

\usage{
dztpois(x, lambda, mean, log = FALSE)
pztpois(q, lambda, mean, lower.tail = TRUE, log.p = FALSE)
qztpois(p, lambda, mean, lower.tail = TRUE, log.p = FALSE)
rztpois(n, lambda, mean)
sztpois(x, lambda, mean, parameter = "lambda", drop = TRUE)
hztpois(x, lambda, mean, parameter = "lambda", drop = TRUE)
mean_ztpois(lambda, mean, drop = TRUE)
var_ztpois(lambda, mean, drop = TRUE)
}

\arguments{
  \item{x}{vector of (positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{lambda}{vector of (non-negative) means of the untruncated
    Poisson distribution. Only one of \code{lambda} or \code{mean}
    should be specified.}
  \item{mean}{vector of means (greater than 1) of the zero-truncated
    Poisson distribution. Only one of \code{lambda} or \code{mean}
    should be specified.}
  \item{log, log.p}{logical. If \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{parameter}{character. Should the derivative with respect to
    \code{"lambda"} or \code{"mean"} be computed?}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE})
    or should the dimension be dropped (\code{drop = TRUE}, the default)?}
}

\details{
The untruncted Poisson distribution has density
\deqn{f(x) = \frac{\lambda^x e^{-\lambda}}{x!}}{f(x) = \lambda^x exp(-\lambda)/x!}
for \eqn{x = 0, 1, 2, \ldots}. The zero-truncted density is then
simply obtained as
\deqn{g(x) = \frac{f(x)}{1 - f(0)}}{g(x) = f(x)/(1 - f(0))}
for \eqn{x = 1, 2, \ldots}.

The zero-truncated distribution has expectation
\eqn{E(X) = \mu = \lambda / (1 - \exp(-\lambda))} and variance
\eqn{Var(X) = \mu \cdot (\lambda + 1 - \mu)}, where \eqn{\lambda}
is the expectation of the untruncated Poisson distribution.

Despite the simple form of the transformation \eqn{\mu(\lambda)} the
inverse \eqn{\lambda(\mu)} has no closed-form solution and is computed
numerically if needed.
}

\value{
  \code{dztpois} gives the (log) density,
  \code{pztpois} gives the (log) distribution function,
  \code{qztpois} gives the quantile function,
  \code{rztpois} generates random deviates, and
  \code{sztpois} gives the score function (= derivative of
  the log-density with respect to lambda or mean).
  \code{hztpois} gives the hessian (= 2nd derivative of
  the log-density with respect to lambda or mean).
  \code{mean_ztpois} and \code{var_ztpois} give the mean
  and the variance, respectively.
}


\seealso{\code{\link{ztpoisson}}, \code{\link[stats]{dpois}}, \code{\link{zerotrunc}}
}

\keyword{distribution}
