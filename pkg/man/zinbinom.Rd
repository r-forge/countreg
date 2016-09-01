\name{zinbinom}
\alias{dzinbinom}
\alias{pzinbinom}
\alias{qzinbinom}
\alias{rzinbinom}
\alias{szinbinom}

\title{The Zero-Inflated Negative Binomial Distribution}

\description{
  Density, distribution function, quantile function, random
  generation and score function for the zero-inflated negative binomial
  distribution with parameters \code{mu} (= mean of the
  uninflated distribution), dispersion parameter \code{theta}
  (or equivalently \code{size}), and inflation probability
  \code{pi} (for structural zeros).
}

\usage{
dzinbinom(x, mu, theta, size, pi, log = FALSE)
pzinbinom(q, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE)
qzinbinom(p, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE)
rzinbinom(n, mu, theta, size, pi)
szinbinom(x, mu, theta, size, pi, parameter = c("mu", "theta", "pi"), drop = TRUE)
}

\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{mu}{vector of non-negative means of the uninflated
    negative binomial distribution.}
  \item{theta, size}{vector of strictly positive dispersion
    parameters (shape parameter of the gamma mixing distribution).
    Only one of \code{theta} or \code{size} must be specified.}
  \item{pi}{vector of zero inflation probabilities for structural
    zeros.}
  \item{log, log.p}{logical. If \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{parameter}{character. Should the derivative with respect to
    \code{"mu"} and/or \code{"theta"} and/or \code{"pi"} be computed?}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE})
    or should the dimension be dropped (\code{drop = TRUE}, the default)?}
}

\details{
The uninflated negative binomial distribution has density
\deqn{f(x) =
  \frac{\Gamma(x + \theta)}{\Gamma(\theta) x!} \cdot \frac{\mu^y \theta^\theta}{(\mu + \theta)^{y + \theta}}}{
  \Gamma(x + \theta)/(\Gamma(\theta) x!) (\mu^y \theta^\theta)/((\mu + \theta)^(y + \theta)}
for \eqn{x = 0, 1, 2, \ldots}. The zero-inflated density is then simply obtained as
\deqn{g(x) = \pi \cdot I_{\{0\}}(x) + (1 - \pi) \cdot f(x)}{g(x) = \pi I(x = 0) + (1 - \pi) f(x)}
where \eqn{I} is the indicator function (for the point mass at zero).
}

\value{
  \code{dzinbinom} gives the (log) density,
  \code{pzinbinom} gives the (log) distribution function,
  \code{qzinbinom} gives the quantile function,
  \code{rzinbinom} generates random deviates, and
  \code{szinbinom} gives the score function (= derivative of
  the log-density with respect to mu and/or theta and/or pi).
}


\seealso{\code{\link[stats]{dnbinom}}, \code{\link{zeroinfl}}
}

\keyword{distribution}
