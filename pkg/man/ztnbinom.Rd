\name{ztnbinom}
\alias{dztnbinom}
\alias{pztnbinom}
\alias{qztnbinom}
\alias{rztnbinom}
\alias{sztnbinom}
\alias{hztnbinom}

\title{The Zero-Truncated Negative Binomial Distribution}

\description{
  Density, distribution function, quantile function, random
  generation and score function for the zero-truncated negative binomial
  distribution with parameters \code{mu} (= mean of the
  untruncated distribution) and dispersion parameter \code{theta}
  (or equivalently \code{size}).
}

\usage{
dztnbinom(x, mu, theta, size, log = FALSE)
pztnbinom(q, mu, theta, size, lower.tail = TRUE, log.p = FALSE)
qztnbinom(p, mu, theta, size, lower.tail = TRUE, log.p = FALSE)
rztnbinom(n, mu, theta, size)
sztnbinom(x, mu, theta, size, parameter = c("mu", "theta"), drop = TRUE)
hztnbinom(x, mu, theta, size, parameter = c("mu", "theta"), drop = TRUE)
}

\arguments{
  \item{x}{vector of (positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{mu}{vector of non-negative means of the untruncated
    negative binomial distribution.}
  \item{theta, size}{vector of strictly positive dispersion
    parameters (shape parameter of the gamma mixing distribution).
    Only one of \code{theta} or \code{size} must be specified.}
  \item{log, log.p}{logical. If \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{parameter}{character. Should the derivative with respect to
    \code{"mu"} and/or \code{"size"} be computed?}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE})
    or should the dimension be dropped (\code{drop = TRUE}, the default)?}
}

\details{
The untruncted negative binomial distribution has density
\deqn{f(x) =
  \frac{\Gamma(x + \theta)}{\Gamma(\theta) x!} \cdot \frac{\mu^y \theta^\theta}{(\mu + \theta)^{y + \theta}}}{
  \Gamma(x + \theta)/(\Gamma(\theta) x!) (\mu^y \theta^\theta)/((\mu + \theta)^(y + \theta)}
for \eqn{x = 0, 1, 2, \ldots}. The zero-truncted density is then simply obtained as
\deqn{g(x) = \frac{f(x)}{1 - f(0)}}{g(x) = f(x)/(1 - f(0))}
for \eqn{x = 1, 2, \ldots}.
}

\value{
  \code{dztnbinom} gives the (log) density,
  \code{pztnbinom} gives the (log) distribution function,
  \code{qztnbinom} gives the quantile function,
  \code{rztnbinom} generates random deviates, and
  \code{sztnbinom} gives the score function (= derivative of
  the log-density with respect to mu and/or theta).
}


\seealso{\code{\link[stats]{dnbinom}}, \code{\link{zerotrunc}}
}

\keyword{distribution}
