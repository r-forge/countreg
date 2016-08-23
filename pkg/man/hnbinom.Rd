\name{hnbinom}
\alias{dhnbinom}
\alias{phnbinom}
\alias{qhnbinom}
\alias{rhnbinom}
\alias{shnbinom}

\title{The Hurdle Negative Binomial Distribution}

\description{
  Density, distribution function, quantile function, random
  generation and score function for the (zero-)hurdle negative binomial
  distribution with parameters \code{mu} (= mean of the
  underlying negative binomial distribution), dispersion parameter \code{theta}
  (or equivalently \code{size}), and hurdle crossing probability
  \code{pi} (i.e., \code{1 - pi} is the probability for observed zeros).
}

\usage{
dhnbinom(x, mu, theta, size, pi, log = FALSE)
phnbinom(q, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE)
qhnbinom(p, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE)
rhnbinom(n, mu, theta, size, pi)
shnbinom(x, mu, theta, size, pi, parameter = c("mu", "theta", "pi"), drop = TRUE)
}

\arguments{
  \item{x}{vector of (positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{mu}{vector of non-negative means of the underlying
    negative binomial distribution.}
  \item{theta, size}{vector of strictly positive dispersion
    parameters (shape parameter of the gamma mixing distribution).
    Only one of \code{theta} or \code{size} must be specified.}
  \item{pi}{vector of hurdle crossing probabilities (i.e., \code{1 - pi}
    is the probability for observed zeros).}
  \item{log, log.p}{logical. If \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{parameter}{character. Should the derivative with respect to
    \code{"mu"} and/or \code{"theta"} and/or \code{"pi"} be computed?}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE})
    or should the dimension be dropped (\code{drop = TRUE}, the default)?}
}

\details{
The underlying negative binomial distribution has density
\deqn{f(x) =
  \frac{\Gamma(x + \theta)}{\Gamma(\theta) x!} \cdot \frac{\mu^y \theta^\theta}{(\mu + \theta)^{y + \theta}}}{
  \Gamma(x + \theta)/(\Gamma(\theta) x!) (\mu^y \theta^\theta)/((\mu + \theta)^(y + \theta)}
for \eqn{x = 0, 1, 2, \ldots}. The hurdle density is then simply obtained as
\deqn{g(x) = \pi \cdot \frac{f(x)}{1 - f(0)}}{g(x) = \pi f(x)/(1 - f(0))}
for \eqn{x = 1, 2, \ldots} and \eqn{g(0) = 1 - \pi}, respectively.
}

\value{
  \code{dhnbinom} gives the (log) density,
  \code{phnbinom} gives the (log) distribution function,
  \code{qhnbinom} gives the quantile function,
  \code{rhnbinom} generates random deviates, and
  \code{shnbinom} gives the score function (= derivative of
  the log-density with respect to mu and/or theta and/or pi).
}


\seealso{\code{\link[stats]{dnbinom}}, \code{\link{hurdle}}
}

\keyword{distribution}
