\name{hpois}
\alias{dhpois}
\alias{phpois}
\alias{qhpois}
\alias{rhpois}
\alias{shpois}

\title{The Hurdle Poisson Distribution}

\description{
  Density, distribution function, quantile function, random
  generation and score function for the (zero-)hurdle Poisson
  distribution with parameters \code{mu} (= mean of the
  underlying Poisson distribution)and hurdle crossing probability
  \code{pi} (i.e., \code{1 - pi} is the probability for observed zeros).
}

\usage{
dhpois(x, lambda, pi, log = FALSE)
phpois(q, lambda, pi, lower.tail = TRUE, log.p = FALSE)
qhpois(p, lambda, pi, lower.tail = TRUE, log.p = FALSE)
rhpois(n, lambda, pi)
shpois(x, lambda, pi, parameter = c("lambda", "pi"), drop = TRUE)
}

\arguments{
  \item{x}{vector of (positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{lambda}{vector of non-negative means of the underlying
    Poisson distribution.}
  \item{pi}{vector of hurdle crossing probabilities (i.e., \code{1 - pi}
    is the probability for observed zeros).}
  \item{log, log.p}{logical. If \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{parameter}{character. Should the derivative with respect to
    \code{"lambda"} and/or \code{"pi"} be computed?}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE})
    or should the dimension be dropped (\code{drop = TRUE}, the default)?}
}

\details{
The underlying Poisson distribution has density
\deqn{f(x) = \frac{\lambda^x e^{-\lambda}}{x!}}{f(x) = \lambda^x exp(-\lambda)/x!}
for \eqn{x = 0, 1, 2, \ldots}. The hurdle density is then simply obtained as
\deqn{g(x) = \pi \cdot \frac{f(x)}{1 - f(0)}}{g(x) = \pi f(x)/(1 - f(0))}
for \eqn{x = 1, 2, \ldots} and \eqn{g(0) = 1 - \pi}, respectively.
}

\value{
  \code{dhpois} gives the (log) density,
  \code{phpois} gives the (log) distribution function,
  \code{qhpois} gives the quantile function,
  \code{rhpois} generates random deviates, and
  \code{shpois} gives the score function (= derivative of
  the log-density with respect to mu and/or theta and/or pi).
}


\seealso{\code{\link[stats]{dpois}}, \code{\link{hurdle}}
}

\keyword{distribution}
