\name{zipois}
\alias{dzipois}
\alias{pzipois}
\alias{qzipois}
\alias{rzipois}
\alias{szipois}

\title{The Zero-Inflated Poisson Distribution}

\description{
  Density, distribution function, quantile function, random
  generation and score function for the zero-inflated Poisson
  distribution with parameters \code{lambda} (= mean of the
  uninflated distribution) and inflation probability
  \code{pi} (for structural zeros).
}

\usage{
dzipois(x, lambda, pi, log = FALSE)
pzipois(q, lambda, pi, lower.tail = TRUE, log.p = FALSE)
qzipois(p, lambda, pi, lower.tail = TRUE, log.p = FALSE)
rzipois(n, lambda, pi)
szipois(x, lambda, pi, parameter = c("lambda", "pi"), drop = TRUE)
}

\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{lambda}{vector of non-negative means of the uninflated
    Poisson distribution.}
  \item{pi}{vector of zero inflation probabilities for structural
    zeros.}
  \item{log, log.p}{logical. If \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{parameter}{character. Should the derivative with respect to
    \code{"mu"} and/or \code{"size"} be computed?}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE})
    or should the dimension be dropped (\code{drop = TRUE}, the default)?}
}

\details{
The uninflated Poisson distribution has density
\deqn{f(x) = \frac{\lambda^x e^{-\lambda}}{x!}}{f(x) = \lambda^x exp(-\lambda)/x!}
for \eqn{x = 0, 1, 2, \ldots}. The zero-inflated density is then simply obtained as
\deqn{g(x) = \pi \cdot I_{\{0\}}(x) + (1 - \pi) \cdot f(x)}{g(x) = \pi I(x = 0) + (1 - \pi) f(x)}
where \eqn{I} is the indicator function (for the point mass at zero).
}

\value{
  \code{dzipois} gives the (log) density,
  \code{pzipois} gives the (log) distribution function,
  \code{qzipois} gives the quantile function,
  \code{rzipois} generates random deviates, and
  \code{szipois} gives the score function (= derivative of
  the log-density with respect to lambda and/or pi).
}


\seealso{\code{\link[stats]{dpois}}, \code{\link{zeroinfl}}
}

\keyword{distribution}
