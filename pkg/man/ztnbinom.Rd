\name{ztnbinom-extensions}
\alias{sztnbinom}
\alias{hztnbinom}
\alias{mean_ztnbinom}
\alias{var_ztnbinom}

\title{Extension of the Zero-Truncated Negative Binomial Distribution}

\description{
  Score function, hessian, mean, and variance
  for the zero-truncated negative binomial
  distribution with parameters \code{mu} (= mean of the
  untruncated distribution) and dispersion parameter \code{theta}
  (or equivalently \code{size}).
}

\usage{
sztnbinom(x, mu, theta, size, parameter = c("mu", "theta", "size"), drop = TRUE)
hztnbinom(x, mu, theta, size, parameter = c("mu", "theta"), drop = TRUE)
mean_ztnbinom(mu, theta, size, drop = TRUE)
var_ztnbinom(mu, theta, size, drop = TRUE)
}

\arguments{
  \item{x}{vector of (positive integer) quantiles.}
  \item{mu}{vector of non-negative means of the untruncated
    negative binomial distribution.}
  \item{theta, size}{vector of strictly positive dispersion
    parameters (shape parameter of the gamma mixing distribution).
    Only one of \code{theta} or \code{size} must be specified.}
  \item{parameter}{character. Should the derivative with respect to
    \code{"mu"} and/or \code{"theta"}/\code{"size"} be computed?}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE})
    or should the dimension be dropped (\code{drop = TRUE}, the default)?}
}

\details{
The untruncated negative binomial distribution has density
\deqn{f(x) =
  \frac{\Gamma(x + \theta)}{\Gamma(\theta) x!} \cdot \frac{\mu^y \theta^\theta}{(\mu + \theta)^{y + \theta}}}{
  \Gamma(x + \theta)/(\Gamma(\theta) x!) (\mu^y \theta^\theta)/((\mu + \theta)^(y + \theta)}
for \eqn{x = 0, 1, 2, \ldots}. The zero-truncated density is then simply obtained as
\deqn{g(x) = \frac{f(x)}{1 - f(0)}}{g(x) = f(x)/(1 - f(0))}
for \eqn{x = 1, 2, \ldots}.
}

\value{
  \code{sztnbinom} gives the score function (= derivative of
  the log-density with respect to mu and/or theta).
  \code{hztnbinom} gives the hessian (= 2nd derivative of
  the log-density with respect to mu and/or theta).
  \code{mean_ztnbinom} and \code{var_ztnbinom} give the
  mean and the variance, respectively.
}


\seealso{\code{\link[distributions3]{dztnbinom}}, \code{\link[stats]{dnbinom}}, \code{\link{zerotrunc}}
}

\keyword{distribution}
