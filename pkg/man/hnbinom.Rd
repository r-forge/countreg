\name{hnbinom-extensions}
\alias{shnbinom}
\alias{hhnbinom}
\alias{mean_hnbinom}
\alias{var_hnbinom}

\title{Extension of the Hurdle Negative Binomial Distribution}

\description{
  Score function, hessian, mean, and, variance
  for the (zero-)hurdle negative binomial
  distribution with parameters \code{mu} (= mean of the
  underlying negative binomial distribution), dispersion parameter \code{theta}
  (or equivalently \code{size}), and hurdle crossing probability
  \code{pi} (i.e., \code{1 - pi} is the probability for observed zeros).
}

\usage{
shnbinom(x, mu, theta, size, pi, parameter = c("mu", "theta", "pi"), drop = TRUE)
hhnbinom(x, mu, theta, size, pi, parameter = c("mu", "theta", "pi"), drop = TRUE)
mean_hnbinom(mu, theta, size, pi, drop = TRUE)
var_hnbinom(mu, theta, size, pi, drop = TRUE)
}

\arguments{
  \item{x}{vector of (positive integer) quantiles.}
  \item{mu}{vector of non-negative means of the underlying
    negative binomial distribution.}
  \item{theta, size}{vector of strictly positive dispersion
    parameters (shape parameter of the gamma mixing distribution).
    Only one of \code{theta} or \code{size} must be specified.}
  \item{pi}{vector of hurdle crossing probabilities (i.e., \code{1 - pi}
    is the probability for observed zeros).}
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
  \code{shnbinom} gives the score function (= derivative of
  the log-density with respect to mu and/or theta and/or pi).
  \code{hhnbinom} gives the hessian (= 2nd derivative of
  the log-density with respect to mu and/or theta and/or pi).
  \code{mean_hnbinom} and \code{var_hnbinom} give the mean
  and the variance, respectively.
}


\seealso{\code{\link[distributions3]{dhnbinom}}, \code{\link[stats]{dnbinom}}, \code{\link{hurdle}}
}

\keyword{distribution}
