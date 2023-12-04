\name{hpois-extensions}
\alias{shpois}
\alias{hhpois}
\alias{mean_hpois}
\alias{var_hpois}

\title{Extension of the Hurdle Poisson Distribution}

\description{
  Score function, hessian, mean, and, variance
  for the (zero-)hurdle Poisson
  distribution with parameters \code{mu} (= mean of the
  underlying Poisson distribution) and hurdle crossing probability
  \code{pi} (i.e., \code{1 - pi} is the probability for observed zeros).
}

\usage{
shpois(x, lambda, pi, parameter = c("lambda", "pi"), drop = TRUE)
hhpois(x, lambda, pi, parameter = c("lambda", "pi"), drop = TRUE)
mean_hpois(lambda, pi, drop = TRUE)
var_hpois(lambda, pi, drop = TRUE)
}

\arguments{
  \item{x}{vector of (positive integer) quantiles.}
  \item{lambda}{vector of non-negative means of the underlying
    Poisson distribution.}
  \item{pi}{vector of hurdle crossing probabilities (i.e., \code{1 - pi}
    is the probability for observed zeros).}
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
  \code{shpois} gives the score function (= derivative of
  the log-density with respect to lambda and/or pi).
  \code{hhpois} gives the hessian (= 2nd derivative of
  the log-density with respect to lambda and/or pi).
  \code{mean_hpois} and \code{var_hpois} give the mean and
  the variance, respectively.
}


\seealso{\code{\link[distributions3]{dhpois}}, \code{\link[stats]{dpois}}, \code{\link{hurdle}}
}

\keyword{distribution}
