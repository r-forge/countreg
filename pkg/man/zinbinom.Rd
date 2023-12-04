\name{zinbinom-extensions}
\alias{szinbinom}

\title{Extension of the Zero-Inflated Negative Binomial Distribution}

\description{
  Score function for the zero-inflated negative binomial
  distribution with parameters \code{mu} (= mean of the
  uninflated distribution), dispersion parameter \code{theta}
  (or equivalently \code{size}), and inflation probability
  \code{pi} (for structural zeros).
}

\usage{
szinbinom(x, mu, theta, size, pi, parameter = c("mu", "theta", "pi"), drop = TRUE)
}

\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{mu}{vector of non-negative means of the uninflated
    negative binomial distribution.}
  \item{theta, size}{vector of strictly positive dispersion
    parameters (shape parameter of the gamma mixing distribution).
    Only one of \code{theta} or \code{size} must be specified.}
  \item{pi}{vector of zero inflation probabilities for structural
    zeros.}
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
  \code{szinbinom} gives the score function (= derivative of
  the log-density with respect to mu and/or theta and/or pi).
}


\seealso{\code{\link[distributions3]{dzinbinom}}, \code{\link[stats]{dnbinom}}, \code{\link{zeroinfl}}
}

\keyword{distribution}
