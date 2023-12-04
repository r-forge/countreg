\name{zipois-extensions}
\alias{szipois}

\title{Extension of the Zero-Inflated Poisson Distribution}

\description{
  Score function for the zero-inflated Poisson
  distribution with parameters \code{lambda} (= mean of the
  uninflated distribution) and inflation probability
  \code{pi} (for structural zeros).
}

\usage{
szipois(x, lambda, pi, parameter = c("lambda", "pi"), drop = TRUE)
}

\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{lambda}{vector of non-negative means of the uninflated
    Poisson distribution.}
  \item{pi}{vector of zero inflation probabilities for structural
    zeros.}
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
  \code{szipois} gives the score function (= derivative of
  the log-density with respect to lambda and/or pi).
}


\seealso{\code{\link[distributions3]{dzipois}}, \code{\link[stats]{dpois}}, \code{\link{zeroinfl}}
}

\keyword{distribution}
