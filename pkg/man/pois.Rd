\name{Poisson-extensions}
\alias{spois}
\alias{hpois}
\alias{mean_pois}
\alias{var_pois}
\title{
Extension of the Poisson Distribution
}
\description{
Score function, hessian, mean, and variance
for the Poisson distribution
with parameter \code{lambda}.
}
\usage{
spois(x, lambda, parameter = "lambda", drop = TRUE)
hpois(x, lambda, parameter = "lambda", drop = TRUE)
mean_pois(lambda, drop =TRUE)
var_pois(lambda, drop =TRUE)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{lambda}{vector of (non-negative) means.}
  \item{parameter}{character. Derivatives are computed wrt this
paramter. Note: Only \code{"lambda"} is implemented.}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE}) or
should the dimension be dropped (\code{drop = TRUE}, the default)?}
}
\details{
The Poisson distribution has density
  \deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{x!}}{p(x) = \lambda^x exp(-\lambda)/x!}
  for \eqn{x = 0, 1, 2, \ldots} .

The score function is
  \deqn{s(\lambda) = \frac{x}{\lambda} - 1}{
  s(lambda) = x/lambda - 1}

The hessian is
  \deqn{h(\lambda) = - \frac{x}{\lambda^2}}{
  h(lambda) = - x/lambda^2}
}
\value{
\code{spois} gives the score function, i.e., the 1st
derivative of the log-density wrt lambda and
\code{hpois} gives the hessian, i.e., the 2nd
derivative of the log-density wrt lambda.
\code{mean} and \code{var} give the mean and
variance, respectively.
}
\seealso{
\link{Poisson} encompassing \code{\link{dpois}}, \code{\link{ppois}},
\code{\link{qpois}} and \code{\link{rpois}}.    
}
\examples{
## Simulate some data
set.seed(123)
y <- rpois(50, lambda = 3)

## Plot log-likelihood function
par(mfrow = c(1,3))
ll <- function(x) {sum(dpois(y, x, log = TRUE))}
curve(sapply(x, ll), 1, 5, xlab = expression(lambda), ylab = "",
      main = "Log-likelihood")
abline(v = 3, lty = 3)

## Plot score function
curve(sapply(x, function(x) sum(spois(y, x))), 1, 5,
      xlab = expression(lambda), ylab = "", main = "Score")
abline(h = 0, lty = 3)
abline(v = 3, lty = 3)

## Plot hessian
curve( sapply(x, function(x) sum(hpois(y, x))), 1, 5,
      xlab = expression(lambda), ylab = "", main = "Hessian")
abline(v = 3, lty = 3)
}
\keyword{distribution}
