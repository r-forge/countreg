\name{Binomial-extensions}
\alias{sbinom}
\alias{hbinom}
\alias{mean_binom}
\alias{var_binom}
\title{
Extension of the Binomial Distribution
}
\description{
Score function, hessian, mean, and variance
for the binomial distribution
with parameters \code{prob} and \code{size}.
}
\usage{
sbinom(x, prob, size, parameter = "prob", drop = TRUE)
hbinom(x, prob, size, parameter = "prob", drop = TRUE)
mean_binom(prob, size, drop = TRUE)
var_binom(prob, size, drop = TRUE)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{prob}{probability of success on each trial.}
  \item{size}{number of trials (zero or more).}
  \item{parameter}{character. Derivatives are computed wrt this
paramter. Note: Only \code{"prob"} is implemented.}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE}) or
should the dimension be dropped (\code{drop = TRUE}, the default)?}
}
\details{
The binomial distribution with \code{size} \eqn{= n} and
\code{prob} \eqn{= p} has density
\deqn{p(x) = {n \choose x} {p}^{x} {(1-p)}^{n-x}}{
  p(x) = choose(n, x) p^x (1-p)^(n-x)}
for \eqn{x = 0, \ldots, n}.

The score function is
\deqn{s(p) = \frac{x}{p} - \frac{n-x}{1-p}}{
  s(p) = x/p - (n-x)/(1-p)}

The hessian is
\deqn{h(p) = - \frac{x}{p^2} - \frac{n-x}{(1-p)^2}}{
  h(p) = - x/p^2 - (n-x)/(1-p)^2}
}
\value{
\code{sbinom} gives the score function, i.e., the 1st
derivative of the log-density wrt prob and
\code{hbinom} gives the hessian, i.e., the 2nd
derivative of the log-density wrt prob.
\code{mean} and \code{var} give the mean and
variance, respectively.
}
\seealso{
\link{Binomial} encompassing \code{\link{dbinom}}, \code{\link{pbinom}},
\code{\link{qbinom}} and \code{\link{rbinom}}.    
}
\examples{
## Simulate some data
set.seed(123)
y <- rbinom(50, size = 1, prob = 0.3)

## Plot log-likelihood function
par(mfrow = c(1,3))
ll <- function(x) {sum(dbinom(y, size = 1, prob = x, log = TRUE))}
curve(sapply(x, ll), xlab = expression(pi), ylab = "", main = "Log-likelihood")
abline(v = 0.3, lty = 3)

## Plot score function
curve(sapply(x, function(x) sum(sbinom(y, size = 1, x))),
      xlab = expression(pi), ylab = "", main = "Score")
abline(h = 0, lty = 3)
abline(v = 0.3, lty = 3)

## Plot hessian
curve(sapply(x, function(x) sum(hbinom(y, size = 1, x))),
      xlab = expression(pi), ylab = "", main = "Hessian")
abline(v = 0.3, lty = 3)
}
\keyword{distribution}
