\name{NegBinomial-extensions}
\alias{snbinom}
\alias{hnbinom}
\alias{mean_nbinom}
\alias{var_nbinom}
\title{
Extension of the Negative Binomial Distribution
}
\description{
Score function, hessian, mean, and variance
for the negative binomial distribution
with parameters \code{mu} and \code{size}.
}
\usage{
snbinom(x, mu, size, parameter = c("mu", "size"), drop = TRUE)
hnbinom(x, mu, size, parameter = c("mu", "size"), drop = TRUE)
mean_nbinom(mu, size, drop = TRUE)
var_nbinom(mu, size, drop = TRUE)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{mu}{_mean_ of distribution.}
  \item{size}{dispersion parameter. Must be strictly positive.}
  \item{parameter}{character. Derivatives are computed wrt this paramter.}
  \item{drop}{logical. Should the result be a matrix (\code{drop = FALSE}) or
should the dimension be dropped (\code{drop = TRUE}, the default)?}
}
\details{
The negative binomial with \code{mu} and \code{size}
(or \code{theta}) has density
\deqn{
f(y | \mu, \theta) = \frac{\Gamma(\theta + y)}{\Gamma({\theta}) \cdot y!} \cdot
 \frac{\mu^y \cdot \theta^\theta}{(\mu + \theta)^{\theta + y}}, \quad y \in \{0, 1, 2, \dots\}
}{
f(y | mu, theta) = Gamma(theta + y) / Gamma(theta) * y! *
 mu^y * theta^theta / (mu + theta)^(\theta + y), y = 0, 1, 2, ... 
}

Derivatives of the log-likelihood \eqn{\ell}{loglik} wrt \eqn{\mu}{mu}:
\deqn{
\frac{\partial \ell}{\partial \mu} = \frac{y}{\mu} - \frac{y + \theta}{\mu + \theta}
}{
d loglik / d mu = y/mu - (y + theta)/(mu + theta)
}

\deqn{
\frac{\partial^2 \ell}{\partial \mu^2} = - \frac{y}{\mu^2}
 + \frac{y + \theta}{(\mu + \theta)^2}
}{
d^2 loglik / d mu^2 = - y / mu^2
 + (y + theta)/(mu + theta)^2
}

Derivatives wrt \eqn{\theta}{theta}:
\deqn{
\frac{\partial \ell}{\partial \theta} = \psi_0(y + \theta) - \psi_0(\theta)
 + \log(\theta) + 1 - \log(\mu + \theta) - \frac{y + \theta}{\mu + \theta}
}{
d loglik / d theta = psi0(y + theta) - psi0(theta)
 + log(theta) + 1 - log(mu + theta) - (y + theta)/(mu + theta)
}

\deqn{
\frac{\partial^2 \ell}{\partial \theta^2} = \psi_1(y + \theta) - \psi_1(\theta)
 + \frac{1}{\theta} - \frac{2}{\mu + \theta} + \frac{y + \theta}{(\mu + \theta)^2}
}{
d^2 loglik / d theta^2 = psi1(y + theta) - psi1(theta)
 + 1/theta - 2/(mu + theta) + (y + theta)/(mu + theta)^2
}

\eqn{\psi_0}{psi0} and \eqn{\psi_1}{psi1} denote the digamma and
trigamma function, respectively.

The derivative wrt \eqn{\mu}{mu} and \eqn{\theta}{theta}:
\deqn{
\frac{\partial^2 \ell}{\partial\mu\partial\theta} =
 = \frac{y - \mu}{(\mu + \theta)^2}
}{
d^2 loglik / d mu d theta
 = (y - mu)/(\mu + \theta)^2
}

}
\value{
\code{snbinom} gives the score function, i.e., the 1st
derivative of the log-density wrt mu or theta and
\code{hnbinom} gives the hessian, i.e., the 2nd
derivative of the log-density wrt mu and/or theta.
\code{mean} and \code{var} give the mean and
variance, respectively.
}
\note{
No parameter \code{prob}---as in \code{\link{dnbinom}}, \code{\link{pnbinom}},
\code{\link{qnbinom}} and \code{\link{rnbinom}}---is implemented in the 
functions \code{snbinom} and \code{hnbinom}.
}
\seealso{
\link{NegBinomial} encompassing \code{\link{dnbinom}}, \code{\link{pnbinom}},
\code{\link{qnbinom}} and \code{\link{rnbinom}}.
}
\examples{
## Simulate some data
set.seed(123)
y <- rnbinom(1000, size = 2, mu = 2)

## Plot log-likelihood function
par(mfrow = c(1, 3))
ll <- function(x) {sum(dnbinom(y, size = x, mu = 2, log = TRUE))}
curve(sapply(x, ll), 1, 4, xlab = expression(theta), ylab = "",
      main = "Log-likelihood")
abline(v = 2, lty = 3)

## Plot score function
curve(sapply(x, function(x) sum(snbinom(y, size = x, mu = 2, parameter = "size"))),
      1, 4, xlab = expression(theta), ylab = "", main = "Score")
abline(h = 0, lty = 3)
abline(v = 2, lty = 3)

## Plot hessian
curve(sapply(x, function(x) sum(hnbinom(y, size = x, mu = 2, parameter = "size"))),
      1, 4, xlab = expression(theta), ylab = "", main = "Hessian")
abline(v = 2, lty = 3)
}
\keyword{distribution}
