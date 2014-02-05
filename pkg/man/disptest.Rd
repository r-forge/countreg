\name{disptest}
\alias{disptest}

\title{Dispersion Tests}

\description{
 Tests the null hypothesis of equidispersion in Poisson GLMs against
 the alternative of overdispersion and/or underdispersion.
}

\usage{
disptest(object, type = c("DL", "DLadj"), trafo = NULL, 
  alternative = c("greater", "two.sided", "less"))
}

\arguments{
 \item{object}{a fitted Poisson GLM of class \code{"glm"} as fitted
    by \code{\link{glm}} with family \code{\link{poisson}}.}
 \item{type}{type of test, one of \code{DL} or \code{DLadj}. See details.}  
 \item{trafo}{a specification of the alternative (see also details),
   can be numeric or a (positive) function or \code{NULL} (the default).}
 \item{alternative}{a character string specifying the alternative hypothesis:
   \code{"greater"} corresponds to overdispersion, \code{"less"} to
   underdispersion and \code{"two.sided"} to either one.}
}

\details{
The standard Poisson GLM models the (conditional) mean
\eqn{\mathsf{E}[y] = \mu}{E[y] = mu} which is assumed to be equal to the
variance \eqn{\mathsf{VAR}[y] = \mu}{VAR[y] = mu}. \code{disptest}
assesses the hypothesis that this assumption holds (equidispersion) against
the alternative that the variance is of the form:
  \deqn{\mathsf{VAR}[y] \quad = \quad \mu \; + \; \alpha \cdot \mathrm{trafo}(\mu).}{VAR[y] = mu + alpha * trafo(mu).}
Overdispersion corresponds to \eqn{\alpha > 0}{alpha > 0} and underdispersion to
\eqn{\alpha < 0}{alpha < 0}. The coefficient \eqn{\alpha}{alpha} can be estimated
by an auxiliary OLS regression and tested with the corresponding t (or z) statistic
which is asymptotically standard normal under the null hypothesis.

Common specifications of the transformation function \eqn{\mathrm{trafo}}{trafo} are
\eqn{\mathrm{trafo}(\mu) = \mu^2}{trafo(mu) = mu^2} or \eqn{\mathrm{trafo}(\mu) = \mu}{trafo(mu) = mu}.
The former corresponds to a negative binomial (NB) model with quadratic variance function
(called NB2 by Cameron and Trivedi, 2005), the latter to a NB model with linear variance
function (called NB1 by Cameron and Trivedi, 2005) or quasi-Poisson model with dispersion 
parameter, i.e.,
  \deqn{\mathsf{VAR}[y] \quad = \quad (1 + \alpha) \cdot \mu = \mathrm{dispersion} \cdot \mu.}{VAR[y] = (1 + alpha) * mu = dispersion * mu.}

By default, for \code{trafo = NULL}, the latter dispersion formulation is used in
\code{dispersiontest}. Otherwise, if \code{trafo} is specified, the test is formulated
in terms of the parameter \eqn{\alpha}{alpha}. The transformation \code{trafo} can either
be specified as a function or an integer corresponding to the function \code{function(x) x^trafo},
such that \code{trafo = 1} and \code{trafo = 2} yield the linear and quadratic formulations
respectively.

Type \code{"DL"} corresponds to the statistic \eqn{T_1} in Dean and Lawless (1989), type \code{"DLadj"} is their \eqn{T_a}.

}

\value{An object of class \code{"htest"}.}

\references{
Cameron AC, Trivedi PK (1990). \dQuote{Regression-based Tests for Overdispersion in the Poisson Model}. \emph{Journal of Econometrics}, \bold{46}, 347--364.

Cameron AC, Trivedi PK (2005). \emph{Microeconometrics: Methods and Applications}.
Cambridge: Cambridge University Press.

Cameron AC, Trivedi PK (2013). \emph{Regression Analysis of Count Data}, 2nd ed. 
Cambridge: Cambridge University Press.

Dean CB (1992). \dQuote{Testing for Overdispersion in Poisson and Binomial Regression Models}. \emph{Journal of the American Statistical Association}, \bold{87}, 451--457.

Dean C, Lawless JF (1989). \dQuote{Tests for Detecting Overdispersion in Poisson Regression Models}. \emph{Journal of the American Statistical Association}, \bold{84}, 467--472.

Jaggia S, Thosar S (1993). \dQuote{Multiple Bids as a Consequence of Target Management Resistance: A Count Data Approach}. \emph{Review of Quantitative Finance and Accounting}, \bold{3}, 447--457.

Lee LF (1986). \dQuote{Specification Test for Poisson Regression Models}. \emph{International Economic Review}, \bold{27}, 689--706.
}

\seealso{\code{\link{glm}}, \code{\link{poisson}}, \code{\link[MASS]{glm.nb}}}

\examples{
data("TakeoverBids")
m <- glm(bids ~ . + I(size^2), data = TakeoverBids, family = poisson)

## Jaggia and Thosar (1993), Table 3
disptest(m, type = "DL", alternative = "two.sided")
disptest(m, type = "DLadj", alternative = "two.sided")


}

\keyword{htest}