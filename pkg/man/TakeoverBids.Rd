\name{TakeoverBids}
\alias{TakeoverBids}
\title{Takeover Bids Data}

\description{
Firms that were targets of takeover bids during the period 1978--1985.  
}

\usage{data("TakeoverBids")}

\format{
A data frame containing 126 observations on 9 variables.
\describe{
  \item{bids}{Number of takeover bids (after the initial bid received by the target firm).}
  \item{legalrest}{factor. Equals \code{"yes"}  if target management responded by lawsuit.}
  \item{realrest}{factor. Equals \code{"yes"} if target management proposed changes in asset structure.}
  \item{finrest}{factor. Equals \code{"yes"} if target management proposed changes in ownership structure.}
  \item{whiteknight}{factor. Equals \code{"yes"} if target management invited friendly third-party bid.}
  \item{bidpremium}{Bid price divided by price 14 working days before bid.}
  \item{insthold}{Percentage of stock held by institutions.}
  \item{size}{Total book value of assets (in billions of USD).}
  \item{regulation}{factor. Equals \code{"yes"} if intervention by federal regulators.}
}
}

\details{
The data were originally used by Jaggia and Thosar (1993), where further details on the variables may be found.
}

\source{
Journal of Applied Econometrics Data Archive for Cameron and Johansson (1997). 

\url{http://qed.econ.queensu.ca/jae/1997-v12.3/cameron-johansson/}
}

\references{
Cameron AC, Johansson P (1997). \dQuote{Count Data Regression Using Series
Expansion: With Applications}, \emph{Journal of Applied Econometrics}, \bold{12}(3), 203--224.

Cameron AC, Trivedi PK (2013). \emph{Regression Analysis of Count Data}, 2nd ed. 
Cambridge: Cambridge University Press.

Jaggia S, Thosar S (1993). \dQuote{Multiple Bids as a Consequence of Target Management Resistance:
A Count Data Approach}, \emph{Review of Quantitative Finance and Accounting}, \bold{3}, 447--457.
}

\examples{
data("TakeoverBids", package = "countreg")

## Poisson model:
## Jaggia and Thosar (1993), Table 3
## Cameron and Johansson (1997), Table IV
tb_p <- glm(bids ~ . + I(size^2), data = TakeoverBids, family = poisson)
summary(tb_p)
logLik(tb_p)

## dispersion tests
## Cameron and Trivedi (2013, p. 185)
AER::dispersiontest(tb_p, alternative = "less", trafo = 2)
AER::dispersiontest(tb_p, alternative = "less", trafo = 1)

## visualization of underdispersion
rootogram(tb_p)
qqrplot(tb_p, range = c(0.05, 0.95))

## Parts of Cameron and Trivedi (2013), Table 5.4
summary(residuals(tb_p, type = "response"))
summary(residuals(tb_p, type = "pearson"))
summary(residuals(tb_p, type = "deviance"))

## hurdle Poisson model mitigates underdispersion
tb_hp <- hurdle(bids ~ . + I(size^2), data = TakeoverBids, dist = "poisson")
AIC(tb_p, tb_hp)
rootogram(tb_hp)
qqrplot(tb_hp, range = c(0.05, 0.95))

}

\keyword{datasets}
