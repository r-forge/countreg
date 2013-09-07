\name{CrabSatellites}
\alias{CrabSatellites}

\title{Horseshoe Crab Mating}

\description{
  Determinants for male satellites to nesting horseshoe crabs.
}

\usage{data("CrabSatellites")}

\format{
A data frame containing 173 observations on 5 variables.
\describe{
  \item{color}{Ordered factor indicating color
    (light medium, medium, dark medium, dark).}
  \item{spine}{Ordered factor indicating spine condition
    (both good, one worn or broken, both worn or broken).}
  \item{width}{Carapace width (cm).}
  \item{weight}{Weight (kg).}
  \item{satellites}{Number of satellites.}
}
}

\details{
Brockmann (1996) investigates horshoe crab mating. The crabs arrive on the beach
in pairs to spawn. Furthermore, unattached males also come to the beach, crowd
around the nesting couples and compete with attached males for fertilizations.
These so-called satellite males form large groups around some couples while
ignoring others. Brockmann (1996) shows that the groupings are not driven
by environmental factors but by properties of the nesting female crabs.
Larger females that are in better condition attract more satellites.

Agresti (2002) reanalyzes the number of satellites using count models.
Explanatory variables are the female crab's color, spine condition, weight,
and carapace width. Color and spine condition are ordered factors but are
treated as numeric in some analyses.
}

\source{
Table 4.3 in Agresti (2002).
}

\references{
Agresti A (2002).
  \emph{Categorical Data Analysis}, 2nd ed.,
  John Wiley & Sons, Hoboken.

Brockmann HJ (1996).
  \dQuote{Satellite Male Groups in Horseshoe Crabs, \emph{Limulus polyphemus}},
  \emph{Ethology}, \bold{102}(1), 1--21.
}

\examples{
## load data and use ordered factors as numeric
data("CrabSatellites", package = "countreg")
CrabSatellites <- transform(CrabSatellites,
  color = as.numeric(color),
  spine = as.numeric(spine)
)  

## exploratory displays for hurdle and counts
par(mfrow = c(2, 2))
plot(factor(satellites == 0) ~ width, data = CrabSatellites, breaks = seq(20, 33.5, by = 1.5))
plot(factor(satellites == 0) ~ color, data = CrabSatellites, breaks = 1:5 - 0.5)
plot(jitter(satellites) ~ width, data = CrabSatellites, subset = satellites > 0, log = "y")
plot(jitter(satellites) ~ factor(color), data = CrabSatellites, subset = satellites > 0, log = "y")

## count data models
cs_p   <-    glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
cs_nb  <- glm.nb(satellites ~ width + color, data = CrabSatellites)
cs_hp  <- hurdle(satellites ~ width + color, data = CrabSatellites, dist = "poisson")
cs_hnb <- hurdle(satellites ~ width + color, data = CrabSatellites, dist = "negbin")
AIC(cs_p, cs_nb, cs_hp, cs_hnb)
BIC(cs_p, cs_nb, cs_hp, cs_hnb)

## rootograms
par(mfrow = c(2, 2))
rootogram(cs_p, max = 15, main = "Poisson")
rootogram(cs_nb, max = 15, main = "Negative Binomial")
rootogram(cs_hp, max = 15, main = "Hurdle Poisson")
rootogram(cs_hnb, max = 15, main = "Hurdle Negative Binomial")

## fitted curves
par(mfrow = c(1, 1))
plot(jitter(satellites) ~ width, data = CrabSatellites)
nd <- data.frame(width = 20:34, color = 2.5)
cs_ag <- glm(satellites ~ width, data = CrabSatellites, family = poisson(link = "identity"),
  start = coef(lm(satellites ~ width, data = CrabSatellites)))
with(CrabSatellites, lines(lowess(width, satellites)))
lines(predict(cs_ag, nd) ~ width, data = nd, col = 4, lwd = 1.5)
lines(predict(cs_p,  nd) ~ width, data = nd, col = 4, lty = 2, lwd = 1.5)
lines(predict(cs_hnb, nd) ~ width, data = nd, col = 2, lwd = 1.5)
legend("topleft", c("Lowess", "Hurdle-NB", "Poisson (id)", "Poisson (log)"),
  col = c(1, 2, 4, 4), lty = c(1, 1, 1, 2), lwd = 1.5, bty = "n")
}

\keyword{datasets}