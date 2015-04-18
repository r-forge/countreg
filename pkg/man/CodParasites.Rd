\name{CodParasites}
\alias{CodParasites}

\encoding{latin1}

\title{Parasite Infections in Cod}

\description{
Data on parasite infection in cod along the coast of Finmark.
}

\usage{data("CodParasites")}

\format{
A data frame containing 1254 observations on 10 variables. 
\describe{
  \item{intensity}{Number of parasites.}
  \item{prevalence}{Factor indicating presence of parasites (i.e.,
    \code{intensity > 0}).}
  \item{area}{Factor indicating sampling area.}
  \item{year}{Factor indicating sampling year.}
  \item{depth}{Depth at which the fish were caught.}
  \item{weight}{Weight of the fish.}
  \item{length}{Length of the fish.}
  \item{sex}{Factor indicating sex of the fish.}
  \item{stage}{Factor indicating stage of the fish.}
  \item{age}{Age of the fish.}
}
}

\details{
The red king crab \emph{Paralithodes camtschaticus} was deliberately introduced to the
Barents Sea in the 1960s and 1970s from its native area in the North Pacific. The carapace
of these crabs is used by the leech \emph{Johanssonia arctica} to deposit its eggs.
The leech in turn is a vector for the blood parasite \emph{Trypanosoma murmanensis} that
can infect marine fish, including cod.

Hemmingsen et al. (2005) examined cod for trypanosome infections during annual cruises
along the coast of Finnmark in North Norway over three successive years and in four
different areas (A1 Sørøya; A2 Magerøya; A3 Tanafjord; A4 Varangerfjord). They show that
trypanosome infections are strongest in the area Varangerfjord where the density of
of red king crabs is highest. Thus, there is evidence that the introduction of the foreign
red king crabs had an indirect detrimental effect on the health of the native cod population.
This situation stands out because it is not an introduced parasite that is dangerous for
a native host, but rather an introduced host that promotes transmission of two endemic
parasites.

Zuur et al. (2009) reanalyze the data using binary and count data regression models
in Chapters 10.2.2, 11.3.2, 11.4.2, 11.5.2.
}

\source{
The data are taken from the online supplements of Zuur et al. (2009).
\url{http://www.highstat.com/book2.htm}
}

\references{
Hemmingsen W, Jansen PA, MacKenzie K (2005).
  \dQuote{Crabs, Leeches and Trypanosomes: An Unholy Trinity?},
  \emph{Marine Pollution Bulletin} \bold{50}(3), 336--339.

Zuur AF, Ieno EN, Walker NJ, Saveliev AA, Smith GM (2009).
  \emph{Mixed Effects Models and Extensions in Ecology with R},
  Springer-Verlag, New York.
}

\examples{
## load data
data("CodParasites", package = "countreg")

## Table 1 from Hemmingsen et al. (2005)
## number of observations
xtabs(~ area + year, data = CodParasites)
## prevalence of parasites (NAs counted as "yes")
tab <- xtabs(~ area + year + factor(is.na(prevalence) | prevalence == "yes"),
  data = CodParasites)
round(100 * prop.table(tab, 1:2)[,,2], digits = 1)

## omit NAs in response
CodParasites <- subset(CodParasites, !is.na(intensity))

## exploratory displays for hurdle and counts
par(mfrow = c(2, 2))
plot(factor(intensity == 0) ~ interaction(year, area), data = CodParasites)
plot(factor(intensity == 0) ~ length, data = CodParasites, breaks = c(15, 3:8 * 10, 105))
plot(jitter(intensity) ~ interaction(year, area), data = CodParasites,
  subset = intensity > 0, log = "y")
plot(jitter(intensity) ~ length, data = CodParasites, subset = intensity > 0, log = "y")

## count data models
cp_p   <-    glm(intensity ~ length + area * year, data = CodParasites, family = poisson)
cp_nb  <- glm.nb(intensity ~ length + area * year, data = CodParasites)
cp_hp  <- hurdle(intensity ~ length + area * year, data = CodParasites, dist = "poisson")
cp_hnb <- hurdle(intensity ~ length + area * year, data = CodParasites, dist = "negbin")
AIC(cp_p, cp_nb, cp_hp, cp_hnb)
BIC(cp_p, cp_nb, cp_hp, cp_hnb)

## rootograms
par(mfrow = c(2, 2))
rootogram(cp_p, max = 50, main = "Poisson")
rootogram(cp_nb, max = 50, main = "Negative Binomial")
rootogram(cp_hp, max = 50, main = "Hurdle Poisson")
rootogram(cp_hnb, max = 50, main = "Hurdle Negative Binomial")
}

\keyword{datasets}
