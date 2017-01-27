\name{pithist}
\alias{pithist}

\title{PIT Histograms for Assessing Goodness of Fit of Probability Models}

\description{
  PIT histograms graphically compare empirical probabilities from
  fitted models with a uniform distribution.
}

\usage{
pithist(object, type = c("random", "proportional"), nsim = 1L,
  breaks = NULL, xlim = c(0, 1), ylim = NULL,
  xlab = "PIT", ylab = "Density", main = NULL,
  border = "black", fill = "lightgray", col = "#B61A51",
  lwd = 2, lty = 1, freq = FALSE, ...)
}

\arguments{
  \item{object}{an object from which a \code{\link{pit}} can be
    extracted.}
  \item{type}{character. In case of discrete distributions should
    the PITs be drawn randomly from the corresponding interval
    or distributed proportionally?}
  \item{nsim}{integer. If \code{type} is \code{"random"} how
    many simulated PITs should be drawn?}
  \item{breaks}{numeric. Breaks for the histogram intervals.}
  \item{xlab, ylab, main}{graphical parameters.}
  \item{xlim, ylim, border, fill, col, lwd, lty}{graphical
    parameters. These may pertain either to the whole plot or just the
    histogram or just the fitted line.}
  \item{freq, \dots}{further arguments passed to
    \code{\link[graphics]{hist}}.}
}

\details{
PIT histograms graphically the probability integral transform (PIT),
i.e., observed probabilities from fitted probability models, with
a uniform distribution. It leverages the \code{\link{pit}} generic
and then essentially draws a \code{\link[graphics]{hist}}.

In case of discrete distributions the PIT is either drawn randomly
from the corresponding interval or distributed proportionally in the
histogram.
}

\seealso{\code{\link{pit}}, \code{\link[graphics]{hist}}}

\examples{
## count data regression models: crab satellites
data("CrabSatellites", package = "countreg")
cs_p   <-    glm(satellites ~     width + color, data = CrabSatellites, family = poisson)
cs_hnb <- hurdle(satellites ~ 1 | width + color, data = CrabSatellites, dist = "negbin")

## PIT histograms
par(mfrow = c(1, 2))
pithist(cs_p, ylim = c(0, 4.2))
pithist(cs_hnb, ylim = c(0, 4.2))
par(mfrow = c(1, 1))

}

\keyword{hplot}
