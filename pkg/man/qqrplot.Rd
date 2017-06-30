\name{qqrplot}
\alias{qqrplot}

\title{Q-Q Plots for Quantile Residuals}

\description{
  Visualize goodness of fit of regression models by Q-Q plots using
  quantile residuals.
}

\usage{
qqrplot(object, type = c("random", "quantile"),
  nsim = 1L, prob = 0.5, range = FALSE, diag = TRUE,
  col = "black", fill = "lightgray", xlim = NULL, ylim = NULL,
  main = "Q-Q residuals plot", xlab = "Theoretical quantiles",
  ylab = "Quantile residuals", \dots)
}

\arguments{
  \item{object}{an object (for which a \code{\link{qresiduals}} method exists).}
  \item{type, nsim, prob}{arguments passed to \code{qresiduals}.}
  \item{range}{logical or quantile specification. Should the range of
    quantiles of the randomized quantile residuals be visualized? If \code{TRUE},
    then \code{range = c(0.01, 0.99)} is used.}
  \item{diag}{logical or color specification. Should a diagonal
    reference line be drawn?}
  \item{col, fill, xlim, ylim, main, xlab, ylab, \dots}{graphical plotting parameters
    passed to \code{\link[graphics]{plot}} or \code{\link[graphics]{points}},
    respectively.}
}

\details{
Q-Q residual plots essentially correspond to carrying out a
\code{qqnorm(qresiduals(object, \dots))}, i.e., a Q-Q plot comparing the quantile
residuals against theoretical quantiles from a standard normal distribution.
The function \code{qqrplot} is provided mainly to offer a few more options
for models with discrete responses where randomization or quantiles are
needed. See below for some examples.
}

\value{
  An list is returned invisibly with:
  \item{normal}{the theoretical normal quantiles,}
  \item{residuals}{the empirical quantile residuals.}
}

\references{
Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile Residuals.}
  \emph{Journal of Computational and Graphical Statistics}, \bold{5}, 1--10.
}

\seealso{\code{\link{qresiduals}}, \code{\link[stats]{qqnorm}}}

\examples{
## count data regression models: crab satellites
data("CrabSatellites", package = "countreg")
cs_p   <-    glm(satellites ~     width + color, data = CrabSatellites, family = poisson)
cs_nb  <- glm.nb(satellites ~     width + color, data = CrabSatellites)
cs_hp  <- hurdle(satellites ~ 1 | width + color, data = CrabSatellites, dist = "poisson")
cs_hnb <- hurdle(satellites ~ 1 | width + color, data = CrabSatellites, dist = "negbin")

## Q-Q residual plots
par(mfrow = c(2, 2))
qqrplot(cs_p, main = "Poisson")
qqrplot(cs_nb, main = "Negative Binomial")
qqrplot(cs_hp, main = "Hurdle Poisson")
qqrplot(cs_hnb, main = "Hurdle Negative Binomial")
par(mfrow = c(1, 1))

## Q-Q residual plots
par(mfrow = c(2, 2))
qqrplot(cs_p, main = "One Random Sample")
qqrplot(cs_p, main = "Median", type = "quantile")
qqrplot(cs_p, main = "10 Random Samples and Range", nsim = 10, range = c(0.005, 0.995))
qqrplot(cs_p, main = "100 Random Samples", nsim = 100, pch = 19, col = gray(0, alpha = 0.01))
par(mfrow = c(1, 1))

}

\keyword{hplot}
