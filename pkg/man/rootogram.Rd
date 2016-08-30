\name{rootogram}
\alias{rootogram}

\alias{plot.rootogram}
\alias{autoplot.rootogram}

\alias{+.rootogram}
\alias{c.rootogram}
\alias{rbind.rootogram}

\alias{rootogram.default}
\alias{rootogram.gam}
\alias{rootogram.gamlss}
\alias{rootogram.glm}
\alias{rootogram.hurdle}
\alias{rootogram.numeric}
\alias{rootogram.zeroinfl}
\alias{rootogram.zerotrunc}

\title{Rootograms for Assessing Goodness of Fit of Probability Models}

\description{
  Rootograms graphically compare (square roots) of empirical frequencies
  with fitted frequencies from a probability model.
}

\usage{
rootogram(object, \dots)

\method{rootogram}{default}(object, fitted, breaks = NULL,
  style = c("hanging", "standing", "suspended"),
  scale = c("sqrt", "raw"), plot = TRUE,
  width = NULL, xlab = NULL, ylab = NULL, main = NULL, \dots)

\method{plot}{rootogram}(x,
  xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, main = NULL,
  border = "black", fill = "lightgray", col = "#B61A51",
  lwd = 2, pch = 19, lty = 1, type = NULL, axes = TRUE, \dots)

\method{autoplot}{rootogram}(object, colour = c("black", "#B61A51"),
  fill = "darkgray", size = c(1.2, 4), \dots)
}

\arguments{
  \item{object}{an object. For the \code{default} method this needs to
    be a vector/table of observed frequencies. Currently, there are also
    methods for \code{numeric}, \code{glm}, \code{gam}, \code{gamlss},
    \code{hurdle}, and \code{zeroinfl} objects.}
  \item{fitted}{an specification of the fitted frequencies. For the
    \code{default} method this needs to be a vector/table. For the
    \code{numeric} method (see below) this can also be a character
    string. For the regression object methods, this is not needed.}
  \item{breaks}{numeric. Breaks for the histogram intervals.}
  \item{style}{character specifying the syle of rootogram (see below).}
  \item{scale}{character specifying whether raw frequencies or their
    square roots (default) should be drawn.}
  \item{plot}{logical. Should the \code{plot} method be called to
    draw the computed rootogram?}
  \item{width}{numeric. Widths of the histogram bars.}
  \item{xlab, ylab, main}{graphical parameters.}
  \item{x}{object of class \code{"rootogram"}.}
  \item{xlim, ylim, border, fill, col, lwd, pch, lty, type, axes}{graphical
    parameters. These may pertain either to the whole plot or just the
    histogram or just the fitted line.}
  \item{colour, size}{graphical parameters passed to \code{geom_line}
    and \code{geom_point}, respectively.}
  \item{\dots}{further graphical parameters passed to the plotting
    function.}
}

\details{
Rootograms compare graphically compare frequencies of empirical
distributions and fitted probability models. For the observed distribution
the histogram is drawn on a square root scale (hence the name) and superimposed
with a line for the fitted frequencies. The histogram can be \code{"standing"}
on the x-axis (as usual), or \code{"hanging"} from the fitted curve, or
a \code{"suspended"} histogram of deviations can be drawn.

\code{rootogram} is the generic function for generating rootograms from data
or fitted model objects. The workhorse function is the default method (that
computes all necessary coordinates based on observed and fitted frequencies and
the breaks for the histogram intervals) and the associating \code{plot} method
that carries out the actual drawing (using base graphics).

There is a wide range of further \code{rootogram} methods that all take
the following approach: based on a fitted probability model observed and
expected frequencies are computed and then the default method is called.
Currently, there are methods for \code{\link[stats]{glm}}, 
\code{\link[mgcv]{gam}}, \code{\link[gamlss]{gamlss}}, \code{\link{hurdle}}, 
\code{\link{zeroinfl}}, and \code{\link{zerotrunc}} models.

Furthermore, there is a \code{numeric} method that uses \code{link[MASS]{fitdistr}}
to obtain a fitted (by maximum likelihood) probability model for a univariate variable. 
For this method, \code{fitted} can either be a character string or a density function
that is passed to \code{fitdistr}. In the latter case, a \code{start} list also
has to be supplied.

In addition to the \code{plot} method for rootogram objects, there are also two
methods that combine two (or more) rootograms: \code{c}/\code{rbind} creates
a set of rootograms that can then be plotted in one go. The \code{+} method adds
up the observed and fitted frequencies from two rootograms (if these use the same
intervals).

The \code{\link[ggplot2]{autoplot}} method creates a \code{\link[ggplot2]{ggplot}}
version of the rootogram.
}

\note{
Note that there is also a \code{\link[vcd]{rootogram}} function in the \pkg{vcd}
package that is similar to the \code{numeric} method provided here. However, it
is much more limited in scope, hence a function has been created here.
}

\value{
  An object of class \code{"rootogram"} inheriting from \code{"data.frame"} with the
  following variables:
  \item{observed}{observed frequencies,}
  \item{expected}{fitted frequencies,}
  \item{x}{histogram interval midpoints on the x-axis,}
  \item{y}{bottom coordinate of the histogram bars,}
  \item{width}{widths of the histogram bars,}
  \item{height}{height of the histogram bars,}
  \item{line}{y-coordinates of the fitted curve.}
  Additionally, \code{style}, \code{scale}, \code{xlab}, \code{ylab}, and
  \code{main} are stored as attributes.
}

\references{
Friendly M (2000), \emph{Visualizing Categorical Data}. SAS Institute, Cary.

Kleiber C, Zeileis A (2016).
  \dQuote{Visualizing Count Data Regressions Using Rootograms.} 
  \emph{The American Statistician}, \bold{70}(3), 296--303.
  \doi{10.1080/00031305.2016.1173590}.

Tukey JW (1977). \emph{Exploratory Data Analysis}. Addison-Wesley, Reading.
}

\seealso{\code{\link[stats]{glm}}, \code{\link[MASS]{glm.nb}},
  \code{\link{hurdle}}, \code{\link{zeroinfl}}
}

\examples{
## different interfaces

## number of deaths by horsekicks in Prussian army (Von Bortkiewicz 1898)
deaths <- rep(0:4, c(109, 65, 22, 3, 1))

## default method: fitted values "by hand"
rootogram(table(deaths), fitted = length(deaths) * dpois(0:4, mean(deaths)))

## numeric method: fitted values via fitdistr()
rootogram(deaths, fitted = "poisson")
rootogram(deaths, fitted = dpois, start = list(lambda = 1),
  breaks = 0:5 - 0.5, width = 0.9)

## glm method: fitted values via glm()
m <- glm(deaths ~ 1, family = poisson)
rootogram(m)

## inspect output (without plotting)
r <- rootogram(m, plot = FALSE)
r

## create ggplot2 version
if(require("ggplot2")) autoplot(r)

#-------------------------------------------------------------------------------

## different styles

## artificial data from negative binomial (mu = 3, theta = 2)
## and Poisson (mu = 3) distribution
set.seed(1090)
y <- rnbinom(100, mu = 3, size = 2)
x <- rpois(100, lambda = 3)

## correctly specified Poisson model fit (mu = 3.34)
par(mfrow = c(2, 3))
rootogram(x, "poisson", style = "standing",  ylim = c(-2.2, 4.8), main = "Standing")
rootogram(x, "poisson", style = "hanging",   ylim = c(-2.2, 4.8), main = "Hanging")
rootogram(x, "poisson", style = "suspended", ylim = c(-2.2, 4.8), main = "Suspended")

## misspecified Poisson model fit (mu = 3.32)
rootogram(y, "poisson", style = "standing",  ylim = c(-2.2, 4.8), main = "Standing")
rootogram(y, "poisson", style = "hanging",   ylim = c(-2.2, 4.8), main = "Hanging")
rootogram(y, "poisson", style = "suspended", ylim = c(-2.2, 4.8), main = "Suspended")
par(mfrow = c(1, 1))

#-------------------------------------------------------------------------------

## artificial data from a t_4 distribution
set.seed(1090)
y <- rt(1000, 4)

## incorrect normal fit (tails too light) and correct t fit
par(mfrow = c(1, 2))
rootogram(y, fitted = "normal", breaks = 40, xlim = c(-6, 6), ylim = c(-2, 14))
rootogram(y, fitted = "t",      breaks = 40, xlim = c(-6, 6), ylim = c(-2, 14))
par(mfrow = c(1, 1))

#-------------------------------------------------------------------------------

## linear regression with normal/Gaussian response: anorexia data
an <- glm(Postwt ~ Prewt + Treat + offset(Prewt), family = gaussian, data = anorexia)
rootogram(an, ylim = c(-1, 4))
abline(h = c(-1, 1), col = "#1E55CE", lty = 2, lwd = 2)

#-------------------------------------------------------------------------------

## count data regression models: crab satellites
data("CrabSatellites", package = "countreg")
cs_p   <-    glm(satellites ~     width + color, data = CrabSatellites, family = poisson)
cs_nb  <- glm.nb(satellites ~     width + color, data = CrabSatellites)
cs_hp  <- hurdle(satellites ~ 1 | width + color, data = CrabSatellites, dist = "poisson")
cs_hnb <- hurdle(satellites ~ 1 | width + color, data = CrabSatellites, dist = "negbin")

## rootograms
par(mfrow = c(2, 2))
rootogram(cs_p, max = 15, main = "Poisson")
rootogram(cs_nb, max = 15, main = "Negative Binomial")
rootogram(cs_hp, max = 15, main = "Hurdle Poisson")
rootogram(cs_hnb, max = 15, main = "Hurdle Negative Binomial")
par(mfrow = c(1, 1))

}

\keyword{hplot}
