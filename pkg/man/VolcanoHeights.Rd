\name{VolcanoHeights}
\alias{VolcanoHeights}

\title{Tukey's Volcano Heights}

\description{
  Heights of 218 volcanos taken from Tukey (1972).
}

\usage{data("VolcanoHeights")}

\format{
A numeric vector of 218 volcano heights (in 1000 feet).
}

\details{
The data are taken from Tukey (1972) who obtained them from
\emph{The World Almanac, 1966} (New York: The New York
World-Telegram and The Sun, 1966), pp. 282--283.
}

\source{
Figure 1 in Tukey (1972).
}

\references{
Tukey JW (1972). \dQuote{Some Graphic and Semigraphic Displays.}
  In Bancroft TA (ed.), \emph{Statistical Papers in Honor of George W. Snedecor},
  pp. 293--316. Iowa State University Press, Ames, IA.
  Reprinted in Cleveland WS (ed.): \emph{The Collected Works of John W. Tukey,
  Volume V. Graphics: 1965--1985}, Wadsworth & Brooks/Cole, Pacific Grove, CA, 1988.
}

\examples{
## Rootograms from Tukey (1972)
## (some 'breaks' don't match exactly)
data("VolcanoHeights", package = "countreg")

## Figure 16
rootogram(VolcanoHeights, "normal", style = "standing",
  breaks = 0:20 - 0.01, col = "transparent")

## Figure 17
rootogram(sqrt(1000 * VolcanoHeights), "normal", style = "standing",
  breaks = 0:17 * 10 - 1.1, col = "transparent")

## Figure 18
rootogram(sqrt(1000 * VolcanoHeights), "normal", style = "hanging",
  breaks = -2:18 * 10 - 1.1)

## Figure 19
rootogram(sqrt(1000 * VolcanoHeights), "normal", style = "suspended",
  breaks = -2:18 * 10 - 1.1, ylim = c(6, -2))
abline(h = c(-1.5, -1, 1, 1.5), lty = c(2, 3, 3, 2))
}

\keyword{datasets}
