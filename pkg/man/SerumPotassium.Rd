\name{SerumPotassium}
\alias{SerumPotassium}

\title{Serum Potassium Levels}

\description{
  Sample of 152 serum potassium levels.
}

\usage{data("SerumPotassium")}

\format{
A numeric vector of 152 serum potassium levels.
}

\details{
The data are taken from Rice (2007) who obtained
the data from Martin, Gudzinowicz and Fanger (1975) 
and reports them rounded to one digit.
}

\source{
Page 350 in Rice (2007).
}

\references{
Rice JA (2007). \emph{Mathematical Statistics and Data Analysis},
  3rd ed. Duxbury, Belmont, CA.

Martin HF, Gudzinowicz BJ, Fanger H (1975).
  \emph{Normal Values in Clinical Chemistry: A Guide to Statistical Analysis of Laboratory Data}.
  Marcel Dekker, New York.

}

\examples{
\dontshow{ if(!requireNamespace("topmodels")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("the topmodels package is required for this example but is not installed")
  } else q() }
}
library("topmodels")
data("SerumPotassium", package = "countreg")

## Figure 9.3a-c from Rice (2007), and actual hanging rootogram
## (note that Rice erroneously refers to suspended rootograms as hanging)
br <- 32:54/10 - 0.05
rootogram(lm(SerumPotassium ~ 1), scale = "raw", style = "standing",
  breaks = br, col = "transparent")
rootogram(lm(SerumPotassium ~ 1), scale = "raw", style = "suspended",
  breaks = br, col = "transparent", ylim = c(2.8, -4))
rootogram(lm(SerumPotassium ~ 1), scale = "sqrt", style = "suspended",
  breaks = br, col = "transparent", ylim = c(1, -1.5))
rootogram(lm(SerumPotassium ~ 1), breaks = br)
}

\keyword{datasets}
