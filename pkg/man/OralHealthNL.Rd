\name{OralHealthNL}
\alias{OralHealthNL}

\title{Oral Health in Children in The Netherlands}

\description{
  Data from a study on oral health status and the
  preventive dental behaviors of 9-year-old children
  in The Netherlands.
}

\usage{data("OralHealthNL")}

\format{
A data frame containing 440 observations on 8 variables.
\describe{
  \item{dmfs}{Numeric index of decayed, missing, and filled surfaces
    (DMFS) in deciduous teeth.}
  \item{education}{Factor indicating whether the highest completed
    education level of the mother is \code{"high"} (senior general
    secondary education, HAVO, or higher) or \code{"low"}.}
  \item{gender}{Factor indicating gender of the child
    (\code{"female"} or \code{"male"}).}
  \item{ethnicity}{Factor indicating whether the mother
    is \code{"immigrant"} (born abroad) or \code{"native"}
    (born in The Netherlands).}
  \item{brushing}{Factor indicating whether the frequency of
    brushing teeth is \code{"< 2"} or \code{">= 2"} times
    per day.}
  \item{breakfast}{Factor indicating whether the frequency of
    having breakfast is \code{"7"} or \code{"< 7"} days per week.}
  \item{fooddrink}{Factor indicating whether the frequency of
    food and drinks in addition to the three main meals is
    \code{"<= 7"} or \code{"> 7"} times per day.}
  \item{corah}{Factor indicating whether Corah's Dental
    Anxiety score is \code{"< 13"} or \code{">= 13"} (see
    also below).}
}
}

\details{
The data are from the study \dQuote{Oral Health in Children and
Adolescents in The Netherlands} (Schuller et al. 2011). The aim of
this study was to describe the oral health status and the
preventive dental behaviors of children from different age
groups (Dusseldorp et al. 2015). Here, the subset of children at
the age of 9 years is provided as analyzed by Hofstetter et al. (2016). 

The data collection consisted of a clinical oral examination and a
questionnaire survey, using a repeated cross-sectional design.
Data contained information about demographic variables (ethnicity
and educational level), nutrition, children's dental attendance,
oral self-care, and dental anxiety. The score on Corah's Dental Anxiety
Questionnaire was used as a measure of dental anxiety.
This questionnaire consists of four questions with answer
categories from 1 (low anxiety) to 5 (high anxiety). A total
Corah score was computed by taking the sum of the four items and
then dichotomized into \sQuote{lower than 13} and \sQuote{higher than
or equal to 13}.
}

\source{
Supplementary materials for Hofstetter et al. (2016).
\doi{10.1159/000448197}
}

\references{
Dusseldorp E, Kamphuis M, Schuller AA (2015).
  \dQuote{Impact of Lifestyle Factors on Caries Experience in Three Different Age Groups: 9, 15, and 21-Year-Olds},
  \emph{Community Dentistry and Oral Epidemiology}, \bold{43}(1), 9--16.
  \doi{10.1111/cdoe.12123}

Hofstetter H, Dusseldorp E, Zeileis A, Schuller AA (2016).
  \dQuote{Modeling Caries Experience: Advantages of the Use of the Hurdle Model},
  \emph{Caries Research}, \bold{50}(6), 517--526.
  \doi{10.1159/000448197}

Schuller AA, Poorterman JHG, van Kempen CPF, Dusseldorp E, van Dommelen P, Verrips GHW (2011).
  \emph{Kies voor tanden: Een onderzoek naar mondgezondheid en preventief tandheelkundig gedrag van jeugdigen.
  Tussenmeting 2009, een vervolg op de reeks TJZ-onderzoeken}. TNO, Leiden.
}

\examples{
## Load data and omit NAs and one dmfs outlier
data("OralHealthNL", package = "countreg")
head(OralHealthNL)
OralHealthNL <- na.omit(subset(OralHealthNL, dmfs < 40))

## Visualization: Is dmfs > 0?
par(mfrow = c(2, 4))
plot(factor(dmfs > 0, levels = c(TRUE, FALSE), labels = c("> 0", "= 0")) ~ .,
  data = OralHealthNL, ylevels = 2:1, ylab = "dmfs")

## Count: How large is log(dmfs) given dmfs > 0?
par(mfrow = c(2, 4))
plot(log(dmfs) ~ ., data = OralHealthNL, subset = dmfs > 0, ylab = "dmfs")

## Relevel the factor variables so that non-risk group is the reference
OralHealthNL <- transform(OralHealthNL,
  ethnicity = relevel(ethnicity, ref = "native"),
  brushing = relevel(brushing, ref = ">= 2"),
  breakfast = relevel(breakfast, ref = "7")
)

## Count regression models
zinb <- zeroinfl(dmfs ~ ., data = OralHealthNL, dist = "negbin")
zip  <- zeroinfl(dmfs ~ ., data = OralHealthNL, dist = "poisson")
hnb  <-   hurdle(dmfs ~ ., data = OralHealthNL, dist = "negbin")
hp   <-   hurdle(dmfs ~ ., data = OralHealthNL, dist = "poisson")

## Model comparisons (Table 3)
## Information criteria
cbind(AIC(hnb, zinb, hp, zip), BIC = BIC(hnb, zinb, hp, zip)[, 2])
## Negative binomial vs. Poisson
if(require("lmtest")) lrtest(hnb, hp)
if(require("lmtest")) lrtest(zinb, zip)
## Zero-inflation vs. hurdle
if(require("nonnest2")) vuongtest(zinb, hnb)

## Coefficients, odds ratios, and rate ratios
## Negative binomial hurdle model (Table 3)
summary(hnb)
exp(confint(hnb))
## Negative binomial zero-inflated model (Table 4)
summary(hnb)
exp(confint(zinb))

## Rootograms (top left: Figure 1)
par(mfrow = c(2, 2))
rootogram(OralHealthNL$dmfs, "normal",
  style = "standing", scale = "raw",
  breaks = 0:23 - 0.5, xlim = c(-0.5, 22.5),
  xlab = "dmfs", main = "Normal distribution")
rootogram(hnb,
  style = "standing", scale = "raw",
  width = 1, xlim = c(-0.5, 22.5),
  xlab = "dmfs", main = "Negative binomial hurdle model")
rootogram(OralHealthNL$dmfs, "normal",
  breaks = 0:23 - 0.5, xlim = c(-0.5, 22.5),
  xlab = "dmfs", main = "Normal distribution")
abline(h = c(-1, 1), lty = 2)
rootogram(hnb,
  width = 1, xlim = c(-0.5, 22.5),
  xlab = "dmfs", main = "Negative binomial hurdle model")
abline(h = c(-1, 1), lty = 2)
par(mfrow = c(1, 1))

## Number of zeros
c(dmfs = sum(OralHealthNL$dmfs == 0),
  ZINB = sum(predict(zinb, type = "prob")[,1]),
  Hurdle = sum(predict(hnb, type = "prob")[,1]))
## Correlation of observations and fitted means
cor(cbind(dmfs = OralHealthNL$dmfs,
  ZINB = fitted(zinb), HNB = fitted(hnb)))

## Bias-reduced logistic regression (due to separation)
if(require("brglm")) {
br <- brglm::brglm(
  factor(dmfs == 0, levels = c(TRUE, FALSE), labels = c("= 0", "> 0")) ~ .,
  data = OralHealthNL)
print(coeftest(br), digits = 1)
}

}

\keyword{datasets}
