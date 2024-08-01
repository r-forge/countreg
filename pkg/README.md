

<!-- README.md is generated from README.qmd via: quarto render README.qmd --to gfm -->

<img src="https://topmodels.R-Forge.R-project.org/countreg/countreg.png" align="right" alt="countreg logo" width="150" />

# countreg: Count Data Regression in R

## Overview

The core of the `countreg` package consists of functions for a variety
of count data regression models:

-   `nbreg()`: *Negative binomial regression* (type 1 and 2), including
    both constant and varying dispersion parameters which can depend on
    regressors.
-   `zeroinfl()` and `hurdle()`: *Zero-inflated and hurdle models* for
    count data based on Poisson, geometric, negative binomial, and
    binomial distributions. These allow for excess zeros compared to the
    unadjusted underlying distributions (or also a deficit of zeros in
    case of the hurdle model).
-   `zerotrunc()`: *Zero-truncated count regression* based on Poisson,
    geometric, and negative binomial distributions. These models can be
    employed for the zero-truncated count part of hurdle models.
-   `FLXMRnegbin()`: *Finite mixtures of negative binomial regressions*
    via driver for the
    [flexmix](https://CRAN.R-project.org/package=flexmix) package.
-   `MBnegbin()`: *Model-based boosting of negative binomial
    regressions* via driver for the
    [mboost](https://CRAN.R-project.org/package=mboost) package.

Previously available functions for *graphical goodness-of-fit
assessment* (rootograms, PIT histograms, Q-Q plots based on randomized
quantile residulas, etc.) are now provided in the
[topmodels](https://topmodels.R-Forge.R-project.org) package.

*Distribution functions* (d/p/q/r) for the different zero-truncated,
zero-inflated, and hurdle count distributions are availabe in the
[distributions3](https://alexpghayes.github.io/distributions3/) package
along with corresponding distribution objects.

## Installation

The latest development version can be installed from R-universe:

``` r
install.packages("countreg", repos = "https://zeileis.R-universe.dev")
```

## License

The package is available under the [General Public License version
3](https://www.gnu.org/licenses/gpl-3.0.html) or [version
2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)

## Illustration

For demonstrating some of the available count regression models we
analyze the number of male satellites around nesting female horseshoe
crabs during mating. Because larger females that are in better condition
tend to attract more satellites we consider their carapace width and
color as (numeric) regressor variables. Package and data can be loaded
and tranformed via:

``` r
library("countreg")
## Loading required package: MASS
data("CrabSatellites", package = "countreg")
cs <- CrabSatellites |>
  transform(color = as.numeric(color)) |>
  subset(select = c("satellites", "width", "color"))
```

The following four count models are compared: a Poisson generalized
linear model via `glm()`, a negative binomial regression with `nbreg()`,
a negative binomial hurdle model with `hurdle()`, and a zero-inflated
negative binomial model with `zeroinfl()`.

``` r
cs_p    <-      glm(satellites ~ ., data = cs, family = poisson)
cs_nb   <-    nbreg(satellites ~ ., data = cs)
cs_hnb  <-   hurdle(satellites ~ ., data = cs, dist = "negbin")
cs_zinb <- zeroinfl(satellites ~ ., data = cs, dist = "negbin")
```

The `glm()` function is from base R, all other models are from
`countreg`. All negative binomial models employ the type 2
parameterization. The hurdle and zero-inflated models employ binary
models for the zero hurdle and zero-inflated part, respectively.

For a quick overview we compare the estimated coefficients from all
parts of the models and compute their BIC.

``` r
cbind(
  poisson = coef(cs_p),
  nb = coef(cs_nb, model = "mu"),
  hnb_count = coef(cs_hnb, model = "count"),
  hnb_zero = coef(cs_hnb, model = "zero"),
  zinb_count = coef(cs_zinb, model = "count"),
  zinb_zero = coef(cs_zinb, model = "zero")
)
##                poisson         nb   hnb_count    hnb_zero zinb_count  zinb_zero
## (Intercept) -2.5199832 -3.2409843 0.428566899 -10.0708390 0.29870287 10.4341839
## width        0.1495725  0.1776543 0.037845165   0.4583097 0.04171687 -0.4890186
## color       -0.1694036 -0.1815656 0.006928678  -0.5090467 0.01866320  0.6068262
BIC(cs_p, cs_nb, cs_hnb, cs_zinb)
##         df      BIC
## cs_p     3 930.9589
## cs_nb    4 769.5455
## cs_hnb   7 736.7985
## cs_zinb  7 736.8958
```

This shows that

-   The number of male satellites increase with `width` and decrease
    with `color` of the female carapace.
-   The effect of the regressors mostly distinguishes between having any
    satellites or not while the expectation in the count part has almost
    no dependency on the regressors.
-   There is overdispersion in the number of satellites because the
    negative binomial models fit much better than the simple Poisson
    model.
-   There is only litte difference between the hurdle and zero-inflated
    model (except that the sign in the zero part is flipped because one
    models the probability to exceed zero while the other models the
    probability of an excess zero).

The summary of the hurdle negative binomial model provides further
details such as the non-significance of the regressors in the count part
and the estimate of the overdispersion parameter `theta`.

``` r
summary(cs_hnb)
## 
## Call:
## hurdle(formula = satellites ~ ., data = cs, dist = "negbin")
## 
## Pearson residuals:
##     Min      1Q  Median      3Q     Max 
## -1.3835 -0.7244 -0.2636  0.5557  3.6080 
## 
## Count model coefficients (truncated negbin with log link):
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) 0.428567   0.941077   0.455    0.649    
## width       0.037845   0.032749   1.156    0.248    
## color       0.006929   0.091078   0.076    0.939    
## Log(theta)  1.527382   0.352950   4.327 1.51e-05 ***
## Zero hurdle model coefficients (binomial with logit link):
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -10.0708     2.8065  -3.588 0.000333 ***
## width         0.4583     0.1040   4.407 1.05e-05 ***
## color        -0.5090     0.2237  -2.276 0.022862 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Theta: count = 4.6061
## Number of iterations in BFGS optimization: 17 
## Log-likelihood: -350.4 on 7 Df
```
