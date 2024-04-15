# countreg 0.3-0

* Major revision of the package, refactoring some crucial infrastructure.
  Rootograms and other visualizations are in `topmodels` (currently on R-Forge),
  see <https://topmodels.R-Forge.R-project.org/articles/topmodels.html>.
  Distribution functions (`d`/`p`/`q`/`r`) are in `distributions3` (on CRAN),
  see <https://www.zeileis.org/news/user2022/>.

* Also, data sets `SerumPotassium` and `VolcanoHeights` are now in `topmodels`.

* Improved `predict()` method with more `type` of predictions, including
  moments (`type = "mean"` or `"variance"`), `"quantile"`, `"density"`,
  cumulative distribution function (`"probability"`) etc.
  By default, the prediction is computed for the `"full"` outcome model
  (e.g., hurdle or zero-inflation) but can also be just the `"count"` or
  `"zero"` component or the distribution `"truncated"` at zero.


# countreg 0.2-1

* Bug fix in `dhpois()` and `dhnbinom()` when `x` is a vector and `pi` is a
  scalar (reported by Andrea Gilardi).

* Try to auto-detect columns in regressor matrices that are aliased, i.e.,
  whose coefficients cannot be estimated. This is based on pivot and rank
  from `qr()` and can detect, e.g., constant or linearly dependent columns.

* Catch errors when inverting the Hessian from optim, returning an `NA`
  matrix of correct dimension instead.

* Conditionally register `autoplot.rootogram()` method.


# countreg 0.2-0

* New `d`/`p`/`q`/`r`/`s` functions for all combinations of Poisson vs. negative
  binomial (NB) and zero-truncated vs. zero-inflated vs. hurdle. The parameter
  names are all greek letters: `lambda` for Poisson mean, `mu` for NB mean,
  `theta` for NB dispersion, `pi` for zero-modification probability (either
  inflation or hurdle crossing). The `theta` parameter may also be called
  `size` for consistency with base R's `nbinom` family.

* New generic function `pit()` for extracting the probability integral
  transform F(y) where F is the predicted cumulative density function
  and y the observed response.
  
* Quantile residuals are now computed based on the `pit()` extractor.

* New graphics function `pithist()` for PIT histograms based on `pit()`.

* New data set `OralHealthNL` accompanying Hofstetter et al. (2016):
  "Modeling Caries Experience: Advantages of the Use of the Hurdle Model",
  _Caries Research_, 50(6), 517-526.
  [doi:10.1159/000448197](https://doi.org/10.1159/000448197)


# countreg 0.1-5

* New `d`/`p`/`q`/`r`/`s` functions for the zero-truncated Poisson distribution
  where the parameter can either be specified in terms of the zero-truncated
  `mean` or the untruncated mean `lambda`.
  
* New `ztpoisson()` family object for the estimation of zero-truncated
  Poisson regression model via `glm()`.


# countreg 0.1-4

* New `qresiduals()` generic for computing (randomized) quantile residuals
  along with methods for various objects. (This is somewhat more flexible
  than `statmod::qresiduals()` which is not generic.)
  
* A Q-Q plot based on quantile residuals is available in the new function
  `qqrplot()`.


# countreg 0.1-3

* The `hurdle()` function now allows for the restriction of the theta
  parameter for the negative binomial distribution across the censored
  zero and truncated count components. The restriction is applied if
  `dist = "negbin"`, `zero.dist = "negbin"`, and `separate = FALSE`.

* Bug fix in the computation of the standard error of log(theta)
  in `zeroinfl()` when with weighted data. The weights were used
  for the estimation of log(theta) but the standard error was not
  scaled accordingly.


# countreg 0.1-2

* The `style` and `scale` are now added as attributes to `rootogram`
  objects so that they can be re-used in the `c()` and `+` methods.


# countreg 0.1-1

* New data sets `SerumPotassium`, `TakeoverBids`, `VolcanoHeights`.

* New function `disptest()` providing several (score) tests for over/underdispersion.

* New `rootogram()` function with a wide range methods for creating rootograms
  based on various fitted model objects (`fitdistr`, `glm`, `hurdle`,
  `zeroinfl`, `zerotrunc`, ...)

* `FLXMRnegbin()` driver for estimating mixtures of negative binomial models
  using `flexmix`.

* `MBnegbin()`, `MBztnegbin()`, `MBztpoisson()` and `MBbinomial()` families for estimating 
  boosted components of hurdle models using `mboost`.


# countreg 0.1-0

* Ported `zeroinfl()` and `hurdle()` along with corresponding vignette
  from `pscl` to `countreg`.

* Added `zerotrunc()` for estimating zero-truncated count regressions.
