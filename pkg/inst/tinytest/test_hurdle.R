## Comparisons against reference
library(Formula)
tol <- 1e-6

# increased tolerance for comparing theta and SE of logtheta in nb2_nb2 because
# the zero-hurdle in nb2_nb2 is poorly conditioned
tol_theta <- tol * 100
tol_se_ltheta <- tol * 1000

data("CrabSatellites")
form <- satellites ~ width + color | width + color

# Reference values for coefs etc. of hurdle models and specs for fitting
ref_list <- list(# logit-poisson
  logit_p = list(formula = form,
                 coef = c(`count_(Intercept)` = 0.52667131775527, count_width = 0.0397109265472577, 
                          count_color.L = 0.107779650124439, count_color.Q = 0.394877726205593, 
                          count_color.C = 0.169369411874262, `zero_(Intercept)` = -11.7555177268425, 
                          zero_width = 0.467955985251572, zero_color.L = -0.958372471323788, 
                          zero_color.Q = -0.589269206154217, zero_color.C = -0.098672167212181),
                 se = c(`count_(Intercept)` = 0.600167451122502, count_width = 0.0222775187986676, 
                        count_color.L = 0.144023240226206, count_color.Q = 0.121096067125091, 
                        count_color.C = 0.0936545055002974, `zero_(Intercept)` = 2.74141393130641, 
                        zero_width = 0.105528725732993, zero_color.L = 0.580325262586876, 
                        zero_color.Q = 0.474722800302477, zero_color.C = 0.337382299253625),
                 ll = structure(-355.646193700518, df = 10L, nobs = 173L, class = "logLik"),
                 dist = "poisson",
                 zero = "binomial",
                 link = "logit"),
  # geometric-poisson
  g_p = list(formula = form,
             coef = c(`count_(Intercept)` = 0.52667131775527, count_width = 0.0397109265472577, 
                      count_color.L = 0.107779650124439, count_color.Q = 0.394877726205593, 
                      count_color.C = 0.169369411874262, `zero_(Intercept)` = -11.7555177268425, 
                      zero_width = 0.467955985251572, zero_color.L = -0.958372471323788, 
                      zero_color.Q = -0.589269206154217, zero_color.C = -0.098672167212181),
             se = c(`count_(Intercept)` = 0.600167451122502, count_width = 0.0222775187986676, 
                    count_color.L = 0.144023240226206, count_color.Q = 0.121096067125091, 
                    count_color.C = 0.0936545055002974, `zero_(Intercept)` = 2.74141393130281, 
                    zero_width = 0.105528725732855, zero_color.L = 0.58032526258691, 
                    zero_color.Q = 0.474722800302479, zero_color.C = 0.337382299253634),
             ll = structure(-355.646193700518, df = 10L, nobs = 173L, class = "logLik"),
             dist = "poisson",
             zero = "geometric",
             link = "log"),
  # poisson-geometric
  p_g = list(formula = form,
             coef = c(`count_(Intercept)` = 0.0651428428002858, count_width = 0.0487416343129737, 
                      count_color.L = 0.1062110550535, count_color.Q = 0.473900811521687, 
                      count_color.C = 0.188547027992332, `zero_(Intercept)` = -7.77983128190527, 
                      zero_width = 0.292745587161981, zero_color.L = -0.609317402102402, 
                      zero_color.Q = -0.473857184589676, zero_color.C = -0.142615188467327),
             se = c(`count_(Intercept)` = 1.48668492513064, count_width = 0.0555008761016575, 
                    count_color.L = 0.376934554794095, count_color.Q = 0.308010553377792, 
                    count_color.C = 0.218467175134894, `zero_(Intercept)` = 1.68297368011838, 
                    zero_width = 0.0635598044517968, zero_color.L = 0.368202398547663, 
                    zero_color.Q = 0.297330489665648, zero_color.C = 0.210653530897402),
             ll = structure(-357.395233883897, df = 10L, nobs = 173L, class = "logLik"),
             dist = "geometric",
             zero = "poisson",
             link = "log"),
  # negbin-negbin (poorly conditioned zero-hurdle)
  nb2_nb2 = list(formula = form,
                 coef = c(`count_(Intercept)` = 0.43853962784324, count_width = 0.0420235024635511, 
                          count_color.L = 0.101310893405628, count_color.Q = 0.414370582115084, 
                          count_color.C = 0.171572744443123, `zero_(Intercept)` = -7.78032721575754, 
                          zero_width = 0.29276682441464, zero_color.L = -0.609379553658299, 
                          zero_color.Q = -0.473885818877084, zero_color.C = -0.142614297219935),
                 se = c(`count_(Intercept)` = 0.851461470439592, count_width = 0.0316794941991245, 
                        count_color.L = 0.20979469324798, count_color.Q = 0.174092532986972, 
                        count_color.C = 0.128649639952385, `zero_(Intercept)` = 1.68431774631175, 
                        zero_width = 0.0636257145756324, zero_color.L = 0.368258280965331, 
                        zero_color.Q = 0.29735794660727, zero_color.C = 0.210668310571957),
                 theta = c(count = 5.44999491584071, zero = 8608.19058947861),
                 se_ltheta = c(count = 0.377383705831662, zero = 140.219767434209),
                 ll = structure(-345.907624112724, df = 12L, nobs = 173L, class = "logLik"),
                 dist = "negbin",
                 zero = "negbin",
                 link = "log")
)



## Helpers
test_mean <- function(model, data, muX, p0_zero, p0_count, tol, add_info) {
  ## Use 'textbook' formulas for computing means, different to implementation
  # hurdle mean
  mu <- ((1 - p0_zero)/(1 - p0_count)) * muX
  expect_equivalent(predict(model, newdata = data, type = "response"), mu, tol = tol,
                    info = paste("Mean prediction of", add_info))
  
  # count mean
  expect_equivalent(predict(model, newdata = data, type = "mean", model = "count"),
                    muX, tol = tol, info = paste("Count mean prediction of", add_info))
  
  # zerotruncated mean
  expect_equivalent(predict(model, newdata = data, type = "mean", model = "truncated"),
                    muX/(1 - p0_count), tol = tol,
                    info = paste("Zerotrunc mean prediction of", add_info))
  
  # zero mean
  expect_equivalent(predict(model, newdata = data, type = "mean", model = "zero"),
                    1 - p0_zero, tol, info = paste("Zero mean prediction of", add_info))
  
  # fitted and prediction
  expect_equivalent(predict(model, newdata = data, type = "response"), fitted(model),
                    tol = tol, info = paste("Fitted and mean prediction of", add_info))
}

test_pred <- function(model, form, data, Y, link, zero_dist, count_dist, tol, add_info) {
  require("topmodels")
  require("Formula")
  
  form <- as.Formula(form)
  
  # use same covars for zero part if only one formula is specified
  if (length(form)[2L] < 2L) form[[3]][3] <- form[[3]][2]
  
  mtX <- terms(form, data = data, rhs = 1L)
  mtZ <- delete.response(terms(form, data = data, rhs = 2L))
  
  X <- model.matrix(mtX, data)
  Z <- model.matrix(mtZ, data)
  
  lin_pred <- X %*% coef(model, "count")
  lin_pred_Z <- Z %*% coef(model, "zero")
  
  if (link == "logit") {
    linkinv <- plogis
  } else if (link == "probit") {
    linkinv <- pnorm
  } else if (link == "log"){
    linkinv <- exp
  } else {
    stop(paste(link, "link not supported."))
  }
  
  # distribution parameters
  muX <- as.numeric(exp(lin_pred))
  muZ <- as.numeric(linkinv(lin_pred_Z))
  
  f_zero <- function(x) {
    if (x > 0) {
      return(0)
    } else {
      return(switch(zero_dist,
                    "binomial" = 1 - muZ,
                    "poisson" = dpois(0, lambda = muZ),
                    "negbin" = dnbinom(0, size = model$theta["zero"], mu = muZ),
                    "geometric" = dnbinom(0, size = 1, mu = muZ)
      ))
    }
  }
  
  
  f_count <- function(x) {
    return(switch(count_dist,
                  "poisson" = dpois(x, lambda = muX),
                  "negbin" = dnbinom(x, size = model$theta["count"], mu = muX),
                  "geometric" = dnbinom(x, size = 1, mu = muX),
                  # "binomial" = dbinom(x, prob = muX/model$size, size = model$size)
    ))
  }
  
  
  # testing only mean predictions
  test_mean(model, data, muX, f_zero(0), f_count(0), tol, add_info)
  
  # other predictions
  y_unique <- 0:max(Y)
  probs <- matrix(NA, nrow = length(muX), ncol = length(y_unique))
  
  for (i in seq_along(y_unique)) {
    ## Use textbook formulas again for testing
    if (i > 1) {
      probs[,i] <- ((1 - f_zero(0))/(1 - f_count(0))) * f_count(i - 1)
    } else {
      probs[,i] <- f_zero(0)
    }
  }
  
  expect_equivalent(predict(model, data, type = "prob", at = y_unique),
                    t(apply(probs, MARGIN = 1, cumsum)), tol,
                    info = paste("cdf predictions of", add_info))
  expect_equivalent(as.matrix(procast(model, data, type = "density", at = y_unique)),
                    probs, tol, info = paste("topmodels probability predictions of",
                                             add_info))
  
}


## Compare with reference values and test predictions
for (ref in ref_list) {
  mod <- hurdle(ref[["formula"]], data = CrabSatellites, dist = ref[["dist"]],
                zero = ref[["zero"]], link = ref[["link"]])
  mod_name <- paste("hurdle" , paste(if (ref[["zero"]] == "binomial") ref[["link"]] else ref[["zero"]], ref[["dist"]], sep = "-"))
  
  # Compare with reference values
  expect_equal(coef(mod), ref[["coef"]], tol, info = paste("Compare coefs to reference for", mod_name))
  expect_equal(sqrt(diag(vcov(mod))), ref[["se"]], tol, info = paste("Compare SEs to reference for", mod_name))
  expect_equal(logLik(mod), ref[["ll"]], tol, info = paste("Compare logLik to reference for", mod_name))
  if (ref[["dist"]] == "negbin" || ref[["zero"]] == "negbin") {
    expect_equal(mod$theta, ref[["theta"]], tol_theta, info = paste("Compare theta to reference for", mod_name))
    expect_equal(mod$SE.logtheta, ref[["se_ltheta"]], tol_se_ltheta, info = paste("Compare SEs of logtheta to reference for", mod_name))
  }
  
  # Test predictions
  test_pred(mod, ref[["formula"]], CrabSatellites, CrabSatellites$satellites, ref[["link"]],
            ref[["zero"]], ref[["dist"]], tol, add_info = mod_name)
}


## Test single-part formula = using same covariates in both parts

# logit-poisson
fm_hp1 <- hurdle(form, data = CrabSatellites, dist = "poisson")
fm_hp11 <- hurdle(satellites ~ width + color, data = CrabSatellites, dist = "poisson")
expect_identical(coef(fm_hp11), coef(fm_hp1),
                 info = "Check if coefs from single-part formula are identical")
expect_identical(vcov(fm_hp11), vcov(fm_hp1),
                 info = "Check if vcov from single-part formula is identical")
expect_identical(logLik(fm_hp11), logLik(fm_hp1),
                 info = "Check if loglik from single-part formula is identical")
expect_identical(residuals(fm_hp11, "pearson"), residuals(fm_hp1, "pearson"),
                 info = "Check if pearson resids from single-part formula are identical")


## Test if logit and geometric models are equivalent
# geometric-poisson
fm_hp2 <- hurdle(form, data = CrabSatellites, zero = "geometric", dist = "poisson")
expect_equal(coef(fm_hp2, model = "zero"), coef(fm_hp1, model = "zero"), tol,
             info = "logit and geometric zero models are identical")
