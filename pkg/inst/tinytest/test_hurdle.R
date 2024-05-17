## Comparisons against reference
library(Formula)
tol <- 1e-6

data("CrabSatellites")
form <- satellites ~ width + color | width + color

fm_hp1 <- hurdle(form, data = CrabSatellites, dist = "poisson")
coef_ref <- dput(coef(fm_hp1))
se_ref <- c(`count_(Intercept)` = 0.600167451122502, count_width = 0.0222775187986676, 
            count_color.L = 0.144023240226206, count_color.Q = 0.121096067125091, 
            count_color.C = 0.0936545055002974, `zero_(Intercept)` = 2.74141393130641, 
            zero_width = 0.105528725732993, zero_color.L = 0.580325262586876, 
            zero_color.Q = 0.474722800302477, zero_color.C = 0.337382299253625)
ll_ref <- structure(-355.646193700518, df = 10L, nobs = 173L, class = "logLik")
expect_equal(coef(fm_hp1), coef_ref, tol, info = "Compare coefs to reference for hurdle logit-poisson")
expect_equal(sqrt(diag(vcov(fm_hp1))), se_ref, tol, info = "Compare SEs to reference for hurdle logit-poisson")
expect_equal(logLik(fm_hp1), ll_ref, tol, info = "Compare logLik to reference for hurdle logit-poisson")


## Test single-part formula = using same covariates in both parts
fm_hp11 <- hurdle(satellites ~ width + color, data = CrabSatellites, dist = "poisson")
expect_identical(coef(fm_hp11), coef(fm_hp1),
                 info = "Check if coefs from single-part formula are identical")
expect_identical(vcov(fm_hp11), vcov(fm_hp1),
                 info = "Check if vcov from single-part formula is identical")
expect_identical(logLik(fm_hp11), logLik(fm_hp1),
                 info = "Check if loglik from single-part formula is identical")
expect_identical(residuals(fm_hp11, "pearson"), residuals(fm_hp1, "pearson"),
                 info = "Check if pearson resids from single-part formula are identical")


## Test predictions
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
                    # "poisson" = dpois(0, lambda = muZ),
                    # "negbin" = dnbinom(0, size = model$theta["zero"], mu = muZ),
                    # "geometric" = dnbinom(0, size = 1, mu = muZ)
                    ))
    }
  }
  
  
  f_count <- function(x) {
    return(switch(count_dist,
                  "poisson" = dpois(x, lambda = muX),
                  # "negbin" = dnbinom(x, size = model$theta["count"], mu = muX),
                  # "geometric" = dnbinom(x, size = 1, mu = muX),
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

test_pred(fm_hp1, form, CrabSatellites, CrabSatellites$satellites, "logit",
          "binomial", "poisson", tol, add_info = "hurdle logit-poisson")
