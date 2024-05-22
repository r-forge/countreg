tol <- 1e-6

## data
data("CrabSatellites", package = "countreg")
cs <- CrabSatellites[, c("satellites", "width", "color")]
cs$color <- as.numeric(cs$color)
cs <- subset(cs, subset = satellites > 0)
form <- satellites ~ .

# Reference values for coefs etc. of hurdle models and specs for fitting
ref_list <- list(# poisson
  zt_p = list(formula = form,
              coef = c(`(Intercept)` = 0.56269906967054, width = 0.0342378636202522, 
                       color = 0.00716550359326728),
              se = c(`(Intercept)` = 0.645438550516392, width = 0.0222274963924805, 
                     color = 0.0666273073717945),
              ll = structure(-267.547433720941, df = 3L, nobs = 111L, class = "logLik"),
              dist = "poisson"),
  # negbin
  zt_nb2 = list(formula = form,
                coef = c(`(Intercept)` = 0.427223921701066, width = 0.0378896014648731, 
                         color = 0.00698505703049874),
                se = c(`(Intercept)` = 0.941130550556333, width = 0.0327507718574901, 
                       color = 0.0910813063753843),
                theta = 4.60546043331554,
                se_ltheta = 0.352936821886573,
                ll = structure(-255.802166814923, df = 4L, nobs = 111L, class = "logLik"),
                dist = "negbin"),
  # geometric
  zt_g = list(formula = form,
              coef = c(`(Intercept)` = 0.0440923450112795, width = 0.0446097834745566, 
                       color = 0.00767870634829046),
              se = c(`(Intercept)` = 1.55756775337772, width = 0.0546318904213941, 
                     color = 0.142937746716148),
              ll = structure(-265.623747300912, df = 3L, nobs = 111L, class = "logLik"),
              dist = "geometric")
)



## Helpers
test_mean <- function(model, data, muX, p0, tol, add_info) {
  ## Use 'textbook' formulas for computing means, different to implementation
  # zerotruncated mean
  mu <- muX/(1 - p0)
  expect_equivalent(predict(model, newdata = data, type = "response"), mu, tol = tol,
                    info = paste("Mean prediction of", add_info))
  
  # count mean
  expect_equivalent(predict(model, newdata = data, type = "count"),
                    muX, tol = tol, info = paste("Count mean prediction of", add_info))
  
  # fitted and prediction
  expect_equivalent(predict(model, newdata = data, type = "response"), fitted(model),
                    tol = tol, info = paste("Fitted and mean prediction of", add_info))
}

test_pred <- function(model, form, data, Y, dist, tol, add_info) {
  require("topmodels")
  
  mtX <- terms(form, data = data, rhs = 1L)
  X <- model.matrix(mtX, data)
  lin_pred <- X %*% coef(model)
  
  # mean of distribution
  muX <- as.numeric(exp(lin_pred))
  
  f_count <- function(x) {
    return(switch(dist,
                  "poisson" = dpois(x, lambda = muX),
                  "negbin" = dnbinom(x, size = model$theta, mu = muX),
                  "geometric" = dnbinom(x, size = 1, mu = muX)
    ))
  }
  
  # testing only mean predictions
  test_mean(model, data, muX, f_count(0), tol, add_info)
  
  # other predictions
  y_unique <- 1:max(Y)
  probs <- matrix(NA, nrow = length(muX), ncol = length(y_unique))
  
  ## Use textbook formulas again for testing
  for (i in y_unique) probs[,i] <- f_count(i)/(1 - f_count(0))
  
  expect_equivalent(predict(model, data, type = "prob", at = y_unique),
                    probs, tol,
                    info = paste("cdf predictions of", add_info))
  expect_equivalent(as.matrix(procast(model, data, type = "density", at = y_unique)),
                    probs, tol, info = paste("topmodels probability predictions of",
                                             add_info))
  
  # zero prob
  expect_equivalent(predict(model, newdata = data, type = "zero"),
                    1 - f_count(0), tol, info = paste("Zero prob prediction of", add_info))
  
}


## Compare with reference values and test predictions
for (ref in ref_list) {
  mod <- zerotrunc(ref[["formula"]], data = cs, dist = ref[["dist"]])
  mod_name <- paste("zerotrunc" , ref[["dist"]])
  
  # Compare with reference values
  expect_equal(coef(mod), ref[["coef"]], tol, info = paste("Compare coefs to reference for", mod_name))
  expect_equal(sqrt(diag(vcov(mod))), ref[["se"]], tol, info = paste("Compare SEs to reference for", mod_name))
  expect_equal(logLik(mod), ref[["ll"]], tol, info = paste("Compare logLik to reference for", mod_name))
  if (ref[["dist"]] == "negbin") {
    expect_equal(mod$theta, ref[["theta"]], tol, info = paste("Compare theta to reference for", mod_name))
    expect_equal(mod$SE.logtheta, ref[["se_ltheta"]], tol, info = paste("Compare SEs of logtheta to reference for", mod_name))
  }
  
  # Test predictions
  test_pred(mod, ref[["formula"]], cs, cs$satellites, ref[["dist"]], tol,
            add_info = mod_name)
}
