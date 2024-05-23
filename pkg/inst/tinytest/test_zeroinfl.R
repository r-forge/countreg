## Comparisons against reference
tol <- 1e-6

data("CrabSatellites")
form <- satellites ~ width + color | width + color
form2 <- satellites ~ width + color | 1

# Reference values for coefs etc. of hurdle models and specs for fitting
ref_list <- list(# logit-poisson
  logit_p = list(formula = form,
                 coef = c(`count_(Intercept)` = 0.495696172809289, count_width = 0.0409135634370457,
                          count_color.L = 0.110499589812118, count_color.Q = 0.392696945634541,
                          count_color.C = 0.165903906770087, `zero_(Intercept)` = 11.7569805701059,
                          zero_width = -0.469398421266328, zero_color.L = 0.967899515144614,
                          zero_color.Q = 0.636388127651717, zero_color.C = 0.111829618595148),
                 se = c(`count_(Intercept)` = 0.597262787683302, count_width = 0.0221520278794397,
                        count_color.L = 0.143989746773683, count_color.Q = 0.120965061737862,
                        count_color.C = 0.0929786820242997, `zero_(Intercept)` = 2.81977803725436,
                        zero_width = 0.108636502110447, zero_color.L = 0.587308567110469,
                        zero_color.Q = 0.483109131180914, zero_color.C = 0.347697538367407),
                 ll = structure(-355.748785994805, df = 10L, nobs = 173L, class = "logLik"),
                 dist = "poisson",
                 link = "logit"),
  chauchit_p = list(formula = form,
                    coef = c(`count_(Intercept)` = 0.494162395584611, count_width = 0.0409821601062638, 
                             count_color.L = 0.11174511314375, count_color.Q = 0.391416562646885, 
                             count_color.C = 0.165125584147918, `zero_(Intercept)` = 11.7669159807711, 
                             zero_width = -0.469276727647141, zero_color.L = 1.27336889369364, 
                             zero_color.Q = 0.569027598291114, zero_color.C = 0.128053804197964),
                    se = c(`count_(Intercept)` = 0.596710674518854, count_width = 0.0221313126663705, 
                           count_color.L = 0.144101221423316, count_color.Q = 0.120942061321052, 
                           count_color.C = 0.0928193621041184, `zero_(Intercept)` = 3.46316085110373, 
                           zero_width = 0.13454780160465, zero_color.L = 0.750950352012352, 
                           zero_color.Q = 0.590637685159934, zero_color.C = 0.374938171937907),
                    ll = structure(-356.34354972, df = 10L, nobs = 173L, class = "logLik"),
                    dist = "poisson",
                    link = "cauchit"),
  probit_nb = list(formula = form,
                   coef = c(`count_(Intercept)` = 0.327220388871631, count_width = 0.0463760066420507,
                            count_color.L = 0.109193394790255, count_color.Q = 0.405597152621354,
                            count_color.C = 0.160719691315573, `zero_(Intercept)` = 7.159944791039,
                            zero_width = -0.28787412648044, zero_color.L = 0.584267443159918,
                            zero_color.Q = 0.419675886613702, zero_color.C = 0.0672576359506783),
                   se = c(`count_(Intercept)` = 0.833551309645838, count_width = 0.0309583216054862,
                          count_color.L = 0.207523792419349, count_color.Q = 0.171966740456552,
                          count_color.C = 0.125850698520329, `zero_(Intercept)` = 1.742842071604,
                          zero_width = 0.0672479947806656, zero_color.L = 0.360016422523919,
                          zero_color.Q = 0.297708130893725, zero_color.C = 0.221419272958228),
                   theta = 5.71132861168444,
                   se_ltheta = 0.373409177087442,
                   ll = structure(-346.750439059584, df = 11L, nobs = 173L, class = "logLik"),
                   dist = "negbin",
                   link = "probit"),
  cloglog_g = list(formula = form2,
                   coef = c(`count_(Intercept)` = -2.94572460920356, count_width = 0.154956965438731, 
                            count_color.L = -0.359988897127733, count_color.Q = 0.183185061777129, 
                            count_color.C = 0.0635082037356864, `zero_(Intercept)` = -2.22312182406824),
                   se = c(`count_(Intercept)` = 1.37962852347718, count_width = 0.0513741323249627, 
                          count_color.L = 0.316729731819911, count_color.Q = 0.259299954460248, 
                          count_color.C = 0.187171779126418, `zero_(Intercept)` = 0.617341672371466),
                   ll = structure(-372.907554863153, df = 6L, nobs = 173L, class = "logLik"),
                   dist = "geometric",
                   link = "cloglog")
)



## Helpers
test_mean <- function(model, data, muX, p0_zero, p0_count, tol, add_info) {
  ## Use 'textbook' formulas for computing means, different to implementation
  # hurdle mean
  mu <- (1 - p0_zero) * muX
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
  ## Raises warning when zero-part is just a constant
  ## because of predict.zeroinfl(..., model = "zero")
  expect_equivalent(predict(model, newdata = data, type = "mean", model = "zero"),
                    p0_zero, tol, info = paste("Zero mean prediction of", add_info))
  
}

test_pred <- function(model, form, data, Y, link, count_dist, tol, add_info) {
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
  
  linkinv <- if (link == "logit") {
    plogis
  } else if (link == "probit") {
    pnorm
  } else if (link == "log") {
    exp
  } else if (link == "cloglog") {
    function(x) -expm1(-exp(x)) 
  } else if (link == "cauchit") {
    pcauchy 
  } else {
    stop(paste(link, "link not supported."))
  }
  
  # distribution parameters
  muX <- as.numeric(exp(lin_pred))
  muZ <- as.numeric(linkinv(lin_pred_Z)) # prob for zero
  
  
  f_count <- function(x) {
    return(switch(count_dist,
                  "poisson" = dpois(x, lambda = muX),
                  "negbin" = dnbinom(x, size = model$theta, mu = muX),
                  "geometric" = dnbinom(x, size = 1, mu = muX)))
  }
  
  
  # testing only mean predictions
  test_mean(model, data, muX, muZ, f_count(0), tol, add_info)
  
  # other predictions
  y_unique <- 0:max(Y)
  probs <- matrix(NA, nrow = length(muX), ncol = length(y_unique))
  
  for (i in seq_along(y_unique)) {
    ## Use textbook formulas again for testing
    if (i > 1) {
      probs[,i] <- (1 - muZ) * f_count(i - 1)
    } else {
      probs[,i] <- muZ + (1 - muZ)*f_count(0)
    }
  }
  
  expect_equivalent(predict(model, data, type = "prob", at = y_unique),
                    t(apply(probs, MARGIN = 1, cumsum)), tol,
                    info = paste("cdf predictions of", add_info))
  
  ## Raises warnings when zero part is a constant
  ## because of predict.zeroinfl(..., type = "zero") used internally in prodist.zeroinfl
  expect_equivalent(as.matrix(procast(model, data, type = "density", at = y_unique)),
                    probs, tol, info = paste("topmodels probability predictions of",
                                             add_info))
  
}


## Compare with reference values and test predictions
for (ref in ref_list) {
  mod <- zeroinfl(ref[["formula"]], data = CrabSatellites, dist = ref[["dist"]],
                  link = ref[["link"]])
  mod_name <- paste("zeroinfl" , paste(ref[["link"]], ref[["dist"]], sep = "-"))
  
  # Compare with reference values
  expect_equal(coef(mod), ref[["coef"]], tol, info = paste("Compare coefs to reference for", mod_name))
  expect_equal(sqrt(diag(vcov(mod))), ref[["se"]], tol, info = paste("Compare SEs to reference for", mod_name))
  expect_equal(logLik(mod), ref[["ll"]], tol, info = paste("Compare logLik to reference for", mod_name))
  if (ref[["dist"]] == "negbin") {
    expect_equal(mod$theta, ref[["theta"]], tol, info = paste("Compare theta to reference for", mod_name))
    expect_equal(mod$SE.logtheta, ref[["se_ltheta"]], tol, info = paste("Compare SEs of logtheta to reference for", mod_name))
  }
  
  # Test predictions
  test_pred(mod, ref[["formula"]], CrabSatellites, CrabSatellites$satellites, ref[["link"]],
            ref[["dist"]], tol, add_info = mod_name)
}
