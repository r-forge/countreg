## Comparisons against reference
library(Formula)
tol <- 1e-6

data("CrabSatellites")
form <- satellites ~ width + color
form_h <- satellites ~ width + color | weight
link <- "log"
link_theta <- "log"


# NB2
fm_nb2 <- nbreg(form, data = CrabSatellites, link = link, link.theta = link_theta)
coef_ref <- c(`mu_(Intercept)` = -3.68546352253491, mu_width = 0.17839041962828,
              mu_color.L = -0.414226407908124, mu_color.Q = 0.130116604025637,
              mu_color.C = 0.0440867559068667, `theta_(Intercept)` = -0.0703175555237518)
se_ref <- c(`mu_(Intercept)` = 1.26893204172812, mu_width = 0.0480652291478379,
            mu_color.L = 0.294290216100736, mu_color.Q = 0.242444381199419,
            mu_color.C = 0.178874473742344, `theta_(Intercept)` = 0.180284018618521)
ll_ref <- structure(-374.297920053537, df = 6L, nobs = 173L, class = "logLik")
expect_equal(coef(fm_nb2), coef_ref, tol, info = "Compare coefs to reference for nbreg NB2")
expect_equal(sqrt(diag(vcov(fm_nb2))), se_ref, tol, info = "Compare SEs to reference for nbreg NB2")
expect_equal(logLik(fm_nb2), ll_ref, tol, info = "Compare logLik to reference for nbreg NB2")

fm_nb2_lower <- nbreg(form, data = CrabSatellites, dist = "nb2", link = link,
                      link.theta = link_theta)
expect_identical(coef(fm_nb2), coef(fm_nb2_lower),
                 info = "Check if coefs of lowercase nb2 are identical")
expect_identical(vcov(fm_nb2), vcov(fm_nb2_lower),
                 info = "Check if vcov of lowercase nb2 is identical")
expect_identical(logLik(fm_nb2), logLik(fm_nb2_lower),
                 info = "Check if loglik of lowercase nb2 is identical")
expect_identical(residuals(fm_nb2, "pearson"), residuals(fm_nb2_lower, "pearson"),
                 info = "Check if pearson resids of lowercase nb2 are identical")


# NBH
fm_nbh <- nbreg(form_h, data = CrabSatellites, dist = "NB2", link = link,
                link.theta = link_theta)
coef_ref <- c(`mu_(Intercept)` = -2.2440670873349, mu_width = 0.125911269094792, 
              mu_color.L = -0.259193239569225, mu_color.Q = 0.0310068596163533, 
              mu_color.C = 0.0326075887875688, `theta_(Intercept)` = -3.97309995040085, 
              theta_weight = 1.55635228341314)
se_ref <- c(`mu_(Intercept)` = 1.11648575682368, mu_width = 0.040472236809787, 
            mu_color.L = 0.29064964663664, mu_color.Q = 0.238930308321897, 
            mu_color.C = 0.175702134821161, `theta_(Intercept)` = 0.870539633894667, 
            theta_weight = 0.349940939749115)
ll_ref <- structure(-363.505261263564, df = 7L, nobs = 173L, class = "logLik")
expect_equal(coef(fm_nbh), coef_ref, tol, info = "Compare coefs to reference for nbreg NBH")
expect_equal(sqrt(diag(vcov(fm_nbh))), se_ref, tol, info = "Compare SEs to reference for nbreg NBH")
expect_equal(logLik(fm_nbh), ll_ref, tol, info = "Compare logLik to reference for nbreg NBH")


# NB1
fm_nb1 <- nbreg(form, data = CrabSatellites, dist = "NB1", link = link,
                link.theta = link_theta)
coef_ref <- c(`mu_(Intercept)` = -3.76573267042693, mu_width = 0.177062322435001,
              mu_color.L = -0.640016895815285, mu_color.Q = -0.193085157579817,
              mu_color.C = -0.0692478243188591, `theta_(Intercept)` = -1.14820467167638)
se_ref <- c(`mu_(Intercept)` = 0.967179850729637, mu_width = 0.035764633133798,
            mu_color.L = 0.306204683640572, mu_color.Q = 0.242492659108847,
            mu_color.C = 0.167208012013685, `theta_(Intercept)` = 0.195036632753825)
ll_ref <- structure(-365.354036365778, df = 6L, nobs = 173L, class = "logLik")
expect_equal(coef(fm_nb1), coef_ref, tol, info = "Compare coefs to reference for nbreg NB1")
expect_equal(sqrt(diag(vcov(fm_nb1))), se_ref, tol, info = "Compare SEs to reference for nbreg NB1")
expect_equal(logLik(fm_nb1), ll_ref, tol, info = "Compare logLik to reference for nbreg NB1")

fm_nb1_lower <- nbreg(form, data = CrabSatellites, dist = "nb1", link = link,
                      link.theta = link_theta)
expect_identical(coef(fm_nb1), coef(fm_nb1_lower),
                 info = "Check if coefs of lowercase nb1 are identical")
expect_identical(vcov(fm_nb1), vcov(fm_nb1_lower),
                 info = "Check if vcov of lowercase nb1 is identical")
expect_identical(logLik(fm_nb1), logLik(fm_nb1_lower),
                 info = "Check if loglik of lowercase nb1 is identical")
expect_identical(residuals(fm_nb1, "pearson"), residuals(fm_nb1_lower, "pearson"),
                 info = "Check if pearson resids of lowercase nb1 are identical")


# NB1-H
fm_nb1h <- nbreg(form_h, data = CrabSatellites, dist = "NB1", link = link,
                 link.theta = link_theta)
coef_ref <- c(`mu_(Intercept)` = -1.97326851123967, mu_width = 0.112488668821603, 
              mu_color.L = -0.552411887406806, mu_color.Q = -0.220354983473709, 
              mu_color.C = -0.0634412493751442, `theta_(Intercept)` = -3.93755772164971, 
              theta_weight = 1.09665147594751)
se_ref <- c(`mu_(Intercept)` = 1.08226233020628, mu_width = 0.0390075410880065, 
            mu_color.L = 0.302060703169686, mu_color.Q = 0.236700736166993, 
            mu_color.C = 0.160739400248805, `theta_(Intercept)` = 1.00554524213124, 
            theta_weight = 0.387054747389924)
ll_ref <- structure(-361.298100848771, df = 7L, nobs = 173L, class = "logLik")
expect_equal(coef(fm_nb1h), coef_ref, tol, info = "Compare coefs to reference for nbreg NB1-H")
expect_equal(sqrt(diag(vcov(fm_nb1h))), se_ref, tol, info = "Compare SEs to reference for nbreg NB1-H")
expect_equal(logLik(fm_nb1h), ll_ref, tol, info = "Compare logLik to reference for nbreg NB1-H")


## Test theta = inf
fm_nb2_inf <- nbreg(form, data = CrabSatellites, link = link,
                    link.theta = link_theta, theta = Inf)
fm_nb1_inf <- nbreg(form, data = CrabSatellites, dist = "NB1", link = link,
                    link.theta = link_theta, theta = Inf)
fm_pois <- glm(form, data = CrabSatellites, family = poisson())

test_theta_inf <- function(nb_model, poisson_model, tol) {
  dist <- nb_model$dist
  stopifnot(dist == "NB2" || dist == "NB1")
  
  expect_equivalent(coef(nb_model)[names(coef(nb_model)) != "theta_(Intercept)"],
                    coef(poisson_model), tol,
                    info = paste("Compare coefs of nbreg", dist, "with theta = Inf to Poisson"))
  expect_equivalent(logLik(nb_model), logLik(poisson_model), tol,
                    info = paste("Compare logLik of nbreg", dist, "with theta = Inf to Poisson"))
  expect_true(all(is.finite(vcov(nb_model))),
              info = paste("Finite vcov for nbreg", dist, "with theta = Inf"))
  
}

test_theta_inf(fm_nb2_inf, fm_pois, tol)
test_theta_inf(fm_nb1_inf, fm_pois, tol)

## Test predictions
test_mean <- function(model, data, reference, add_info) {
  expect_equal(predict(model, newdata = data, type = "response"),
               reference, tol = tol,
               info = paste("Mean prediction of", add_info))
  expect_equivalent(predict(model, newdata = data, type = "response"), fitted(model),
                    info = paste("Fitted and mean prediction of", add_info))
}

test_pred <- function(model, form, data, Y, link, link_theta, tol, add_info){
  dist <- model$dist
  stopifnot(dist == "NB2" || dist == "NB1")
  
  form <- as.Formula(form)
  
  # constant theta if no covars
  if (length(form)[2L] < 2L) form <- as.Formula(formula(form), ~ 1)
  
  mtX <- terms(form, data = data, rhs = 1L)
  mtZ <- delete.response(terms(form, data = data, rhs = 2L))
  
  X <- model.matrix(mtX, data)
  Z <- model.matrix(mtZ, data)
  
  lin_pred <- X %*% coef(model, "mu")
  lin_pred_Z <- Z %*% coef(model, "theta")
  
  if (link == "log") {
    linkinv <- exp
  } else {
    stop(paste(link, "link not supported."))
  }
  
  if (link_theta == "log") {
    linkinv_theta <- exp
  } else {
    stop(paste(link_theta, "link not supported."))
  }
  
  # distribution parameters
  mu <- as.numeric(linkinv(lin_pred))
  theta <- as.numeric(linkinv(lin_pred_Z))
  
  # testing only mean predictions
  test_mean(model, data, mu, add_info)
  
  # other predictions
  y_unique <- 0:max(Y)
  probs <- matrix(NA, nrow = length(mu), ncol = length(y_unique))
  size <- if (dist == "NB2") theta else mu * theta
  for (i in seq_along(y_unique)) probs[,i] <- dnbinom(y_unique[i], mu = mu,
                                                      size = size)
  
  expect_equivalent(predict(model, data, type = "prob"), probs, tol,
                    info = paste("Probability predictions of", add_info))
  expect_equal(predict(model, data, type = "theta"), theta, tol,
               info = paste("Predicted theta of", add_info))
  expect_equal(predict(model, data, type = "parameters"),
               data.frame(mu = mu, theta = theta), tol,
               info = paste("Predicted parameters of", add_info))
}

test_pred(fm_nb2, form, CrabSatellites, CrabSatellites$satellites,
          link, link_theta, tol, add_info = "nbreg NB2")
test_pred(fm_nb1, form, CrabSatellites, CrabSatellites$satellites,
          link, link_theta, tol, add_info = "nbreg NB1")
test_pred(fm_nbh, form_h, CrabSatellites, CrabSatellites$satellites,
          link, link_theta, tol, add_info = "nbreg NBH")
test_pred(fm_nb1h, form_h, CrabSatellites, CrabSatellites$satellites,
          link, link_theta, tol, add_info = "nbreg NB1-H")

