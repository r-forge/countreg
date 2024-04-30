# Comparisons against reference
library(Formula)
tol <- 1e-6

data("CrabSatellites")
form <- satellites ~ width + color
link <- "log"
link_theta <- "log"

## TODO : Add NB1 and NBH model

fm_nb2 <- nbreg(form, data = CrabSatellites, link = link, link.theta = link_theta)

coef_ref <- c(`mu_(Intercept)` = -3.68546352253491, mu_width = 0.17839041962828,
              mu_color.L = -0.414226407908124, mu_color.Q = 0.130116604025637,
              mu_color.C = 0.0440867559068667, `theta_(Intercept)` = -0.0703175555237518)
se_ref <- c(`(Intercept)` = 1.26893204172812, width = 0.0480652291478379,
            color.L = 0.294290216100736, color.Q = 0.242444381199419,
            color.C = 0.178874473742344, `(Intercept)` = 0.180284018618521)
ll_ref <- structure(-374.297920053537, df = 6L, nobs = 173L, class = "logLik")

expect_equal(coef(fm_nb2), coef_ref, tol, info = "Compare coefs to reference")
expect_equal(sqrt(diag(vcov(fm_nb2))), se_ref, tol, info = "Compare SEs to reference")
expect_equal(logLik(fm_nb2), ll_ref, tol, info = "Compare logLik to reference")


# Test predictions
test_mean <- function(model, data, reference, add_info) {
  expect_equal(predict(model, newdata = data, type = "response"),
               reference, tol = tol,
               info = paste("Mean prediction of", add_info, "."))
  expect_equivalent(predict(model, newdata = data, type = "response"), fitted(model),
                    info = paste("Fitted and mean prediction of", add_info))
}

test_pred <- function(model, form, data, Y, link, link_theta, tol, add_info){
  dist <- model$dist
  
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
  for (i in seq_along(y_unique)) probs[,i] <- dnbinom(y_unique[i], mu = mu,
                                                      size = theta)
  
  expect_equivalent(predict(model, data, type = "prob"), probs, tol,
                    info = paste("Probability predictions of", add_info))
  expect_equal(predict(model, data, type = "theta"), theta, tol,
               info = paste("Predicticted theta of", add_info))
  expect_equal(predict(model, data, type = "parameters"),
               data.frame(mu = mu, theta = theta), tol,
               info = paste("Predicticted parameters of", add_info))
}

test_pred(fm_nb2, form, CrabSatellites, CrabSatellites$satellites,
          link, link_theta, tol, add_info = paste("nbreg NB2"))
