# Comparisons against reference
tol <- 1e-6

data("CrabSatellites")
fm_nb2 <- nbreg(satellites ~ width + color, data = CrabSatellites)

coef_ref <- c(`mu_(Intercept)` = -3.68546352253491, mu_width = 0.17839041962828, 
              mu_color.L = -0.414226407908124, mu_color.Q = 0.130116604025637, 
              mu_color.C = 0.0440867559068667, `theta_(Intercept)` = -0.0703175555237518)
se_ref <- c(`(Intercept)` = 1.26893204172812, width = 0.0480652291478379, 
            color.L = 0.294290216100736, color.Q = 0.242444381199419,
            color.C = 0.178874473742344, `(Intercept)` = 0.180284018618521)
ll_ref <- structure(-374.297920053537, df = 6L, nobs = 173L, class = "logLik")

expect_equal(coef(fm_nb2), coef_ref, tol, info = "Compare coefs to reference.")
expect_equal(sqrt(diag(vcov(fm_nb2))), se_ref, tol, info = "Compare SEs to reference.")
expect_equal(logLik(fm_nb2), ll_ref, tol, info = "Compare logLik to reference.")
