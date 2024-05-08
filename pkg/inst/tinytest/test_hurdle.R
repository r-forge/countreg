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
