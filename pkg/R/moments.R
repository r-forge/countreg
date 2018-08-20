## binom: Binomial
mean_binom <- function(prob, size, drop = TRUE) {
    if(drop) prob * size else cbind("mean" = prob * size)
}

var_binom <- function(prob, size, drop = TRUE) {
    if(drop) size * prob * (1 - prob) else cbind("var" = size * prob * (1 - prob))
}


## pois: Poisson
mean_pois <- function(lambda, drop = TRUE) {
    if(drop) lambda else cbind("mean" = lambda)
}

var_pois <- function(lambda, drop = TRUE) {
    if(drop) lambda else cbind("var" = lambda)
}

## nbinom: Negative binomial
mean_nbinom <- function(mu, size, drop = TRUE) {
    if(drop) mu else cbind("mean" = mu)
}

var_nbinom <- function(mu, size, drop = TRUE) {
    if(drop) mu + mu^2 / size else cbind("var" = mu + mu^2 / size)
}

