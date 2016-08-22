## binom: Binomial
sbinom <- function(x, prob, size, parameter = "prob", drop = TRUE) {
  s <- x/prob - (size - x)/(1 - prob)
  if(drop) s else cbind("prob" = s)
}

## nbinom: Negative binomial
snbinom <- function(x, mu, size, parameter = c("mu", "size"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "size")))
  s <- cbind(
    if("mu" %in% parameter) x/mu - (x + size)/(mu + size) else NULL,
    if("size" %in% parameter) digamma(x + size) - digamma(size) +
      log(size) + 1 - log(mu + size) - (x + size) / (mu + size) else NULL
  )
  colnames(s) <- c("mu", "size")[c("mu", "size") %in% parameter]
  if(drop) drop(s) else s
}

## ## ztpoisson: zero-truncated Poisson
## ## -> more flexible implementation in ztpois.R now
## dztpois <- function(x, lambda, log = FALSE) {
##   rval <- dpois(x, lambda = lambda, log = TRUE) -
##     ppois(0, lambda = lambda, lower.tail = FALSE, log.p = TRUE)
##   if(log) rval else exp(rval)
## }
## 
## sztpois <- function(x, lambda, parameter = "lambda", drop = TRUE) {
##   s <- x/lambda - 1 - exp(-lambda)/(1 - exp(-lambda))
##   if(drop) s else cbind("lambda" = s)
## }

## ztnbinom: Zero-truncated negative binomial
dztnbinom <- function(x, mu, size, log = FALSE) {
  rval <- dnbinom(x, mu = mu, size = size, log = TRUE) -
    pnbinom(0, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
  if(log) rval else exp(rval)
}

sztnbinom <- function(x, mu, size, parameter = c("mu", "size"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "size")))
  s <- snbinom(x, mu = mu, size = size, parameter = parameter, drop = FALSE)
  logratio <- pnbinom(0, mu = mu, size = size, log.p = TRUE) -
    pnbinom(0, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
  if("mu" %in% parameter) s[, "mu"] <- s[, "mu"] -
    exp(logratio + log(size) - log(mu + size))
  if("size" %in% parameter) s[, "size"] <- s[, "size"] +
    exp(logratio) * (log(size) - log(mu + size) + 1 - size/(mu + size))
  if(drop) drop(s) else s
}

