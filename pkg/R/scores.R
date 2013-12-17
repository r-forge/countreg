## binom: Binomial
sbinom <- function(y, prob, size, parameter = "prob", drop = TRUE) {
  s <- lchoose(size, y) + y/prob - (size - y)/(1 - prob)
  if(drop) s else cbind("prob" = s)
}

## nbinom: Negative binomial
snbinom <- function(y, mu, size, parameter = c("mu", "size"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "size")))
  s <- cbind(
    if("mu" %in% parameter) y/mu - (y + size)/(mu + size) else NULL,
    if("size" %in% parameter) digamma(y + size) - digamma(size) +
      log(size) + 1 - log(mu + size) - (y + size) / (mu + size) else NULL
  )
  colnames(s) <- c("mu", "size")[c("mu", "size") %in% parameter]
  if(drop) drop(s) else s
}

## ztnbinom: Zero-truncated negative binomial
dztnbinom <- function(y, mu, size, log = FALSE) {
  rval <- dnbinom(y, mu = mu, size = size, log = TRUE) -
    pnbinom(0, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
  if(log) rval else exp(rval)
}

sztnbinom <- function(y, mu, size, parameter = c("mu", "size"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "size")))
  s <- snbinom(y, mu = mu, size = size, parameter = parameter, drop = FALSE)
  logratio <- pnbinom(0, mu = mu, size = size, log.p = TRUE) -
    pnbinom(0, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
  if("mu" %in% parameter) s[, "mu"] <- s[, "mu"] -
    exp(logratio + log(size) - log(mu + size))
  if("size" %in% parameter) s[, "size"] <- s[, "size"] +
    exp(logratio) * (log(size) - log(mu + size) + 1 - size/(mu + size))
  if(drop) drop(s) else s
}
