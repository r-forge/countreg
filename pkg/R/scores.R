## binom: Binomial
sbinom <- function(x, prob, size, parameter = "prob", drop = TRUE) {
  s <- x/prob - (size - x)/(1 - prob)
  s[(x < 0) | (x > size) | (abs(x - round(x)) > sqrt(.Machine$double.eps))] <- 0
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
  s[(x < 0) | (abs(x - round(x)) > sqrt(.Machine$double.eps)), ] <- 0
  if(drop & NCOL(s) < 2L) drop(s) else s
}

## pois: Poisson
spois <- function(x, lambda, parameter = "lambda", drop = TRUE) {
  s <- x/lambda - 1
  s[(x < 0) | (abs(x - round(x)) > sqrt(.Machine$double.eps))] <- 0
  if(drop) s else cbind("lambda" = s)
}
