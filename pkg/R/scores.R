## binom: Binomial
sbinom <- function(x, prob, size, parameter = "prob", drop = TRUE) {
  s <- x/prob - (size - x)/(1 - prob)
  s[(x < 0) | (x > size) | (abs(x - round(x)) > sqrt(.Machine$double.eps))] <- 0
  if(drop) s else cbind("prob" = s)
}

hbinom <- function(x, prob, size, parameter = "prob", drop = TRUE) {
  h <- - x/prob^2 - (size - x)/(1 - prob)^2
  h[(x < 0) | (x > size) | (abs(x - round(x)) > sqrt(.Machine$double.eps))] <- 0
  if(drop) h else cbind("prob" = h)
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

hnbinom <- function(x, mu, size, parameter = c("mu", "size"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "size",
                      "mu.size", "size.mu")))
  h <- cbind(
    if("mu" %in% parameter) - x/mu^2 + (x + size)/(mu + size)^2 else NULL,
    if("size" %in% parameter) trigamma(x + size) - trigamma(size) +
      1/size - 2/(mu + size) + (x + size) / (mu + size)^2 else NULL,
    if(any(c("mu.size", "size.mu") %in% parameter)) (x - mu)/(mu + size)^2
  )
  colnames(h) <- c(if("mu" %in% parameter) "mu",
                   if("size" %in% parameter) "size",
                   if(any(c("mu.size", "size.mu") %in% parameter)) "mu.size")
  h[(x < 0) | (abs(x - round(x)) > sqrt(.Machine$double.eps)), ] <- 0
  if(drop & NCOL(h) < 2L) drop(h) else h
}



## pois: Poisson
spois <- function(x, lambda, parameter = "lambda", drop = TRUE) {
  s <- x/lambda - 1
  s[(x < 0) | (abs(x - round(x)) > sqrt(.Machine$double.eps))] <- 0
  if(drop) s else cbind("lambda" = s)
}

hpois <- function(x, lambda, parameter = "lambda", drop = TRUE) {
  h <- - x/lambda^2
  h[(x < 0) | (abs(x - round(x)) > sqrt(.Machine$double.eps))] <- 0
  if(drop) h else cbind("lambda" = h)
}
