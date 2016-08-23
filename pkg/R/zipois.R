## zipois: Zero-inflated Poisson
dzipois <- function(x, lambda, pi, log = FALSE) {
  rval <- log(1 - pi) + dpois(x, lambda = lambda, log = TRUE)
  if(any(x0 <- (x == 0L))) rval[x0] <- log(exp(rval[x0]) + pi)
  rval[pi < 0 | pi > 1] <- NaN
  if(log) rval else exp(rval)
}

pzipois <- function(q, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  rval <- log(1 - pi) + ppois(q, lambda = lambda, lower.tail = lower.tail, log.p = TRUE)
  if(any(q0 <- (is.finite(rval) & (lower.tail | q < 0)))) rval[q0] <- log(exp(rval[q0]) + pi)
  rval[pi < 0 | pi > 1] <- NaN
  if(log.p) rval else exp(rval)
}

qzipois <- function(p, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p <- pmax(0, (p - pi)/(1 - pi))
  rval <- qpois(p, lambda = lambda, lower.tail = TRUE, log.p = FALSE)
  rval[pi < 0 | pi > 1] <- NaN
  rval
}

rzipois <- function(n, lambda, pi) {
  rval <- rpois(n, lambda = lambda)
  rval[runif(n) < pi] <- 0
  rval
}

szipois <- function(x, lambda, pi, parameter = c("lambda", "pi"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("lambda", "pi")))
  s <- 0

  ## FIXME ##

  if(drop) drop(s) else s
}

