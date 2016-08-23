## hpois: Hurdle Poisson
dhpois <- function(x, lambda, pi, log = FALSE) {
  rval <- dpois(x, lambda = lambda, log = TRUE) -
    ppois(0, lambda = lambda, lower.tail = FALSE, log.p = TRUE) + 
    log(pi)
  rval[x == 0L] <- log(1 - pi)[x == 0L]
  if(log) rval else exp(rval)
}

phpois <- function(q, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  rval <- log(pmax(0, ppois(q, lambda = lambda) - dpois(0L, lambda = lambda))) -
    ppois(0, lambda = lambda, lower.tail = FALSE, log.p = TRUE) + 
    log(pi)
  rval <- exp(rval)
  rval[q == 0L] <- 0
  rval <- rval + (1 - pi)
  rval[q < 0] <- 0
  if(!lower.tail) rval <- 1 - rval
  if(log.p) log(rval) else rval
}

qhpois <- function(p, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p <- log(pmax(0, p + pi - 1)) - log(pi)
  rval <- qztpois(p, lambda = lambda, lower.tail = TRUE, log.p = TRUE)
  rval[!is.finite(p)] <- 0L
  rval[pi < 0 | pi > 1] <- NaN
  rval
}

rhpois <- function(n, lambda, pi) {
  qhpois(runif(n), lambda = lambda, pi = pi)
}

shpois <- function(x, lambda, pi, parameter = c("lambda", "pi"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("lambda", "pi")))
  n <- max(length(x), length(lambda), length(pi))
  if("pi" %in% parameter) {
    sp <- rep(1/pi, length.out = n)
    sp[x == 0L] <- -1/(1 - pi)[x == 0L]
    sp[(x < 0) | (abs(x - round(x)) > sqrt(.Machine$double.eps))] <- 0
  }
  s <- cbind(
    if("lambda" %in% parameter) sztpois(rep(x, length.out = n), lambda = lambda) else NULL,
    if("pi" %in% parameter) sp else NULL
  )
  colnames(s) <- c("lambda", "pi")[c("lambda", "pi") %in% parameter]
  if(drop) drop(s) else s
}
