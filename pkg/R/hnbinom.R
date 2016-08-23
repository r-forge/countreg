## hnbinom: Hurdle Poisson
dhnbinom <- function(x, mu, theta, size, pi, log = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- dnbinom(x, mu = mu, size = theta, log = TRUE) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE) + 
    log(pi)
  rval[x == 0L] <- log(1 - pi)
  if(log) rval else exp(rval)
}

phnbinom <- function(q, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(pmax(0, pnbinom(q, mu = mu, size = theta) - dnbinom(0L, mu = mu, size = theta))) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE) + 
    log(pi)
  rval <- exp(rval)
  rval[q == 0L] <- 0
  rval <- rval + (1 - pi)
  rval[q < 0] <- 0
  if(!lower.tail) rval <- 1 - rval
  if(log.p) log(rval) else rval
}

qhnbinom <- function(p, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p <- log(pmax(0, p + pi - 1)) - log(pi)
  rval <- qztnbinom(p, mu = mu, theta = theta, lower.tail = TRUE, log.p = TRUE)
  rval[!is.finite(p)] <- 0L
  rval[pi < 0 | pi > 1] <- NaN
  rval
}

rhnbinom <- function(n, mu, theta, size, pi) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  qhnbinom(runif(n), mu = mu, theta = theta, pi = pi)
}

shnbinom <- function(x, mu, theta, size, pi, parameter = c("mu", "theta", "pi"), drop = TRUE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "pi")))
  s <- cbind(
    if("mu" %in% parameter) x/mu - 1 - exp(-mu)/(1 - exp(-mu)) else NULL,
    if("pi" %in% parameter) ifelse(x == 0L, -1, 0) else NULL
  )
  colnames(s) <- c("mu", "theta", "pi")[c("mu", "theta", "pi") %in% parameter]
  if(drop) drop(s) else s
}
