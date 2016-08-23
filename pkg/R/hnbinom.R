## hnbinom: Hurdle Poisson
dhnbinom <- function(x, mu, theta, size, pi, log = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- dnbinom(x, mu = mu, size = theta, log = TRUE) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE) + 
    log(pi)
  rval[x == 0L] <- log(1 - pi)[x == 0L]
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
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "theta", "pi")))
  n <- max(length(x), length(mu), length(theta), length(pi))
  s <- if(any(parameter %in% c("mu", "theta"))) {
    sztnbinom(rep(x, length.out = n), mu = mu, theta = theta, parameter = parameter[parameter != "pi"], drop = FALSE)
  } else {
    NULL
  }
  if("pi" %in% parameter) {
    sp <- rep(1/pi, length.out = n)
    sp[x == 0L] <- -1/(1 - pi)[x == 0L]
    sp[(x < 0) | (abs(x - round(x)) > sqrt(.Machine$double.eps))] <- 0
    s <- cbind(s, "pi" = sp)
  }
  if(drop & (NCOL(s) < 2L)) drop(s) else s
}
