## zinbinom: Zero-inflated negative binomial
dzinbinom <- function(x, mu, theta, size, pi, log = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(1 - pi) + dnbinom(x, mu = mu, size = theta, log = TRUE)
  if(any(x0 <- (x == 0L))) rval[x0] <- log(exp(rval[x0]) + pi)
  rval[pi < 0 | pi > 1] <- NaN
  if(log) rval else exp(rval)
}

pzinbinom <- function(q, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(1 - pi) + pnbinom(q, mu = mu, size = theta, lower.tail = lower.tail, log.p = TRUE)
  if(any(q0 <- (is.finite(rval) & (lower.tail | q < 0)))) rval[q0] <- log(exp(rval[q0]) + pi)
  rval[pi < 0 | pi > 1] <- NaN
  if(log.p) rval else exp(rval)
}

qzinbinom <- function(p, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p <- pmax(0, (p - pi)/(1 - pi))
  rval <- qnbinom(p, mu = mu, size = theta, lower.tail = TRUE, log.p = FALSE)
  rval[pi < 0 | pi > 1] <- NaN
  rval
}

rzinbinom <- function(n, mu, theta, size, pi) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- rnbinom(n, mu = mu, size = theta)
  rval[runif(n) < pi] <- 0
  rval
}

szinbinom <- function(x, mu, theta, size, pi, parameter = c("mu", "theta", "pi"), drop = TRUE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "theta", "pi")))
  clp0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
  p0 <- pi * (x < 1) + exp(log(1 - pi) + clp0)
  x1 <- x > 0L
  if("mu" %in% parameter) {
    sm <- -exp(-log(p0) + log(1 - pi) + clp0 + log(theta) - log(mu + theta))
    sm[x1] <- (x/mu - (x + theta)/(mu + theta))[x1]
  }
  if("theta" %in% parameter) {
    st <- exp(-log(p0) + log(1 - pi) + clp0) * (log(theta) - log(mu + theta) + 1 - theta/(mu + theta))
    st[x1] <- (digamma(x + theta) - digamma(theta) + log(theta) - log(mu + theta) + 1 - (x + theta)/(mu + theta))[x1]
  }
  if("pi" %in% parameter) {
    sp <- (1 - exp(clp0))/p0
    sp[x1] <- -1/(1 - pi)
  }
  s <- cbind(
    if("mu" %in% parameter) sm else NULL,
    if("theta" %in% parameter) st else NULL,
    if("pi" %in% parameter) sp else NULL
  )
  colnames(s) <- c("mu", "theta", "pi")[c("mu", "theta", "pi") %in% parameter]
  s[(x < 0) | (abs(x - round(x)) > sqrt(.Machine$double.eps)), ] <- 0
  if(drop & NCOL(s) < 2L) drop(s) else s
}
