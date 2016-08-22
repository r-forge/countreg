## ztnbinom: Zero-truncated negative binomial
dztnbinom <- function(x, mu, theta, size, log = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- dnbinom(x, mu = mu, size = theta, log = TRUE) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  if(log) rval else exp(rval)
}

sztnbinom <- function(x, mu, theta, size, parameter = c("mu", "theta"), drop = TRUE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "theta")))
  s <- snbinom(x, mu = mu, size = theta, parameter = parameter, drop = FALSE)
  logratio <- pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  if("mu" %in% parameter) s[, "mu"] <- s[, "mu"] -
    exp(logratio + log(theta) - log(mu + theta))
  if("theta" %in% parameter) s[, "theta"] <- s[, "theta"] +
    exp(logratio) * (log(theta) - log(mu + theta) + 1 - theta/(mu + theta))
  if(drop) drop(s) else s
}

