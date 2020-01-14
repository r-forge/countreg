## ztnbinom: Zero-truncated negative binomial
dztnbinom <- function(x, mu, theta, size, log = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- dnbinom(x, mu = mu, size = theta, log = TRUE) - pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  rval[x < 1] <- -Inf
  rval[mu <= 0] <- 0
  if(log) rval else exp(rval)
}

pztnbinom <- function(q, mu, theta, size, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(pnbinom(q, mu = mu, size = theta, lower.tail = lower.tail, log.p = FALSE) - dnbinom(0, mu = mu, size = theta)) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  rval[q < 1] <- if(lower.tail) -Inf else 0
  if(log.p) rval else exp(rval)
}

qztnbinom <- function(p, mu, theta, size, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  p_orig <- p
  p <- if(log.p) p else log(p)
  p <- p + pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  p <- exp(p) + dnbinom(0, mu = mu, size = theta)
  rval <- qnbinom(p, mu = mu, size = theta, lower.tail = lower.tail, log.p = FALSE)
  if(lower.tail) rval[p_orig < dztnbinom(1, mu = mu, theta = theta, log = log.p)] <- 1
  rval
}

rztnbinom <- function(n, mu, theta, size) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  qztnbinom(runif(n), mu = mu, theta = theta)
}

sztnbinom <- function(x, mu, theta, size, parameter = c("mu", "theta", "size"), drop = TRUE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "theta", "size")))
  s <- snbinom(x, mu = mu, size = theta, parameter = ifelse(parameter == "theta", "size", parameter), drop = FALSE)
  colnames(s) <- parameter
  logratio <- pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  if("mu" %in% parameter) s[, "mu"] <- s[, "mu"] -
    exp(logratio + log(theta) - log(mu + theta))
  if("theta" %in% parameter) s[, "theta"] <- s[, "theta"] +
    exp(logratio) * (log(theta) - log(mu + theta) + 1 - theta/(mu + theta))
  s[(x < 1) | (abs(x - round(x)) > sqrt(.Machine$double.eps)), ] <- 0
  if(drop & NCOL(s) < 2L) drop(s) else s
}

hztnbinom <- function(x, mu, theta, size, parameter = c("mu", "theta"), drop = TRUE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size

  parameter <- unique(ifelse(parameter == "theta.mu", "mu.theta", parameter))
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "theta", "mu.theta")))

  para2 <- ifelse(parameter == "theta", "size", parameter)
  para2 <- ifelse(para2 == "mu.theta", "mu.size", para2)
  h <- hnbinom(x, mu = mu, size = theta, parameter = para2, drop = FALSE)
  colnames(h) <- parameter

  logratio <- pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)

  if("mu" %in% parameter) h[, "mu"] <- h[, "mu"] +
    exp(logratio - 2*log(mu + theta) + log(theta) + log(1 + theta)) +
    exp(2*logratio + 2*log(theta) - 2*log(mu + theta))

  if("theta" %in% parameter) {
    term <- (log(theta) - log(mu + theta) + 1 - theta/(mu + theta))^2
    h[, "theta"] <- h[, "theta"] +
      exp(logratio) * ( 1/theta - 2/(mu + theta) + theta/(mu + theta)^2 + term) +
      exp(2*logratio) * term
  }

  if("mu.theta" %in% parameter) {
    ratio <- theta/(mu + theta)
    term  <- ratio*log(ratio)
    h[, "mu.theta"] <- h[, "mu.theta"] -
      exp(logratio) * (term + mu * (theta + 1)/(mu + theta)^2) -
      exp(2*logratio) * (term + mu * theta / (mu + theta)^2)
  }

  h[(x < 1) | (abs(x - round(x)) > sqrt(.Machine$double.eps)), ] <- 0
  if(drop & NCOL(h) < 2L) drop(h) else h
}

mean_ztnbinom <- function(mu, theta, size, drop = TRUE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size

  if(drop) {
    mu / pnbinom(0, mu = mu, size = theta, lower.tail = FALSE)
  } else {
    cbind("mean" = mu / pnbinom(0, mu = mu, size = theta, lower.tail = FALSE))
  }
}

var_ztnbinom <- function(mu, theta, size, drop = TRUE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size

  mean <- mu / pnbinom(0, mu = mu, size = theta, lower.tail = FALSE)
  if(drop) {
    mean * (1 + mu/theta + mu - mean)
  } else {
    cbind("var" = mean * (1 + mu/theta + mu - mean))
  }
}

