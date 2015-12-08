## binom: Binomial
sbinom <- function(x, prob, size, parameter = "prob", drop = TRUE) {
  s <- x/prob - (size - x)/(1 - prob)
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
  if(drop) drop(s) else s
}

## ## ztpoisson: zero-truncated Poisson
## ## -> more flexible implementation in ztpois.R now
## dztpois <- function(x, lambda, log = FALSE) {
##   rval <- dpois(x, lambda = lambda, log = TRUE) -
##     ppois(0, lambda = lambda, lower.tail = FALSE, log.p = TRUE)
##   if(log) rval else exp(rval)
## }
## 
## sztpois <- function(x, lambda, parameter = "lambda", drop = TRUE) {
##   s <- x/lambda - 1 - exp(-lambda)/(1 - exp(-lambda))
##   if(drop) s else cbind("lambda" = s)
## }

## ztnbinom: Zero-truncated negative binomial
dztnbinom <- function(x, mu, size, log = FALSE) {
  rval <- dnbinom(x, mu = mu, size = size, log = TRUE) -
    pnbinom(0, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
  if(log) rval else exp(rval)
}

sztnbinom <- function(x, mu, size, parameter = c("mu", "size"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("mu", "size")))
  s <- snbinom(x, mu = mu, size = size, parameter = parameter, drop = FALSE)
  logratio <- pnbinom(0, mu = mu, size = size, log.p = TRUE) -
    pnbinom(0, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
  if("mu" %in% parameter) s[, "mu"] <- s[, "mu"] -
    exp(logratio + log(size) - log(mu + size))
  if("size" %in% parameter) s[, "size"] <- s[, "size"] +
    exp(logratio) * (log(size) - log(mu + size) + 1 - size/(mu + size))
  if(drop) drop(s) else s
}

## logseries: Log-series
## Mean: prob / ((prob - 1) * log(1 - prob) ## prob/(1 - prob)  * (-1/log(1 - prob)))
## Var: -prob * (prob + log(1 - prob)) / ((prob - 1)^2 * log(1 - prob)^2)
dlogseries <- function(x, prob = 0.5, log = FALSE) { 
  if(any(prob < 0) | any(prob > 1))  stop("'prob' must be in [0, 1]")
  rval <- x * log(prob) - log(abs(x)) - log(-log(1 - prob))
  x <- rep(x, length.out = length(rval))
  rval[x < 1 | x != as.integer(x)] <- -Inf
  if(log) rval else exp(rval)
}

slogseries <- function(x, prob = 0.5, parameter = "prob", drop = TRUE) {
  if(any(prob < 0) | any(prob > 1))  stop("'prob' must be in [0, 1]")
  s <- x/prob + 1/(log(1 - prob) * (1 - prob))
  if(drop) s else cbind("prob" = s)
}

plogseries <- function(q, prob = 0.5, lower.tail = TRUE, log.p = FALSE) {
  if(any(prob < 0) | any(prob > 1))  stop("'prob' must be in [0, 1]")
  uprob <- sort(unique(prob))
  pr <- lapply(uprob, function(p) cumsum(dlogseries(0L:ceiling(max(c(1, q))), prob = p)))
  pr <- do.call("cbind", pr)
  pr <- pr[cbind(pmax(floor(q), 0L) + 1L, match(prob, uprob))]
  if(!lower.tail) pr <- 1 - pr
  if(log.p) log(pr) else pr
}

qlogseries <- function(p, prob = 0.5, lower.tail = TRUE, log.p = FALSE, round = TRUE) {
  if(any(prob < 0) | any(prob > 1))  stop("'prob' must be in [0, 1]")

  ## probabilities
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p	
  nan <- p < 0 | p > 1

  ## cumulative probabilities for varying prob
  uprob <- sort(unique(prob))
  m <- 999L
  pr <- do.call("cbind", lapply(uprob, function(p)
    cumsum(dlogseries(0L:m, prob = p))))
  while(any(p > max(pr[m + 1L, ]) & p < 1)) {
    pr <- rbind(pr, do.call("cbind", lapply(uprob, function(p)
      cumsum(dlogseries((m + 1L):(m + 1000L), prob = p)))))
    m <- m + 1000L
  }
  q <- matrix(NA, nrow = length(p), ncol = length(uprob))
  for(i in 1L:ncol(q)) q[,i] <- approx(pr[,i], 0L:m, xout = p)$y
  q <- q[cbind(seq_along(p), match(prob, uprob))]
  if(round) q <- as.integer(round(q))
  q[p >= 1] <- Inf
  if(any(nan)) {
    q[nan] <- NaN
    warning("NaNs produced")
  }
  q
} 

rlogseries <- function(n, prob = 0.5) {
  if(any(prob < 0) | any(prob > 1))  stop("'prob' must be in [0, 1]")
  as.integer(ceiling(qlogseries(runif(n), prob = prob, round = FALSE)))
}
