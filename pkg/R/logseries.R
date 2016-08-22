## Log-series distribution (d/p/q/r)

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
