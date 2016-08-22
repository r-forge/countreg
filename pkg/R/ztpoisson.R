.ztpois_lambda_to_mean <- function(lambda) ifelse(lambda <= 0, 1, lambda / (1 - exp(-lambda)))

.make_.ztpois_mean_to_lambda <- function(n = 1000, max = 1000) {
  ## sampling points (equidistant on log-scale)
  x <- exp(seq(log(1/n),  log(max), length.out = n))

  ## set up inverse function
  f <- splinefun(c(1, .ztpois_lambda_to_mean(x)), c(0, x))
  rval <- function(mu) ifelse(mu <= max, f(mu), mu)
}

.ztpois_mean_to_lambda <- .make_.ztpois_mean_to_lambda()

dztpois <- function(x, lambda, mean, log = FALSE) {
  if(!missing(lambda) & !missing(mean)) stop("only 'lambda' or 'mean' may be specified")
  if(!missing(mean)) lambda <- .ztpois_mean_to_lambda(mean)
  rval <- dpois(x, lambda, log = TRUE) - ppois(0, lambda, lower.tail = FALSE, log.p = TRUE)
  rval[x < 1] <- -Inf
  rval[lambda <= 0] <- 0
  if(log) rval else exp(rval)
}

pztpois <- function(q, lambda, mean, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(lambda) & !missing(mean)) stop("only 'lambda' or 'mean' may be specified")
  if(!missing(mean)) lambda <- .ztpois_mean_to_lambda(mean)
  rval <- log(ppois(q, lambda, lower.tail = lower.tail, log.p = FALSE) - dpois(0, lambda)) -
    ppois(0, lambda, lower.tail = FALSE, log.p = TRUE)
  rval[q < 1] <- if(lower.tail) -Inf else 0
  if(log.p) rval else exp(rval)
}

qztpois <- function(p, lambda, mean, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(lambda) & !missing(mean)) stop("only 'lambda' or 'mean' may be specified")
  if(!missing(mean)) lambda <- .ztpois_mean_to_lambda(mean)
  p_orig <- p
  p <- if(log.p) p else log(p)
  p <- p + ppois(0, lambda, lower.tail = FALSE, log.p = TRUE)
  p <- exp(p) + dpois(0, lambda)
  rval <- qpois(p, lambda, lower.tail = lower.tail, log.p = FALSE)
  if(lower.tail) rval[p_orig < dztpois(1, lambda, log = log.p)] <- 1
  rval
}

rztpois <- function(n, lambda, mean) {
  if(!missing(lambda) & !missing(mean)) stop("only 'lambda' or 'mean' may be specified")
  if(!missing(mean)) lambda <- .ztpois_mean_to_lambda(mean)
  qztpois(runif(n), lambda)
}

sztpois <- function(x, lambda, mean, parameter = "lambda", drop = TRUE) {
  if(!missing(lambda) & !missing(mean)) stop("only 'lambda' or 'mean' may be specified")
  if(!missing(mean)) {
    lambda <- .ztpois_mean_to_lambda(mean)
  } else if(parameter == "mean") {
    mean <- .ztpois_lambda_to_mean(lambda)
  }
  parameter <- match.arg(parameter, c("lambda", "mean"))
  s <- x/lambda - 1 - exp(-lambda)/(1 - exp(-lambda))
  if(parameter == "mean") s * lambda / (mean * (lambda + 1 - mean))
  if(drop) s else matrix(s, dimnames = list(names(x), parameter))
}

ztpoisson <- function() {
  ## theta = eta
  ## lambda = exp(eta) = exp(theta)
  ## mean = lambda / (1 - exp(-lambda)) = exp(theta) / (1 - exp(-exp(theta)))

  ## link: "log" but in terms of lambda
  link <- "ztlog"
  stats <- structure(list(
      linkfun = function(mu) log(.ztpois_mean_to_lambda(mu)),
      linkinv = function(eta) .ztpois_lambda_to_mean(exp(eta)),
      mu.eta = function(eta) {
        lambda <- exp(eta)
	mu <- .ztpois_lambda_to_mean(lambda)
	mu * (1 - (mu - lambda))
      },
      valideta = function(eta) TRUE,
      name = "ztlog"
    ), class = "link-glm")

  variance <- function(mu) {
    lambda <- .ztpois_mean_to_lambda(mu)
    mu * (1 + lambda - mu)
  }
  validmu <- function(mu) all(mu > 1)
  dev.resids <- function(y, mu, wt) {
    -2 * wt * (dztpois(y, mean = mu, log = TRUE) - dztpois(y, mean = y, log = TRUE))
  }
  aic <- function(y, n, mu, wt, dev) {
    -2 * sum(dztpois(y, mean = mu, log = TRUE) * wt)
  }
  initialize <- expression({
    if (any(y < 1)) stop("zero or negative values not allowed for the zero-truncated Poisson family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  simfun <- function(object, nsim) {
    wts <- object$prior.weights
    if(any(wts != 1)) warning("ignoring prior weights")
    ftd <- fitted(object)
    rztpois(nsim * length(ftd), mean = ftd)
  }

  .make_ztpois_dlink <- function(n = 1e5, max = 100, order = 2) {
    m <- seq(1 + 1/n, max, by = max/n)
    f <- splinefun(m[-(1:order)], diff(log(.ztpois_mean_to_lambda(m)), differences = order) * (n/max)^order)
    dlog <- function(mu, order = 2) switch(as.character(order),
      "2" = -1/mu^2,
      "3" = 2/mu^3,
      "4" =  -6/mu^4
    )
    function(mu) ifelse(mu <= max, f(mu), dlog(mu, order = order))
  }

  dvar <- function(mu) {
    lambda <- .ztpois_mean_to_lambda(mu)
    V <- mu * (1 + lambda - mu)
    lambda + 1 - 2 * mu + mu * lambda / V
  }
  d2var <- function(mu) {
    lambda <- .ztpois_mean_to_lambda(mu)
    V <- mu * (1 + lambda - mu)
    -mu * lambda/V^2 * (mu * lambda / V + lambda - 2 * mu) + 2 * lambda / V - 2
  }
  .make_ztpois_d3var <- function(n = 1e5, max = 100) {
    m <- seq(1 + 1/n, max, by = max/n)
    f <- splinefun(m[-1], diff(d2var(m)) * n/max)
    function(mu) ifelse(mu <= max, f(mu), 0)
  }
  
  structure(list(
    family = "ztpoisson", ## poisson to get dispersion = 1
    link = link,
    linkfun = stats$linkfun, 
    linkinv = stats$linkinv,
    variance = variance,
    dev.resids = dev.resids, 
    aic = aic,
    mu.eta = stats$mu.eta,
    initialize = initialize, 
    validmu = validmu,
    valideta = stats$valideta,
    simulate = simfun,
    d2link = .make_ztpois_dlink(order = 2, max = 100),
    d3link = .make_ztpois_dlink(order = 3, max = 30),
    d4link = .make_ztpois_dlink(order = 4, max = 10),
    dvar = dvar,
    d2var = d2var, 
    d3var = .make_ztpois_d3var(),
    ls = function(y, w, n, scale) {
      c(sum(dztpois(y, mean = y + 1e-8, log = TRUE) * w), 0, 0)
    }),
    class = "family")
}

