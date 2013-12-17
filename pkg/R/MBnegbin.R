MBbinomial <- function(link = "logit")
{
  stopifnot(require("mboost"))
  if(!inherits(link, "link-glm")) link <- make.link(link)
  mboost::Family(
    ngradient = function(y, f, w = 1) w * sbinom(y, prob = link$linkinv(f), size = 1L) * link$mu.eta(f),
    loss = function(y, f) -dbinom(y, prob = link$linkinv(f), size = 1L, log = TRUE),
    offset = function(y, w) link$linkfun(weighted.mean(y, w)),
    response = function(f) link$linkinv(f),
    rclass = function(f) (f > 0) + 1 ,
    check_y = function(y) {
      if(is.factor(y)) y <- as.integer(y) - 1L
      stopifnot(all(y >= 0L) & all(y <= 1L))
      y
    },
    name = sprintf("Binomial (%s link)", link$name)
  )
}

MBnegbin <- function(theta = NULL, link = "log",
  control = list(reltol = .Machine$double.eps^(1/1.5), maxit = 500))
{
  stopifnot(require("mboost"))
  
  ## currently used value of theta (possibly fixed)
  if(is.null(theta)) {
    theta <- 1
    fix <- FALSE
  } else {
    fix <- TRUE
  }
  
  ## link function
  if(!inherits(link, "link-glm")) link <- make.link(link)

  ## negative gradient that estimates theta (if necessary)
  ngradient <- function(y, f, w = 1) {
    if(!fix) {
      nll <- function(par) -sum(w * dnbinom(y, mu = link$linkinv(f), size = exp(par), log = TRUE))
      gr <- function(par) -sum(w * snbinom(y, mu = link$linkinv(f), size = exp(par), parameter = "size") * exp(par))
      theta <<- exp(optim(par = log(theta), fn = nll, gr = gr, method = "BFGS", control = control)$par)
    }
    w * snbinom(y, mu = link$linkinv(f), size = theta, parameter = "mu") * link$mu.eta(f)
  }

  mboost::Family(
    ngradient = ngradient,
    loss = function(y, f) -dnbinom(y, mu = link$linkinv(f), size = theta, log = TRUE),
    check_y = function(y) {
      stopifnot(all(y >= 0))
      y
    },
    nuisance = function() return(theta),
    name = "Negative binomial",
    response = function(f) link$linkinv(f)
  )
}

MBztnegbin <- function(theta = NULL, link = "log",
  control = list(reltol = .Machine$double.eps^(1/1.5), maxit = 500))
{
  stopifnot(require("mboost"))

  ## currently used value of theta (possibly fixed)
  if(is.null(theta)) {
    theta <- 1
    fix <- FALSE
  } else {
    fix <- TRUE
  }

  ## link function
  if(!inherits(link, "link-glm")) link <- make.link(link)

  ## negative gradient that estimates theta (if necessary)
  ngradient <- function(y, f, w = 1) {
    if(!fix) {
      nll <- function(par) -sum(w * dztnbinom(y, mu = link$linkinv(f), size = exp(par), log = TRUE))
      gr <- function(par) -sum(w * sztnbinom(y, mu = link$linkinv(f), size = exp(par), parameter = "size") * exp(par))
      theta <<- exp(optim(par = log(theta), fn = nll, gr = gr, method = "BFGS", control = control)$par)
    }
    w * sztnbinom(y, mu = link$linkinv(f), size = theta, parameter = "mu") * link$mu.eta(f)
  }

  mboost::Family(
    ngradient = ngradient,
    loss = function(y, f) -dztnbinom(y, mu = link$linkinv(f), size = theta, log = TRUE),
    check_y = function(y) {
      stopifnot(all(y > 0))
      y
    },
    nuisance = function() return(theta),
    name = "Zero-truncated negative binomial",
    response = function(f) link$linkinv(f)
  )
}
