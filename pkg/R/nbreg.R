nbreg <- function(formula, data, subset, na.action, weights, offset,
  theta = NULL, link = "log", control = nbreg.control(...),
  model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## link function
  if(is.character(link)) link <- make.link(link)
  stopifnot(inherits(link, "link-glm"))

  ## distribution
  fix <- !is.null(theta)
  llfun <- function(parms) {
    mu <- as.vector(link$linkinv(X %*% parms[1L:k] + offset))
    if(!fix) theta <- exp(parms[k+1])
    sum(weights * suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE)))
  }
  grad <- function(parms) {
    eta <- as.vector(X %*% parms[1L:k] + offset)
    mu <- link$linkinv(eta)    
    if(!fix) theta <- exp(parms[k+1])
    gr <- weights * countreg:::snbinom(Y, mu = mu, size = theta,
      parameter = c("mu", if(fix) NULL else "size"), drop = FALSE)
    gr <- cbind(
      gr[, 1] * link$mu.eta(eta) * X,
      if(fix) NULL else gr[, 2] * theta
    )
    colSums(gr)
  }

  ## call and formula
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  ## call model.frame()
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract terms, model matrices, response
  mt <- terms(formula, data = data)
  X <- model.matrix(mt, mf)
  Y <- model.response(mf, "numeric")

  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(any(Y < 0)) stop("invalid dependent variable, negative counts")
  if(!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001)))))
    stop("invalid dependent variable, non-integer values")
  Y <- as.integer(round(Y + 0.001))
  
  ## convenience variables
  n <- length(Y)
  k <- NCOL(X)

  ## weights and offset
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- 0
  if(length(offset) == 1) offset <- rep(offset, n)
  offset <- as.vector(offset)

  ## starting values
  start <- control$start
  start_theta <- NULL
  if(!is.null(start)) {
    if(length(start) != k + !fix) {
      warning("invalid starting values, model coefficients not correctly specified")
      start <- NULL
    }
  }
  
  ## default starting values
  if(is.null(start)) start <- c(
    glm.fit(X, Y, family = poisson(), weights = weights, offset = offset)$coefficients,
    if(fix) NULL else 0
  )

  ## model fitting
  ## control parameters
  method <- control$method
  hessian <- control$hessian
  ocontrol <- control
  control$method <- control$hessian <- control$start <- NULL

  ## ML estimation
  fit <- optim(fn = llfun, gr = grad, par = start,
    method = method, hessian = hessian, control = control)

  ## coefficients
  cf <- fit$par[1L:k]
  if(!fix) theta <- as.vector(exp(fit$par[k + 1L]))
  names(cf) <- colnames(X)

  ## covariances
  vc <- -solve(as.matrix(fit$hessian))
  if(fix) {
    SE.logtheta <- NULL
  } else {
    SE.logtheta <- as.vector(sqrt(diag(vc)[k + 1L]))
    vc <- vc[-(k+1L), -(k+1L), drop = FALSE]
  }
  colnames(vc) <- rownames(vc) <- colnames(X)

  ## fitted and residuals
  mu <- link$linkinv(X %*% cf + offset)[,1L]
  res <- sqrt(weights) * (Y - mu)

  ## effective observations
  nobs <- sum(weights > 0)

  rval <- list(coefficients = cf,
    residuals = res,
    fitted.values = mu,
    optim = fit,
    method = method,
    control = control,
    start = start,
    weights = if(identical(as.vector(weights), rep(1, n))) NULL else weights,
    offset = if(identical(offset, rep(0, n))) NULL else offset,
    n = nobs,
    df.null = nobs - 1L - !fix,
    df.residual = nobs - k - !fix,
    terms = mt,
    theta = theta,
    SE.logtheta = SE.logtheta,
    loglik = fit$value,
    vcov = vc,
    fixed = fix,
    link = link,
    converged = fit$convergence < 1,
    call = cl,
    formula = formula,
    levels = .getXlevels(mt, mf),
    contrasts = attr(X, "contrasts")
  )
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- X
      
  class(rval) <- "nbreg"
  return(rval)
}

nbreg.control <- function(method = "BFGS", maxit = 10000, start = NULL, ...) {
  rval <- list(method = method, maxit = maxit, start = start)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- -1
  if(!is.null(rval$hessian)) warning("hessian must not be modified")
  rval$hessian <- TRUE
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.6)
  rval
}

coef.nbreg <- function(object, ...) {
  object$coefficients
}

vcov.nbreg <- function(object, ...) {
  object$vcov
}

logLik.nbreg <- function(object, ...) {
  structure(object$loglik, df = object$n - object$df.residual, nobs = object$n, class = "logLik")
}

nobs.nbreg <- function(object, ...) object$n

print.nbreg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    cat(paste("Coefficients (negative binomial with ", x$link$name, " link):\n", sep = ""))
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    cat(paste(ifelse(x$fixed, "Theta (fixed) =", "Theta ="),
      round(x$theta, digits), "\n\n"))
  }
  
  invisible(x)
}



## ---------------------------



summary.nbreg <- function(object,...)
{
  ## deviance residuals
  object$residuals <- residuals(object, type = "deviance")

  ## compute z statistics
  cf <- object$coefficients
  se <- sqrt(diag(object$vcov))
  k <- length(cf)
  
  if(object$dist == "negbin") {
    cf <- c(cf, "Log(theta)" = as.vector(log(object$theta)))
    se <- c(se, object$SE.logtheta)
  }
  zstat <- cf/se
  pval <- 2*pnorm(-abs(zstat))
  cf <- cbind(cf, se, zstat, pval)
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients <- cf

  ## number of iterations
  object$iterations <- tail(na.omit(object$optim$count), 1)
  
  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.nbreg"
  object
}

print.summary.nbreg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {

    if(x$dist == "negbin-fixed") {
      dist <- "negbin"
      fixed <- TRUE
    } else {
      dist <- x$dist
      fixed <- FALSE
    }

    cat("Deviance residuals:\n")
    print(structure(quantile(x$residuals),
      names = c("Min", "1Q", "Median", "3Q", "Max")), digits = digits, ...)  

    cat(paste("\nCoefficients (truncated ", dist, " with log link):\n", sep = ""))
    printCoefmat(x$coefficients, digits = digits, ...)
  
    if(dist == "negbin") cat(paste(ifelse(fixed, "\nTheta (fixed) =", "\nTheta ="),
      round(x$theta, digits)))
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits), "on", x$n - x$df.residual, "Df\n")
  }
  
  invisible(x)
}

terms.nbreg <- function(x, ...) {
  x$terms
}

model.frame.nbreg <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  NextMethod()
}

model.matrix.nbreg <- function(object, ...) {
  rval <- if(!is.null(object$x)) object$x
    else model.matrix(object$terms, model.frame(object), contrasts = object$contrasts)
  return(rval)
}

predict.nbreg <- function(object, newdata, type = c("response", "prob", "count", "zero"),
  na.action = na.pass, ...)
{
    type <- match.arg(type)

    ## if no new data supplied
    if(missing(newdata)) {
      if(type != "response") {
        if(!is.null(object$x)) {
	  X <- object$x
	} else if(!is.null(object$model)) {
          X <- model.matrix(object$terms, object$model, contrasts = object$contrasts)
	} else {
	  stop("predicted probabilities cannot be computed with missing newdata")
	}
	offset <- if(is.null(object$offset)) rep(0, NROW(X)) else object$offset
      } else {
        return(object$fitted.values)
      }
    } else {
      mf <- model.frame(delete.response(object$terms), newdata, na.action = na.action, xlev = object$levels)
      X <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)
      offset <- if(!is.null(off.num <- attr(object$terms, "offset"))) 
          eval(attr(object$terms, "variables")[[off.num + 1]], newdata)
        else if(!is.null(object$offset)) eval(object$call$offset, newdata)
      if(is.null(offset)) offset <- rep(0, NROW(X))
    }

    mu <- exp(X %*% object$coefficients + offset)[,1]
    p0 <- if(object$dist == "poisson") ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE)
      else pnbinom(0, size = object$theta, mu = mu, lower.tail = FALSE, log.p = TRUE)
    
    if(type == "response") rval <- exp(log(mu) - p0)
    if(type == "count") rval <- mu
    if(type == "zero") rval <- exp(p0)

    ## predicted probabilities
    if(type == "prob") {
      y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
      yUnique <- 1:max(y)
      nUnique <- length(yUnique)
      rval <- matrix(NA, nrow = length(mu), ncol = nUnique)
      dimnames(rval) <- list(rownames(X), yUnique)
      
      if(object$dist == "poisson") {
        for(i in 1:nUnique) rval[,i] <- exp(dpois(yUnique[i], lambda = mu, log = TRUE) - p0)
      } else {
        for(i in 1:nUnique) rval[,i] <- exp(dnbinom(yUnique[i], mu = mu, size = object$theta, log = TRUE) - p0)
      }
    }
   
    rval
}

fitted.nbreg <- function(object, ...) {
  object$fitted.values
}

residuals.nbreg <- function(object, type = c("deviance", "pearson", "response"), ...) {

  type <- match.arg(type)
  res <- object$residuals
  wts <- object$weights
  if(is.null(wts)) wts <- 1
  
  switch(type,
  
  "response" = {
    return(res)
  },
  
  "pearson" = {
    mu <- predict(object, type = "count")
    p0 <- mu/object$fitted.values
    theta1 <- 1/object$theta
    vv <- (mu + (1 + 1/object$theta - 1/p0) * mu^2) / p0
    return(res/sqrt(vv))
  },
  
  "deviance" = {
    yhat <- object$fitted.values
    y <- yhat + sqrt(wts) * object$residuals
    mu <- predict(object, type = "count")
    theta <- object$theta
    
    if(object$dist == "poisson") {
      mu2y <- function(mu) mu / ifelse(mu > 0, ppois(0, lambda = mu, lower.tail = FALSE), 1)
      y2mu <- function(y) {
        yunique <- sort(unique(y))
	munique <- sapply(yunique,
	  function(z) uniroot(function(mu) z - mu2y(mu), interval = c(0, z))$root)
	munique[factor(y, levels = yunique)]
      }      
      ll <- function(mu) {
  	dpois(y, lambda = mu, log = TRUE) -
        ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE)
      }
    } else {
      mu2y <- function(mu) mu / ifelse(mu > 0, pnbinom(0, size = theta, mu = mu, lower.tail = FALSE), 1)
      y2mu <- function(y) {
        yunique <- sort(unique(y))
	munique <- sapply(yunique,
	  function(z) uniroot(function(mu) z - mu2y(mu), interval = c(0, z))$root)
	munique[factor(y, levels = yunique)]
      }
      ll <- function(mu) {
        dnbinom(y, size = theta, mu = mu, log = TRUE) - 
  	pnbinom(0, size = theta, mu = mu, lower.tail = FALSE, log.p = TRUE)
      }
    }
    
    return(sqrt(wts) * sign(y - yhat) * sqrt(2 * abs(ll(y2mu(y)) - ll(mu))))
  })
}

predprob.nbreg <- function(obj, ...){
    predict(obj, type = "prob", ...)
}

extractAIC.nbreg <- function(fit, scale = NULL, k = 2, ...) {
  c(attr(logLik(fit), "df"), AIC(fit, k = k))
}

estfun.nbreg <- function(x, ...) {
  ## extract data
  Y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  X <- model.matrix(x)
  beta <- coef(x)
  theta <- x$theta
  offset <- if(is.null(x$offset)) 0 else x$offset
  wts <- weights(x)
  if(is.null(wts)) wts <- 1

  ## count component: working residuals
  eta <- as.vector(X %*% beta + offset)
  mu <- exp(eta)

  wres <- as.numeric(Y > 0) * switch(x$dist,
    "poisson" = {
      (Y - mu) - exp(ppois(0, lambda = mu, log.p = TRUE) -
        ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE) + eta)    
    },
    "geometric" = {
      (Y - mu * (Y + 1)/(mu + 1)) - exp(pnbinom(0, mu = mu, size = 1, log.p = TRUE) -
        pnbinom(0, mu = mu, size = 1, lower.tail = FALSE, log.p = TRUE) - log(mu + 1) + eta)
    },
    "negbin" = {
      (Y - mu * (Y + theta)/(mu + theta)) - exp(pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
        pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE) +
	log(theta) - log(mu + theta) + eta)
    })

  ## compute gradient from data
  rval <- cbind(wres * wts * X)
  colnames(rval) <- names(beta)
  rownames(rval) <- rownames(X)
  return(rval)
}

getSummary.nbreg <- function(obj, alpha = 0.05, ...) {
  ## extract summary object
  s <- summary(obj)
  
  ## coefficient matrix and confidence interval
  ## compute confidence interval manually to include Log(theta)
  cf <- cbind(s$coefficients,
    s$coefficients[, 1] + qnorm(alpha/2) * s$coefficients[, 2],
    s$coefficients[, 1] + qnorm(1 - alpha/2) * s$coefficients[, 2])
  colnames(cf) <- c("est", "se", "stat", "p", "lwr", "upr")

  ## further summary statistics
  sstat <- c(
    "theta" = s$theta,
    "N" = nobs(obj),
    "logLik" = as.vector(logLik(obj)),
    "AIC" = AIC(obj),
    "BIC" = BIC(obj))

  ## return everything
  return(list(
    coef = cf,
    sumstat = sstat,
    contrasts = obj$contrasts,
    xlevels = obj$levels,
    call = obj$call
  ))
}
