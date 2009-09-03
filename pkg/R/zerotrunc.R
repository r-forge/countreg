## hurdle/zeroinfl todo:
##   print methods
##   deviance residuals
##   weights in residuals
##   object$weights = 0 for rep(1, n)

zerotrunc <- function(formula, data, subset, na.action, weights, offset,
                      dist = c("poisson", "negbin", "geometric"), theta = Inf,
	   	      control = zerotrunc.control(...),
		      model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## set up likelihood components
  llPoisson <- function(parms) {
    mu <- as.vector(exp(X %*% parms + offset))
    sum(weights * (
      dpois(Y, lambda = mu, log = TRUE) - ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE)
    ))
  }

  llNegBin <- function(parms) {
    ## parameters
    mu <- as.vector(exp(X %*% parms[1:k] + offset))
    theta <- exp(parms[k+1])
    sum(weights * suppressWarnings(
      dnbinom(Y, size = theta, mu = mu, log = TRUE) - 
      pnbinom(0, size = theta, mu = mu, lower.tail = FALSE, log.p = TRUE)
    ))
  }

  llGeom <- function(parms) llNegBin(c(parms, 0))

  llNBFixed <- function(parms) llNegBin(c(parms, log(theta)))

  gradPoisson <- function(parms) {
    eta <- as.vector(X %*% parms + offset)
    mu <- exp(eta)
    colSums(((Y - mu) - exp(ppois(0, lambda = mu, log.p = TRUE) -
      ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE) + eta)) * weights * X)
  }
  
  gradGeom <- function(parms) {
    eta <- as.vector(X %*% parms + offset)
    mu <- exp(eta)      
    colSums(((Y - mu * (Y + 1)/(mu + 1)) -
      exp(pnbinom(0, mu = mu, size = 1, log.p = TRUE) -
        pnbinom(0, mu = mu, size = 1, lower.tail = FALSE, log.p = TRUE) -
	log(mu + 1) + eta)) * weights * X)
  }

  gradNegBin <- function(parms) {
    eta <- as.vector(X %*% parms[1:k] + offset)
    mu <- exp(eta)      
    theta <- exp(parms[k+1])
    logratio <- pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
        pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
    rval <- colSums(((Y - mu * (Y + theta)/(mu + theta)) -
      exp(logratio + log(theta) - log(mu + theta) + eta)) * weights * X)
    rval2 <- sum((digamma(Y + theta) - digamma(theta) +    
      log(theta) - log(mu + theta) + 1 - (Y + theta)/(mu + theta) +
      exp(logratio) * (log(theta) - log(mu + theta) + 1 - theta/(mu + theta))) * weights) * theta
    c(rval, rval2)
  }  
  
  gradNBFixed <- function(parms) {
    eta <- as.vector(X %*% parms + offset)
    mu <- exp(eta)      
    logratio <- pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
        pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
    colSums(((Y - mu * (Y + theta)/(mu + theta)) -
      exp(logratio + log(theta) - log(mu + theta) + eta)) * weights * X)
  }

  ## dist/theta processing
  if(!missing(theta)) {
    if(!missing(dist)) warning("supply either 'theta' or 'dist' but not both, 'dist' is ignored")
    if(!is.null(theta) && theta <= 0) stop("theta has to be > 0")
    dist <- if(is.null(theta)) {
      "negbin"
    } else if(!is.finite(theta)) {
      "poisson"
    } else if(isTRUE(all.equal(theta, 1))) {
      "geometric"
    } else {
      "negbin-fixed"
    }
  } else {
    dist <- match.arg(dist)
    theta <- switch(dist,
      "poisson" = Inf,
      "geometric" = 1,
      "negbin" = NULL)
  }

  ## collect likelihood components
  llfun <- switch(dist,
                  "poisson" = llPoisson,
	          "geometric" = llGeom,
	          "negbin" = llNegBin,
		  "negbin-fixed" = llNBFixed)
  grad <- switch(dist,
                 "poisson" = gradPoisson,
		 "geometric" = gradGeom,
		 "negbin" = gradNegBin,
		 "negbin-fixed" = gradNBFixed)

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
  if(any(Y == 0)) stop("invalid dependent variable, no zeros allowed")  
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
    if(dist != "negbin") {
      start_theta <- NULL
      if(length(start) != k) {
        warning("invalid starting values, model coefficients not correctly specified")
        start <- NULL
      }
    } else {
      if(length(start) == k + 1) {
        start_theta <- start[k+1]
        start <- start[1:k]
      } else if(length(start) == k) {
        start_theta <- NULL
      } else {
   	warning("invalid starting values, model coefficients not correctly specified")
   	start <- start_theta <- NULL  
      }
    }
  }
  
  ## default starting values
  if(is.null(start)) start <- glm.fit(X, Y, family = poisson(), weights = weights, offset = offset)$coefficients
  if(is.null(start_theta)) start_theta <- 1
  if(dist == "negbin") start <- c(start, log(start_theta))

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
  cf <- fit$par[1:k]
  if(dist == "negbin") theta <- as.vector(exp(fit$par[k + 1]))
  names(cf) <- colnames(X)

  ## covariances
  vc <- -solve(as.matrix(fit$hessian))
  if(dist == "negbin") {
    SE.logtheta <- as.vector(sqrt(diag(vc)[k + 1]))
    vc <- vc[-(k+1), -(k+1), drop = FALSE]
  } else {
    SE.logtheta <- NULL
  }
  colnames(vc) <- rownames(vc) <- colnames(X)

  ## fitted and residuals
  mu <- exp(X %*% cf + offset)[,1]
  p0 <- if(dist == "poisson") ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE)
    else pnbinom(0, size = theta, mu = mu, lower.tail = FALSE, log.p = TRUE)
  Yhat <- exp(log(mu) - p0)
  res <- sqrt(weights) * (Y - Yhat)

  rval <- list(coefficients = cf,
    residuals = res,
    fitted.values = Yhat,
    optim = fit,
    method = method,
    control = control,
    start = start,
    weights = if(identical(as.vector(weights), rep(1, n))) NULL else weights,
    offset = if(identical(offset, rep(0, n))) NULL else offset,
    n = n,
    df.null = n - 1 - (dist == "negbin"),     ## effective: df.null     - sum(weights == 0)
    df.residual = n - k - (dist == "negbin"), ## effective: df.residual - sum(weights == 0)
    terms = mt,
    theta = theta,
    SE.logtheta = SE.logtheta,
    loglik = fit$value,
    vcov = vc,
    dist = dist,
    converged = fit$convergence < 1,
    call = cl,
    formula = formula,
    levels = .getXlevels(mt, mf),
    contrasts = attr(X, "contrasts")
  )
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- X
      
  class(rval) <- "zerotrunc"
  return(rval)
}

zerotrunc.control <- function(method = "BFGS", maxit = 10000, start = NULL, ...) {
  rval <- list(method = method, maxit = maxit, start = start)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- -1
  if(!is.null(rval$hessian)) warning("hessian must not be modified")
  rval$hessian <- TRUE
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.6)
  rval
}

coef.zerotrunc <- function(object, ...) {
  object$coefficients
}

vcov.zerotrunc <- function(object, ...) {
  object$vcov
}

logLik.zerotrunc <- function(object, ...) {
  structure(object$loglik, df = object$n - object$df.residual, class = "logLik")
}

print.zerotrunc <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  if(x$dist == "negbin-fixed") {
    dist <- "negbin"
    fixed <- TRUE
  } else {
    dist <- x$dist
    fixed <- FALSE
  }

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    cat(paste("Coefficients (truncated ", dist, " with log link):\n", sep = ""))
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    if(dist == "negbin") cat(paste(ifelse(fixed, "Theta (fixed) =", "Theta ="),
      round(x$theta, digits), "\n\n"))
  }
  
  invisible(x)
}

summary.zerotrunc <- function(object,...)
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
  class(object) <- "summary.zerotrunc"
  object
}

print.summary.zerotrunc <- function(x, digits = max(3, getOption("digits") - 3), ...)
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

terms.zerotrunc <- function(x, ...) {
  x$terms
}

model.frame.zerotrunc <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  NextMethod()
}

model.matrix.zerotrunc <- function(object, ...) {
  rval <- if(!is.null(object$x)) object$x
    else model.matrix(object$terms, model.frame(object), contrasts = object$contrasts)
  return(rval)
}

predict.zerotrunc <- function(object, newdata, type = c("response", "prob", "count", "zero"),
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

fitted.zerotrunc <- function(object, ...) {
  object$fitted.values
}

residuals.zerotrunc <- function(object, type = c("deviance", "pearson", "response"), ...) {

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

predprob.zerotrunc <- function(obj, ...){
    predict(obj, type = "prob", ...)
}

extractAIC.zerotrunc <- function(fit, scale = NULL, k = 2, ...) {
  c(attr(logLik(fit), "df"), AIC(fit, k = k))
}

