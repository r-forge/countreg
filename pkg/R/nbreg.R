nbreg <- function(formula, data, subset, na.action, weights, offset,
  theta = NULL, dist = "NB2", link = "log", link.theta = "log", control = nbreg.control(...),
  model = TRUE, y = TRUE, x = FALSE, z = FALSE, hessA = TRUE,...) # hessA argument to use analytical hessian
{
  ## link function
  if(is.character(link)) {
    linkobj <- make.link2(link)
    stopifnot(inherits(linkobj, "link-glm"))
    link <- linkobj
  }

  ## link function for theta
  if(is.character(link.theta)) {
    linkobj <- make.link2(link.theta)
    stopifnot(inherits(linkobj, "link-glm"))
    link.theta <- linkobj
  }
  
  ## allow lowercase dist
  dist <- toupper(dist)

  switch(dist,
         NB2 = {
           densFun <- dnbinom
           derivFun <- snbinom
           hessFun <- hnbinom
           varFun <- var_nbinom
           pFun <- pnbinom
           dev.resids <- function(y, mu, theta, wt) # inspired by negbin.R in MASS
             2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) *
                         log((y + theta)/ (mu + theta)))
           # why pmax(1, y)? if y = 0 --> then choose 1
           },
         NB1 = {
           densFun <- dnbinom1
           derivFun <- snbinom1
           hessFun <- hnbinom1
           varFun <- var_nbinom1
           pFun <- pnbinom1
           dev.resids <- function(y, mu, theta, wt){NULL} # does not work, defined to avoid errors
         }, stop("invalid model"))

  
  ## distribution
  fix <- !is.null(theta)

  llfun <- function(parms) {
    mu <- link$linkinv(linPred(X, parms[idxMu], offsetx))
    if(!fix) theta <- link.theta$linkinv(linPred(Z, parms[idxEta], offsetz))
    sum(weights * suppressWarnings(densFun(Y, size = theta, mu = mu, log = TRUE)))
  }
  
  grad <- function(parms) {
    eta <- linPred(X, parms[idxMu], offsetx)
    etaTheta <- linPred(Z, parms[idxEta], offsetz)
    mu <- link$linkinv(eta)
    if(!fix) theta <- link.theta$linkinv(etaTheta) # modified to include covariates
    gr <- weights * derivFun(Y, mu = mu, size = theta,
      parameter = c("mu", if(fix) NULL else "size"), drop = FALSE)
    gr <- cbind(
      gr[, 1] * link$mu.eta(eta) * X,
      if(fix) NULL else gr[, 2] * link.theta$mu.eta(etaTheta) * Z
    )
    colSums(gr)
  }
  
  hess <- function(parms) {
    eta <- linPred(X, parms[idxMu], offsetx)
    etaTheta <- linPred(Z, parms[idxEta], offsetz)
    mu <- link$linkinv(eta)
    if(!fix) theta <- link.theta$linkinv(etaTheta)

    gr <- derivFun(Y, mu = mu, size = theta,
                  parameter = c("mu", if(fix) NULL else "size"), drop = FALSE)
    hm <- hessFun(Y, mu = mu, size = theta,
                            parameter = c("mu", if(fix) NULL else "size",
                                          if(fix) NULL else "mu.size"),
                            drop = FALSE)

    if(fix){
      return(hessMat1(hm[,1], gr[,1], X, eta = eta, dlinkX = link$mu.eta, ddlinkX = link$dmu.eta))
    } else {
      return(hessMat2(hm, gr, X, Z, etaX = eta, etaZ = etaTheta, fix = fix,
                      dlinkX = link$mu.eta, ddlinkX = link$dmu.eta,
                      dlinkZ = link.theta$mu.eta, ddlinkZ = link.theta$dmu.eta,
                      weights = weights))
    }
  }

  
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE

  ## formula
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1) # specify theta as constant if no covariates specified
  } else {
    if(!is.null(theta)) stop("invalid formula, theta is already prespecified")
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      stop("formula must not have more than two RHS parts")
    }
  }
  mf$formula <- formula

  ## call model.frame()
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrices, response
  mt <- terms(formula, data = data, dot = control$dot)
  mtX <- terms(formula, data = data, rhs = 1L, dot = control$dot)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L, dot = control$dot))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)

  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(any(Y < 0)) stop("invalid dependent variable, negative counts")
  if(!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001)))))
    stop("invalid dependent variable, non-integer values")
  Y <- as.integer(round(Y + 0.001))

  ## convenience variables
  n <- length(Y)
  k <- NCOL(X)
  if(is.null(Z)) {
    q <- 1L
    Z <- matrix(1, ncol = q, nrow = n)
    colnames(Z) <- "(Intercept)"
    rownames(Z) <- rownames(x)
  } else {
    q <- NCOL(Z)
    if(q < 1L) stop("theta regression needs to have at least one parameter")
  }
  idxMu <- 1L:k
  idxEta <- (k+1L):(k+q)

  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)

  ## offsets
  offsetx <- model_offset_2(mf, terms = mtX, offset = TRUE) # model_offset_2 is defined in hurdle.R
  if(is.null(offsetx)) offsetx <- 0
  if(length(offsetx) == 1) offsetx <- rep.int(offsetx, n)
  offsetx <- as.vector(offsetx)
  offsetz <- model_offset_2(mf, terms = mtZ, offset = FALSE)
  if(is.null(offsetz)) offsetz <- 0
  if(length(offsetz) == 1) offsetz <- rep.int(offsetz, n)
  offsetz <- as.vector(offsetz)

  ## starting values
  start <- c(control$start$mu, control$start$theta)
  if(!is.null(start)) {
    if(length(start) != k + (!fix)*q) {
      warning("invalid starting values, model coefficients not correctly specified")
      start <- NULL
    }
  }

  ## default starting values
  if(is.null(start)) start <- c(
    glm.fit(X, Y, family = poisson(), weights = weights, offset = offsetx)$coefficients,
    if(fix) NULL else rep(0, ncol(Z))
    ## TODO: Fix starting values, currently all coefs of theta param set to 0
    ## TODO: If Z is a constant, use starting value provided in Cameron & Trivedi 2013, p.76, eq. (3.20)
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
  coefTheta <- if(fix) link.theta$linkfun(theta) else as.vector(fit$par[(k + 1L):(k + q)])
  names(cf) <- colnames(X)
  names(coefTheta) <- colnames(Z)

  ## covariances
  if(hessA) {
    vc <- -solve(hess(fit$par))
    hMat <- hess(fit$par)
  } else{
    if (hessian) {
      vc <- -solve(as.matrix(fit$hessian))
      hMat <- fit$hessian
    } else {
      vc <- hMat <- matrix(NA_real_, nrow = k + q, ncol = k + q)
    }

  }
  if(fix) {
    colnames(vc) <- rownames(vc) <- colnames(X)
  } else {
    colnames(vc) <- rownames(vc) <- c(paste("mu", colnames(X), sep = "_"),
                                      paste("theta", colnames(Z), sep = "_"))

  }

  ## fitted and residuals
  mu <- link$linkinv(X %*% cf + offsetx)[,1L]
  res <- sqrt(weights) * (Y - mu)

  ## effective observations
  nobs <- sum(weights > 0)

  rval <- list(coefficients = cf,
    coefficients.theta = coefTheta,
    residuals = res,
    fitted.values = mu,
    optim = fit,
    method = method,
    control = control,
    start = start,
    weights = if(identical(as.vector(weights), rep(1, n))) NULL else weights,
    offset = list(mu = if(identical(offsetx, rep.int(0, n))) NULL else offsetx,
                  theta = if(identical(offsetz, rep.int(0, n))) NULL else offsetz),
    n = nobs,
    df.null = nobs - 1L - !fix,
    df.residual = nobs - k - (!fix)*q,
    terms = list(mu = mtX, theta = mtZ, full = mt), # modified to include Z
    vcov = vc,
    fixed = fix,
    link = link,
    link.theta = link.theta,
    loglik = fit$value,
    converged = fit$convergence < 1,
    call = cl,
    formula = formula,
    levels = .getXlevels(mt, mf),
    contrasts = list(mu = attr(X, "contrasts"), theta = attr(Z, "contrasts")),
    dist = dist,
    densFun = densFun,
    derivFun = derivFun,
    hessFun = hessFun,
    varFun = varFun,
    pFun = pFun,
    dev.resids = dev.resids,
    hMat = hMat

  )
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- X
  if(z) rval$z <- Z

  class(rval) <- "nbreg"
  return(rval)
}

nbreg.control <- function(method = "BFGS", maxit = 10000, start = NULL,
                          hessian = TRUE, dot = "separate", ...) {
  rval <- list(method = method, maxit = maxit, start = start, hessian = hessian)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- -1 # tells optim to maximize
  # if(!is.null(rval$hessian)) warning("hessian must not be modified")
  # rval$hessian <- TRUE
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.6)
  rval
}

coef.nbreg <- function(object, model = c("full", "mu", "theta"), ...) {
  model <- match.arg(model)
  rval <- switch(model,
                 "full" = structure(c(object$coefficients, object$coefficients.theta),
                                    .Names = c(paste("mu", names(object$coefficients), sep = "_"),
                                               paste("theta", names(object$coefficients.theta), sep = "_"))),
                 "mu" = object$coefficients,
                 "theta" = object$coefficients.theta)
  rval
}


vcov.nbreg <- function(object, model = c("full", "mu", "theta"), ...) {
  model <- match.arg(model)
  rval <- object$vcov
  if(model == "full") return(rval)
  
  cf <- if(model == "mu") object$coefficients else object$coefficients.theta
  wi <- seq_along(object$coefficients)
  rval <- if(model == "mu") rval[wi, wi, drop = FALSE] else rval[-wi, -wi, drop = FALSE]
  colnames(rval) <- rownames(rval) <- names(cf)
  return(rval)
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
    cat(paste("Theta Coefficients with ", x$link.theta$name, " link:\n", sep = ""))
    print.default(format(x$coefficients.theta, digits = digits), print.gap = 2, quote = FALSE)
  }

  invisible(x)
}

summary.nbreg <- function(object,...)
{
  ## pearson residuals
  object$residuals <- residuals(object, type = "pearson")

  ## compute z statistics
  kc <- length(object$coefficients)
  kz <- length(object$coefficients.theta)
  se <- sqrt(diag(object$vcov))
  coef <- c(object$coefficients, if(object$fixed) NULL else object$coefficients.theta)

  zstat <- coef/se
  pval <- 2*pnorm(-abs(zstat))
  coef <- cbind(coef, se, zstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients <- coef[1:kc,,drop = FALSE]
  object$coefficients.theta <- (if(object$fixed){
    outMat <- as.matrix(object$coefficients.theta) # convert to matrix so that printCoefMat works in print.summary.nbreg
    colnames(outMat) <- "Fixed value"
    outMat
  } else coef[(kc+1):(kc+kz),,drop = FALSE])

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

    cat("Pearson residuals:\n") # modified to pearson
    print(structure(quantile(x$residuals),
      names = c("Min", "1Q", "Median", "3Q", "Max")), digits = digits, ...)

    cat(paste("\nCoefficients (", dist, " with ", x$link$name ," link):\n", sep = ""))
    printCoefmat(x$coefficients, digits = digits, ...)

    cat(paste(ifelse(fixed, "\nTheta (fixed) =\n",
                     paste("\nTheta coefficients ", "(", x$link.theta$name, " link):\n", sep = ""))))
    printCoefmat(x$coefficients.theta, digits = digits, ...)
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits), "on", x$n - x$df.residual, "Df\n")
  }

  invisible(x)
}

terms.nbreg <- function(x, model = c("full", "mu", "theta"), ...) {
  x$terms[[match.arg(model)]]
}

model.frame.nbreg <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  NextMethod()
}

model.matrix.nbreg <- function(object, model = c("mu", "theta"), ...) {
  model <- match.arg(model)
  if(!is.null(object$x)) rval <- object$x[[model]]
  else if(!is.null(object$model)) rval <- model.matrix(object$terms[[model]],
                                                       object$model,
                                                       contrasts = object$contrasts[[model]])
  else stop("not enough information in fitted model to return model.matrix")
  return(rval)
}

# inspired by predict.zeroinfl
predict.nbreg <- function(object, newdata, type = c("response", "prob", "theta", "parameters"),
                          na.action = na.pass, ...)
{
  type <- match.arg(type)
  
  ## if no new data supplied
  if(missing(newdata)) {
    if(type != "response") {
      if(!is.null(object$x)) {
        X <- object$x
      } else if(!is.null(object$model)) {
        X <- model.matrix(object$terms$mu, object$model,
                          contrasts = object$contrasts$mu)
      } else {
        stop("predicted probabilities cannot be computed with missing newdata")
      }
      if(!is.null(object$z)) {
        Z <- object$z
      } else if(!is.null(object$model)) {
        Z <- model.matrix(delete.response(object$terms$theta), object$model,
                          contrasts = object$contrasts$theta)
      }
      offsetx <- if(is.null(object$offset$mu)) rep.int(0, NROW(X)) else object$offset$mu
      offsetz <- if(is.null(object$offset$theta))  rep.int(0, NROW(Z)) else object$offset$theta
    } else {
      return(object$fitted.values)
    }
  } else {
    mf <- model.frame(delete.response(object$terms$full), newdata, na.action = na.action, xlev = object$levels)
    X <- model.matrix(delete.response(object$terms$mu), mf, contrasts = object$contrasts$mu)
    Z <- model.matrix(delete.response(object$terms$theta), mf, contrasts = object$contrasts$theta)
    offsetx <- model_offset_2(mf, terms = object$terms$mu, offset = FALSE)
    offsetz <- model_offset_2(mf, terms = object$terms$theta,  offset = FALSE)
    if(is.null(offsetx)) offsetx <- rep.int(0, NROW(X))
    if(is.null(offsetz)) offsetz <- rep.int(0, NROW(Z))
    if(!is.null(object$call$offset)) offsetx <- offsetx + eval(object$call$offset, newdata)
    # overwrite if we have custom offset
    if(!is.null(object$offset)){eval(object$call$offset, newdata)}
    else{
      if(!is.null(off.num <- attr(object$terms$mu, "offset"))) {
        for(j in off.num){
          offsetx <- offsetx +
            eval(attr(object$terms$mu, "variables")[[j + 1]], newdata)
        }
      }
      if(!is.null(off.num <- attr(object$terms$theta, "offset"))){
        for(j in off.num){
          offsetz <- offsetz +
            eval(attr(object$terms$theta, "variables"))[[j + 1]]
        }
      }
    }
  }
  
  mu <- object$link$linkinv(linPred(X, object$coefficients, offsetx))
  theta <- object$link.theta$linkinv(linPred(Z, object$coefficients.theta, offsetz))
  
  if(type == "response") rval <- mu
  if(type == "theta") rval <- theta
  if(type == "parameters") rval <- data.frame(mu = mu, theta = theta)
  
  ## predicted probabilities
  if(type == "prob") {
    y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
    yUnique <- 0:max(y)
    nUnique <- length(yUnique)
    rval <- matrix(NA, nrow = length(mu), ncol = nUnique)
    dimnames(rval) <- list(rownames(X), yUnique)
    
    for(i in 1:nUnique) rval[,i] <- exp(object$densFun(yUnique[i], mu = mu, size = theta, log = TRUE))
  }
  
  rval
}

fitted.nbreg <- function(object, ...) {
  object$fitted.values
}

residuals.nbreg <- function(object, type = c("pearson","deviance", "response"), ...) {
  
  type <- match.arg(type)
  res <- object$residuals
  wts <- object$weights
  if(is.null(wts)) wts <- 1
  
  switch(type,
         "response" = {
           return(res)
         },
         "pearson" = {
           mu <- predict(object, type = "response")
           theta <- predict(object, type = "theta")
           vv <- object$varFun(mu = mu, size = theta)
           return(res/sqrt(vv))
         },
         "deviance" = {
           if(object$dist == "NB1") stop("Deviance residuals not supported for negative binomial type 1 (NB1)")
           ll <- function(y, mu, theta){
             suppressWarnings(object$densFun(y, size = theta, mu = mu, log = TRUE))
           }
           yhat <- object$fitted.values
           y <- yhat + object$residuals / sqrt(wts) # recover y from residuals and fitted values
           mu <- predict(object, type = "response")
           theta <- predict(object, type = "theta")
           return(sign(y - yhat) * sqrt(pmax(object$dev.resids(y = y, mu = mu, theta = theta, wt = wts), 0)))
           # why pmax(devresids, 0)? Can get dev >= 0 by mu = muhat if dev < 0.
         })
}

prodist.nbreg <- function(object, ...) {
  stopifnot(requireNamespace("distributions3"))
  p <- predict(object, type = "parameters", ...)
  if(object$dist == "NB1") p$theta <- p$theta * p$mu
  distributions3::NegativeBinomial(mu = p$mu, size = p$theta)
}

extractAIC.nbreg <- function(fit, scale = NULL, k = 2, ...) {
  c(attr(logLik(fit), "df"), AIC(fit, k = k))
}

getSummary.nbreg <- function(obj, alpha = 0.05, ...) {
  ## extract coefficient summary
  cf <- list(summary(obj)$coefficients, summary(obj)$coefficients.theta)
  ## augment with confidence intervals
  cval <- qnorm(1 - alpha/2)
  for(i in seq_along(cf)) cf[[i]] <- cbind(cf[[i]],
                                           cf[[i]][, 1] - cval * cf[[i]][, 2],
                                           cf[[i]][, 1] + cval * cf[[i]][, 2])
  ## collect in array
  nam <- unique(unlist(lapply(cf, rownames)))
  acf <- array(dim = c(length(nam), 6, length(cf)),
               dimnames = list(nam, c("est", "se", "stat", "p", "lwr", "upr"), names(cf)))
  for(i in seq_along(cf)) acf[rownames(cf[[i]]), , i] <- cf[[i]]
  
  ## further summary statistics
  sstat <- c(
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

