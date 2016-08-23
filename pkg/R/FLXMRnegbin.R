FLXMRnegbin <- function(formula = . ~ ., theta = NULL, offset = NULL,
  control = list(reltol = .Machine$double.eps^(1/1.5), maxit = 500))
{
  .theta <- theta
  
  nbrefit <- function(x, y, w) {
    fit <- c(glm.fit(x, y, weights = w, offset = offset, family = MASS::negative.binomial(theta)),
	     list(call = sys.call(), offset = offset,
		  control = eval(formals(glm.fit)$control),
		  method = "weighted.glm.fit"))
    fit$df.null <- sum(w) + fit$df.null - fit$df.residual - fit$rank - is.null(.theta)
    fit$df.residual <- sum(w) - fit$rank - is.null(.theta)
    fit$x <- x
    fit
  }
	      
  z <- methods::new("FLXMRglm", weighted = TRUE, formula = formula,
	   name = "FLXMRglm: negative.binomial", offset = offset,
	   family = "negative.binomial", refit = nbrefit)

  z@preproc.y <- function(x){
    if (ncol(x) > 1L)
      stop(paste("for the", family, "family y must be univariate"))
    x
  }

  z@defineComponent <- if(is.null(.theta)) {
  
    expression({
      predict <- function(x, ...) {
        dotarg <- list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset
        p <- x %*% coef
        if (!is.null(offset)) p <- p + offset
        exp(p)
      }
      
      logLik <- function(x, y, ...)
        suppressWarnings(dnbinom(y, mu = predict(x, ...), size = theta, log = TRUE))
    
      methods::new("FLXcomponent",
        parameters = list(coef = coef, theta = theta),
        logLik = logLik,
	predict = predict,
        df = df)
    })
    
  } else {
  
    as.expression(substitute({
      predict <- function(x, ...) {
        dotarg <- list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset
        p <- x %*% coef
        if (!is.null(offset)) p <- p + offset
        exp(p)
      }
      
      logLik <- function(x, y, ...)
        suppressWarnings(dnbinom(y, mu = predict(x, ...), size = theta, log = TRUE))
    
      methods::new("FLXcomponent",
        parameters = list(coef = coef),
        logLik = logLik,
	predict = predict,
        df = df)
    }, as.environment(list(theta = .theta))))
    
  }
      
  z@fit <- function(x, y, w, component){
    if(is.null(component$theta)) {
      df <- ncol(x)
      theta <- if(is.null(.theta)) 1 else .theta
      cf <- glm.fit(x, y, weights = w, family = MASS::negative.binomial(theta),
        offset = offset, start = component$coef)$coefficients
    } else {
      ## degrees of freedom
      df <- ncol(x) + 1
      
      ## offset
      if(is.null(offset)) offset <- 0
      
      ## objective function
      nll <- function(par) {
        beta <- par[-df]
        theta <- exp(par[df])
	mu <- exp(drop(x %*% beta + offset))
	suppressWarnings(-sum(w * dnbinom(y, mu = mu, size = theta, log = TRUE)))
      }

      ## corresponding gradients
      gr <- function(par) {
        beta <- par[-df]
        theta <- exp(par[df])
        mu <- exp(drop(x %*% beta + offset))
	gr <- drop(y - mu * (y + theta)/(mu + theta))
	colSums(-w * cbind(gr * x, theta * (digamma(y + theta) - digamma(theta) +
	  log(theta) + 1 - log(mu + theta) - (y + theta)/(mu + theta))))
      }

      ## starting values from previous iteration
      start <- c(component$coef, component$theta)
      ## if not available: use geometric GLM
      if(length(start) < df) start <- c(
        glm.fit(x, y, weights = w, family = MASS::negative.binomial(1), offset = offset)$coefficients,
	0
      )

      ## BFGS optimization
      opt <- optim(par = start, fn = nll, gr = gr, method = "BFGS", control = control)
      
      ## estimated parameters
      cf <- opt$par[-df]
      theta <- exp(opt$par[df])
    }
    
    with(list(coef = cf, theta = theta, df = ncol(x) + is.null(.theta)),
      eval(z@defineComponent))
  }

  z
}
