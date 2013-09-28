FLXMRnegbin <- function(formula = . ~ ., theta = NULL, offset = NULL)
{
  stopifnot(require("flexmix"))
  .theta <- theta
  
  nbrefit <- function(x, y, w) {
    fit <- c(glm.fit(x, y, weights = w, offset = offset, family = negative.binomial(theta)),
	     list(call = sys.call(), offset = offset,
		  control = eval(formals(glm.fit)$control),
		  method = "weighted.glm.fit"))
    fit$df.null <- sum(w) + fit$df.null - fit$df.residual - fit$rank - is.null(.theta)
    fit$df.residual <- sum(w) - fit$rank - is.null(.theta)
    fit$x <- x
    fit
  }
	      
  z <- new("FLXMRglm", weighted = TRUE, formula = formula,
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
        dnbinom(y, mu = predict(x, ...), size = theta, log = TRUE)
    
      new("FLXcomponent",
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
        dnbinom(y, mu = predict(x, ...), size = theta, log = TRUE)
    
      new("FLXcomponent",
        parameters = list(coef = coef),
        logLik = logLik,
	predict = predict,
        df = df)
    }, as.environment(list(theta = .theta))))
    
  }
      
  z@fit <- function(x, y, w, component){
    theta <- if(is.null(component@parameters$theta)) {
      if(is.null(.theta)) 1 else .theta
    } else {
      nll <- function(th) with(
        glm.fit(x, y, weights = w, family = negative.binomial(th), offset = offset, start = component@parameters$coef),
	-sum(w * dnbinom(y, mu = fitted.values, size = th, log = TRUE)))
      theta <- optimize(nll, c(1/500, 500))$minimum ## FIXME: improve speed
    }
    fit <- glm.fit(x, y, weights = w, family = negative.binomial(theta), offset = offset, start = component@parameters$coef)
    with(list(coef = coef(fit), theta = theta, df = ncol(x) + is.null(.theta)),
	 eval(z@defineComponent))
  }

  z
}


if(FALSE) {
set.seed(1)
d <- data.frame(x = runif(500, -1, 1))
d$cluster <- rep(1:2, each = 250)
d$y <- rnbinom(500, mu = exp(c(1, -1)[d$cluster] + c(0, 3)[d$cluster] * d$x), size = 1)

library("flexmix")
fm1 <- flexmix(y ~ x, data = d, k = 2, model = FLXMRnegbin(theta = 1))
fm0 <- flexmix(y ~ x, data = d, k = 2, model = FLXMRnegbin())

parameters(fm1)
parameters(fm0)

plot(fm1)
plot(fm0)

table(clusters(fm0), clusters(fm1))

plot(posterior(fm0)[,1] ~ factor(cluster), data = d)
plot(posterior(fm1)[,1] ~ factor(cluster), data = d)
plot(posterior(fm0)[,1] ~ posterior(fm1)[,1])

rf1 <- lapply(1:2, function(i) glm(y ~ x, data = d, family = negative.binomial(1), weights = posterior(fm1)[,i]))
rootogram(rf1[[1]])
rootogram(rf1[[2]])

rf0 <- lapply(1:2, function(i) glm.nb(y ~ x, data = d, weights = posterior(fm0)[,i]))
rootogram(rf0[[1]])
rootogram(rf0[[2]])

"+.rootogram" <- function(e1, e2) {
  xlab <- paste(unique(c(attr(e1, "xlab"), attr(e2, "xlab"))), collapse = " / ")
  main <- paste(unique(c(attr(e1, "main"), attr(e2, "main"))), collapse = " / ")
  e1 <- as.data.frame(e1)
  e2 <- as.data.frame(e2)
  e <- e1[e1$x %in% e2$x, ] + e2[e2$x %in% e1$x, ]
  rootogram.default(structure(e$observed, .Names = e$x/2), e$expected,
    main = main, xlab = xlab)
}

rootogram(rf1[[1]]) + rootogram(rf1[[2]])
rootogram(rf0[[1]]) + rootogram(rf0[[2]])
}
