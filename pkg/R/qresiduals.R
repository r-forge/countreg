qresiduals <- function(object, ...) {
  UseMethod("qresiduals")
}

qresiduals.lm <- function(object, ...)
{
  ## obtain preprocessed response and fitted means
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  mu <- fitted(object)
  
  ## compute probabilities from response distribution
  pr <- pnorm(y, mean = mu, sd = summary(object)$sigma)

  ## call default method
  qresiduals(pr, ...)
}

qresiduals.glm <- function(object, ...)
{
  ## obtain preprocessed response and fitted means
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  weights <- weights(object)
  nobs <- nobs(object)
  etastart <- NULL
  mustart <- NULL
  n <- NULL
  eval(object$family$initialize)
  mu <- fitted(object)
  
  ## compute probabilities from response distribution
  pr <- switch(object$family$family,
    "gaussian" = {
      pnorm(y, mean = mu, sd = sqrt(sum((y - mu)^2) / df.residual(object)))
    },
    "poisson" = {
      cbind(ppois(y - 1L, lambda = mu), ppois(y, lambda = mu))
    },
    "binomial" = {
      y <- y * weights / n
      if(!isTRUE(all.equal(as.numeric(y), as.numeric(round(y))))) {
        stop("binomial quantile residuals require integer response")
      }
      y <- round(y)
      cbind(pbinom(y - 1L, size = n, prob = mu), pbinom(y, size = n, prob = mu))
    },
    stop("not implemented yet")
  )
  
  ## call default method
  qresiduals(pr, ...)
}

qresiduals.negbin <- function(object, ...)
{
  ## response and fitted means
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  mu <- fitted(object)
  
  ## probabilities from response distribution
  pr <- cbind(pnbinom(y - 1L, mu = mu, size = object$theta), pnbinom(y, mu = mu, size = object$theta))  

  ## call default method
  qresiduals(pr, ...)
}

qresiduals.zeroinfl <- qresiduals.hurdle <- function(object, ...)
{
  ## response
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  n <- length(y)
  
  ## corresponding probabilities
  pr <- predict(object, type = "prob", max = max(y))
  pr <- cbind(0, t(apply(pr, 1L, cumsum)))
  pr <- cbind(pr[cbind(1L:n, y + 1L)], pr[cbind(1L:n, y + 2L)])

  ## call default method
  qresiduals(pr, ...)
}

qresiduals.zerotrunc <- function(object, ...)
{
  ## response
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  n <- length(y)
  
  ## corresponding probabilities
  pr <- predict(object, type = "prob", max = max(y))
  pr <- cbind(0, t(apply(pr, 1L, cumsum)))
  pr <- cbind(pr[cbind(1L:n, y)], pr[cbind(1L:n, y + 1L)])

  ## call default method
  qresiduals(pr, ...)
}

qresiduals.default <- function(object, aggregate = c("median", "mean", "range", "none"),
  nsim = 10L, ...)
{
  ## type of aggregation (if any)
  aggregate <- match.arg(aggregate)

  ## preprocess supplied probabilities
  nc <- NCOL(object)
  nr <- NROW(object)
  if(nc > 2L) stop("quantiles must either be 1- or 2-dimensional")
  if(nc == 2L & aggregate == "median") {
    object <- rowMeans(object)
    nc <- 1L
  }
  if(aggregate == "range") {
    if(nc == 1L) object <- cbind(object, object)
    colnames(object) <- c("lower", "upper")
  }
  if(!is.null(dim(object)) & nc == 1L) object <- drop(object)

  ## draw random probabilities (if necessary)
  if(nc == 2L & aggregate %in% c("mean", "none")) {
    object <- matrix(
      runif(nr * nsim, min = rep(object[, 1L], nsim), max = rep(object[, 2L], nsim)),
      nrow = nr, ncol = nsim, dimnames = list(rownames(object), paste0("s", 1L:nsim))
    )
  }

  ## compute quantile residuals  
  qres <- qnorm(object)
  if(aggregate == "mean") qres <- rowMeans(qres)  
  return(qres)
}
