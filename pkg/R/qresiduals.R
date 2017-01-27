qresiduals <- function(object, ...) {
  UseMethod("qresiduals")
}

qresiduals.default <- function(object, type = c("random", "quantile"), nsim = 1L, prob = 0.5, ...)
{
  ## type of residual for discrete distribution (if any)
  type <- match.arg(type)

  ## if 'object' is not a vector/matrix, apply pit() method
  if(is.object(object) | !is.numeric(object)) {
    object <- try(pit(object))
    if(inherits(object, "try-error")) stop("could not obtain probability integral transform from 'object'")
  }

  ## preprocess supplied probabilities
  nc <- NCOL(object)
  nr <- NROW(object)
  if(nc > 2L) stop("quantiles must either be 1- or 2-dimensional")
  if(nc == 2L) {
    if(type == "random") {
      object <- matrix(
        runif(nr * nsim, min = rep(object[, 1L], nsim), max = rep(object[, 2L], nsim)),
        nrow = nr, ncol = nsim, dimnames = list(rownames(object), paste("r", 1L:nsim, sep = "_"))
      )
    } else {
      nam <- rownames(object)
      object <- object[, 1L]  %*% t(1 - prob) + object[, 2L] %*% t(prob)
      dimnames(object) <- list(nam, paste("q", prob, sep = "_"))
    }
    nc <- NCOL(object)
  }
  if(!is.null(dim(object)) & nc == 1L) object <- drop(object)

  ## compute quantile residuals  
  qnorm(object)
}

pit <- function(object, ...) {
  UseMethod("pit")
}

pit.lm <- function(object, ...)
{
  ## obtain preprocessed response and fitted means
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  mu <- fitted(object)
  
  ## compute probabilities from response distribution
  pr <- pnorm(y, mean = mu, sd = summary(object)$sigma)
  return(pr)
}

pit.glm <- function(object, ...)
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
  return(pr)
}

pit.negbin <- function(object, ...)
{
  ## response and fitted means
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  mu <- fitted(object)
  
  ## probabilities from response distribution
  pr <- cbind(pnbinom(y - 1L, mu = mu, size = object$theta), pnbinom(y, mu = mu, size = object$theta))  
  return(pr)
}

pit.zeroinfl <- pit.hurdle <- function(object, ...)
{
  ## response
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  n <- length(y)
  
  ## corresponding probabilities
  pr <- predict(object, type = "prob", max = max(y))
  pr <- cbind(0, t(apply(pr, 1L, cumsum)))
  pr <- cbind(pr[cbind(1L:n, y + 1L)], pr[cbind(1L:n, y + 2L)])
  return(pr)
}

pit.zerotrunc <- function(object, ...)
{
  ## response
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  n <- length(y)
  
  ## corresponding probabilities
  pr <- predict(object, type = "prob", max = max(y))
  pr <- cbind(0, t(apply(pr, 1L, cumsum)))
  pr <- cbind(pr[cbind(1L:n, y)], pr[cbind(1L:n, y + 1L)])
  return(pr)
}
