# Note: family(object)[2] gives "Poisson Zero Inflated" for both ZIP and ZAP,
# also "Sichel" for both SI and SICHEL
rootogram.gamlss <- function(object, newdata = NULL, breaks = NULL,
                             max = NULL, xlab = NULL, main = NULL, width = NULL, ...) 
{
  family <- substr(family(object)[2], 1L, 30L)
  if(!(family %in% c("Negative Binomial type I", "Negative Binomial type II", "Poisson", 
                     "Poisson.Inverse.Gaussian", "Sichel", "Yule", "Delaporte", "Logarithmic", 
                     "Zero Adjusted Logarithmic", "Poisson Zero Inflated", "Zero Inflated Poisson 2"))) {
    stop("family currently not supported")
  }
  
  ## observed response
  mt <- terms(object)
  mf <- if(is.null(newdata)) {
    model.frame(object)
  } else {
    model.frame(mt, newdata, na.action = na.omit)
  }
  y <- model.response(mf)
  w <- object$weights
  if(is.null(w)) w <- rep(1, NROW(y))
  mu <- predict(object, what = "mu", newdata = newdata, type = "response", na.action = na.omit)
  sigma <- predict(object, what = "sigma", newdata = newdata, type = "response", na.action = na.omit)
  nu <- predict(object, what = "nu", newdata = newdata, type = "response", na.action = na.omit)
  
  if(family == "gaussian") {
    ## estimated standard deviation (ML)
    s <- sqrt(weighted.mean(residuals(object)^2, w))
    
    ## breaks
    if(is.null(breaks)) breaks <- "Sturges"
    breaks <- hist(y[w > 0], plot = FALSE, breaks = breaks)$breaks
    obsrvd <- as.vector(xtabs(w ~ cut(y, breaks, include.lowest = TRUE)))
    
    ## expected frequencies
    p <- matrix(NA, nrow = length(y), ncol = length(breaks) - 1L)
    for(i in 1L:ncol(p)) p[, i] <- pnorm(breaks[i + 1L], mean = mu, sd = s) -
      pnorm(breaks[i], mean = mu, sd = s)
    expctd <- colSums(p * w)
  } else if(family == "binomial") {
    ## successes and failures
    if(NCOL(y) < 2L) y <- cbind(y, 1L - y)
    
    ## number of attempts
    size <- unique(rowSums(y))
    if(length(size) > 1L) stop("rootogram only applicable to binomial distributions with same size")
    at <- 0L:size
    breaks <- -1L:size + 0.5
    
    ## observed and expected
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = at)))
    p <- matrix(NA, length(mu), length(at))
    for(i in at) p[, i + 1L] <- dbinom(i, prob = mu, size = size)
    expctd <- colSums(p * w)
  } else {
    ## observed frequencies
    max0 <- if(is.null(max)) max(1.5 * max(y[w > 0]), 20L) else max  
    obsrvd <- as.vector(xtabs(w ~ factor(y, levels = 0L:max0)))
    
    ## expected frequencies
    at <- 0L:max0
    p <- matrix(NA, length(mu), length(at))
    if(family == "Poisson") {
      for(i in at) p[, i + 1L] <- dpois(i, lambda = mu)
    } 
    if(family == "Poisson.Inverse.Gaussian") {
      for(i in at) p[, i + 1L] <- gamlss.dist::dPIG(i, mu = mu, sigma = sigma)
    } 
    if(family == "Negative Binomial type II") {
      for(i in at) p[, i + 1L] <- gamlss.dist::dNBII(i, mu = mu, sigma = sigma)
    } 
    if(family == "Negative Binomial type I") {
      for(i in at) p[, i + 1L] <- gamlss.dist::dNBI(i, mu = mu, sigma = sigma)
    } 
    if(family == "YULE") {
      for(i in at) p[, i + 1L] <- gamlss.dist::dYULE(i, mu = mu)
    } 
    if(family == "SICHEL") {
      for(i in at) p[, i + 1L] <- gamlss.dist::dSICHEL(i, mu = mu, sigma = sigma, nu = nu)
    } 
    if(family == "DEL") {
      for(i in at) p[, i + 1L] <- gamlss.dist::dDEL(i, mu = mu, sigma = sigma, nu = nu)
    } 
    expctd <- colSums(p * w)
    
    ## try to guess a good maximum
    if(is.null(max)) {
      max <- if(all(expctd >= 1L)) max0 else max(ceiling(mean(y)), min(which(expctd < 1L)) - 1L)
      max <- min(max, length(expctd) - 1L)
    }
    breaks <- -1L:max + 0.5
    
    ## observed and expected frequencies
    obsrvd <- obsrvd[1L:(max + 1L)]
    expctd <- expctd[1L:(max + 1L)]
  }
  
  if(is.null(xlab)) xlab <- as.character(attr(mt, "variables"))[2L]
  if(is.null(main)) main <- deparse(substitute(object))
  rootogram.default(obsrvd, expctd, breaks = breaks,
                    xlab = xlab, main = main,
                    width = if(family == "gaussian") 1 else 0.9, ...)  
}

