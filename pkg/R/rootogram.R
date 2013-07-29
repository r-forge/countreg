rootogram <- function(object, ...) {
  UseMethod("rootogram")
}

rootogram.default <- function(object, fitted, breaks = NULL,
  style = c("hanging", "standing", "deviation"),
  scale = c("sqrt", "raw"), width = NULL,
  xlab = NULL, ylab = NULL, main = NULL)
{
  ## rectangle style
  scale <- match.arg(scale)
  style <- match.arg(style)

  ## default annotation
  if(is.null(xlab)) {
    xlab <- if(is.null(names(dimnames(object)))) {
      deparse(substitute(object))
    } else {
      names(dimnames(object))[1L]
    }
  }
  if(is.null(ylab)) {
    ylab <- if(scale == "raw") "Frequency" else "sqrt(Frequency)" 
  }
  if(is.null(main)) main <- deparse(substitute(fitted))
  
  ## breaks, midpoints, widths
  if(is.null(breaks)) {
    x <- as.numeric(names(object))
    if(length(x) < 1L) x <- 0L:(length(object) - 1L)
    breaks <- (head(x, -1L) + tail(x, -1L))/2
    breaks <- c(2 * head(x, 1L) - head(breaks, 1L), breaks,
      2 * tail(x, 1L) - tail(breaks, 1L))
    if(is.null(width)) width <- 0.9
  } else {
    x <- (head(breaks, -1L) + tail(breaks, -1L))/2
    if(is.null(width)) width <- 1
  }

  ## raw vs. sqrt scale
  if(scale == "sqrt") {
    obsrvd <- sqrt(as.vector(object))
    expctd <- sqrt(as.vector(fitted))
  } else {
    obsrvd <- as.vector(object)
    expctd <- as.vector(fitted)
  }

  ## height/position of rectangles
  y <- if(style == "hanging") expctd - obsrvd else 0
  height <- if(style == "deviation") expctd - obsrvd else obsrvd

  ## return everything as data.frame
  rval <- data.frame(observed = as.vector(object), expected = as.vector(fitted),
    x = x, y = y, width = diff(breaks) * width, height = height,
    line = expctd)
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  class(rval) <- c("rootogram", "data.frame")
  return(rval)
}

c.rootogram <- rbind.rootogram <- function(...)
{
  ## list of rootograms
  rval <- list(...)
  
  ## group sizes
  for(i in seq_along(rval)) {
    if(is.null(rval[[i]]$group)) rval[[i]]$group <- 1L
  }
  n <- lapply(rval, function(r) table(r$group))

  ## labels
  xlab <- unlist(lapply(rval, function(r) attr(r, "xlab")))
  ylab <- unlist(lapply(rval, function(r) attr(r, "ylab")))
  nam <- names(rval)
  main <- if(is.null(nam)) {
    as.vector(sapply(rval, function(r) attr(r, "main")))
  } else {
    make.unique(rep.int(nam, sapply(n, length)))
  }
  n <- unlist(n)

  ## combine and return
  rval <- do.call("rbind.data.frame", rval)
  rval$group <- if(length(n) < 2L) NULL else rep.int(seq_along(n), n)
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  class(rval) <- c("rootogram", "data.frame")
  return(rval)
}

print.rootogram <- function(x, plot = TRUE, ...) {
  if(plot) {
    plot(x, ...)
  } else {
    print.data.frame(x, ...)
  }
  invisible(x)
}

plot.rootogram <- function(x,
  xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, main = NULL,
  border = "black", fill = "lightgray", col = "#B61A51",
  lwd = 2, pch = 19, lty = 1, type = NULL, ...)
{
  ## handling of groups
  if(is.null(x$group)) x$group <- 1L
  n <- max(x$group)
  if(is.null(type)) type <- ifelse(any(table(x$group) > 15L), "l", "b")

  ## annotation
  if(is.null(xlab)) xlab <- TRUE
  if(is.null(ylab)) ylab <- TRUE
  if(is.null(main)) main <- TRUE
  xlab <- rep(xlab, length.out = n)
  ylab <- rep(ylab, length.out = n)
  main <- rep(main, length.out = n)
  if(is.logical(xlab)) xlab <- ifelse(xlab, attr(x, "xlab"), "")
  if(is.logical(ylab)) ylab <- ifelse(ylab, attr(x, "ylab"), "")
  if(is.logical(main)) main <- ifelse(main, attr(x, "main"), "")

  ## plotting function
  rootogram1 <- function(d, ...) {
    ## rect elements
    xleft <- d$x - d$width/2
    xright <- d$x + d$width/2
    ybottom <- d$y
    ytop <- d$y + d$height
    j <- unique(d$group)
    
    ## defaults
    if(is.null(xlim)) xlim <- range(c(xleft, xright))
    if(is.null(ylim)) ylim <- range(c(ybottom, ytop, d$line))

    ## draw rootogram
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
      xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...)
    axis(1)
    axis(2)
    rect(xleft, ybottom, xright, ytop, border = border, col = fill)
    abline(h = 0, col = border)
    lines(d$x, d$line,
      col = col, pch = pch, type = type, lty = lty, lwd = lwd)
   }
   
   ## draw plots
   if(n > 1L) par(mfrow = n2mfrow(n))
   for(i in 1L:n) rootogram1(x[x$group == i, ], ...)
}

rootogram.numeric <- function(object, fitted, breaks = NULL,
  start = NULL, width = NULL, xlab = NULL, ylab = NULL, main = NULL, ...)
{
  ## distribution to be fitted
  dist <- fitted
  if(!is.character(fitted)) fitted <- deparse(substitute(fitted))
  if(substr(fitted, 1L, 1L) != "d") {
    fitted <- match.arg(tolower(fitted),
      c("beta", "cauchy", "chi-squared", "chisquared", "exponential", "f",
        "gamma", "geometric", "log-normal", "lognormal", "logistic", "negative binomial",
	"negbin", "normal", "gaussian", "poisson", "t", "weibull"))
    fitted <- switch(fitted,
      "chisquared" = "chi-squared",
      "lognormal" = "log-normal",
      "negbin" = "negative binomial",
      "gaussian" = "normal",
      fitted)
    if(is.character(dist)) dist <- fitted      
  }

  ## labels
  if(is.null(xlab)) xlab <- deparse(substitute(object))
  if(is.null(main)) main <- sprintf('fitdistr(%s, "%s")', deparse(substitute(object)), fitted)

  ## call MASS::fitdistr
  xfit <- suppressWarnings(try(fitdistr(object, dist, start = start), silent = TRUE))
  if(!inherits(xfit, "fitdistr")) stop("could not obtain fitted distribution")

  ## fitted probability distribution function
  pdist <- switch(fitted,
    "beta" = pbeta,
    "cauchy" = pcauchy, 
    "chi-squared" = pchisq,
    "exponential" = pexp,
    "f" = pf, 
    "gamma" = pgamma,
    "geometric" = pgeom,
    "log-normal" = plnorm, 
    "logistic" = plogis,
    "negative binomial" = pnbinom, 
    "normal" = pnorm,
    "poisson" = ppois,
    "t" = function(x, m, s, df) pt((x - m)/s, df),
    "weibull" = pweibull, 
    paste("p", substr(fitted, 2L, nchar(fitted)), sep = ""))
  if(is.character(pdist)) pdist <- try(get(pdist), silent = TRUE)
  if(!is.function(pdist)) stop("invalid specification of fitted distribution")

  ## second argument should be the full parameter vector
  f <- formals(pdist)
  args <- names(f)
  m <- match(names(xfit$estimate), args)
  formals(pdist) <- c(f[c(1, m)], f[-c(1, m)])
  pfun <- function(x, parm, ...) pdist(x, parm, ...)
  l <- length(xfit$estimate)
  if(l > 1L) {
    body(pfun) <- parse(text = paste("pdist(x, ", paste("parm[", 1L:l, "]", collapse = ", "), ", ...)"))
  }

  ## different default breaks for discrete distributions
  if(is.null(breaks)) {
    if(tolower(fitted) %in% c("geometric", "negative binomial", "poisson", "binomial")) {
      breaks <- -1L:max(object) + 0.5
      if(is.null(width)) width <- 0.9
    } else {
      breaks <- "Sturges"
    }
  }

  ## observed and expected frequencies
  xhist <- hist(object, plot = FALSE, breaks = breaks)
  expctd <- xfit$n * (pfun(tail(xhist$breaks, -1L), xfit$estimate) -
    pfun(head(xhist$breaks, -1L), xfit$estimate))

  ## call base rootogram function
  rootogram.default(xhist$counts, expctd, breaks = xhist$breaks,
    xlab = xlab, main = main, width = width, ...)
}

rootogram.zeroinfl <- rootogram.hurdle <- function(object, newdata = NULL,
  max = NULL, xlab = NULL, main = NULL, width = 0.9, ...)
{
  ## observed response
  mt <- terms(object)
  mf <- if(is.null(newdata)) {
    model.frame(object)
  } else {
    model.frame(mt, newdata, na.action = na.omit)
  }
  y <- model.response(mf)
  
  ## observed and expected frequencies
  max0 <- if(is.null(max)) max(1.5 * max(y), 20L) else max  
  obsrvd <- table(factor(y, levels = 0L:max0))
  expctd <- if(is.null(newdata)) {
    colSums(predict(object, type = "prob", at = 0L:max0))  
  } else {
    colSums(predict(object, newdata = newdata, type = "prob", at = 0L:max0, na.action = na.omit))
  }

  ## try to guess a good maximum
  if(is.null(max)) {
    max <- if(all(expctd >= 1L)) max0 else max(ceiling(mean(y)), min(which(expctd < 1L)) - 1L)
    max <- min(max, length(expctd) - 1L)
  }

  ## observed and expected frequencies
  obsrvd <- obsrvd[1L:(max + 1L)]
  expctd <- expctd[1L:(max + 1L)]
  
  if(is.null(xlab)) xlab <- as.character(attr(mt, "variables"))[2L]
  if(is.null(main)) main <- deparse(substitute(object))
  rootogram.default(obsrvd, expctd, breaks = -1L:max + 0.5,
    xlab = xlab, main = main, width = width, ...)  
}

rootogram.glm <- function(object, newdata = NULL, breaks = NULL,
  max = NULL, xlab = NULL, main = NULL, width = NULL, ...) 
{
  family <- substr(family(object)$family, 1L, 17L)
  if(!(family %in% c("negbin", "Negative Binomial", "poisson", "binomial", "gaussian"))) {
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
  mu <- predict(object, newdata = newdata, type = "response", na.action = na.omit)

  if(family == "gaussian") {
    ## estimated standard deviation (ML)
    s <- sqrt(mean(residuals(object)^2))

    ## breaks
    if(is.null(breaks)) breaks <- "Sturges"
    yhist <- hist(y, plot = FALSE, breaks = breaks)
    breaks <- yhist$breaks
    obsrvd <- yhist$count

    ## expected frequencies
    p <- matrix(NA, nrow = length(y), ncol = length(breaks) - 1L)
    for(i in 1L:ncol(p)) p[, i] <- pnorm(yhist$breaks[i + 1L], mean = mu, sd = s) -
      pnorm(yhist$breaks[i], mean = mu, sd = s)
    expctd <- colSums(p)
  } else if(family == "binomial") {
    ## successes and failures
    if(NCOL(y) < 2L) y <- cbind(y, 1L - y)

    ## number of attempts
    size <- unique(rowSums(y))
    if(length(size) > 1L) stop("rootogram only applicable to binomial distributions with same size")
    at <- 0L:size
    breaks <- -1L:size + 0.5
    
    ## observed and expected
    obsrvd <- table(factor(y[, 1L], levels = at))
    p <- matrix(NA, length(mu), length(at))
    for(i in at) p[, i + 1L] <- dbinom(i, prob = mu, size = size)
    expctd <- colSums(p)
  } else {
    ## observed frequencies
    max0 <- if(is.null(max)) max(1.5 * max(y), 20L) else max  
    obsrvd <- table(factor(y, levels = 0L:max0))

    ## expected frequencies
    at <- 0L:max0
    p <- matrix(NA, length(mu), length(at))
    if(family == "poisson") {
      for(i in at) p[, i + 1L] <- dpois(i, lambda = mu)
    } else {
      theta <- object$theta
      if(is.null(theta)) theta <- get(".Theta", environment(family(object)$variance))
      for(i in at) p[, i + 1L] <- dnbinom(i, mu = mu, size = theta)
    }
    expctd <- colSums(p)

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

autoplot.rootogram <- function(object,
  colour = c("black", "#B61A51"), fill = "darkgray", size = c(1.2, 4), ...)
{
  ## determine grouping
  class(object) <- "data.frame"
  if(is.null(object$group)) object$group <- 1L
  n <- max(object$group)
  object$group <- factor(object$group, levels = 1L:n, labels = attr(object, "main"))

  ## rectangles and fitted lines
  rval <- ggplot(object, aes(xmin = x - width/2, xmax = x + width/2, ymin = y, ymax = y + height, x = x, y = line)) +
    geom_rect(colour = colour[1L], fill = fill) + geom_line(colour = colour[2L], size = size[1L]) +
    geom_hline(yintercept = 0)
  if(all(table(object$group) <= 15L)) rval <- rval + geom_point(colour = colour[2L], size = size[2L])

  ## grouping (if any)
  if(n > 1L) rval <- rval + facet_grid(group ~ .)
  
  ## annotation
  rval <- rval + xlab(paste(unique(attr(object, "xlab")), collapse = "/")) +
    ylab(paste(unique(attr(object, "ylab")), collapse = "/"))

  ## return with annotation
  rval
}
