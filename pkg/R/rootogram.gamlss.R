#--------------------------------------------------------------------
# Reviced:  Mikis Stasinopoulos 7-Sep-2017
# It uses the fist gamlss.family name which is unique
# there were mistakes in the second family names 
# (mistakes are are corrected in the latest gamlss.dist)
#--------------------------------------------------------------------
# NOTES 
# Only count distributions are working here 
# Binomial type should be added
# No attemp for continuous is done here 
# QUESTIONS
# Can we have rootogram agaist explanatory variables? 
#-------------------------------------------------------------------
rootogram.gamlss <- function (object, newdata = NULL, breaks = NULL, max = NULL, 
          xlab = NULL, main = NULL, width = NULL, ...) 
{
  family <- substr(family(object)[1], 1L, 30L) # fist family
  if (!(family %in% c( "GEOM", "GEOMo","LG","PO", "YULE", "ZIPF",
                       "NBI", "NBII", "PIG",  "WARING",  "ZALG", 
                       "ZAP", "ZAZIPF", "ZIP", "ZIP2", "GPO", 
                       "DPO", "BNB", "NBF", "DEL", "SICHEL", "SI",
                       "ZANBI", "ZAPIG", "ZINBI", "ZIPIG", "ZANBF",
                       "ZABNB", "ZASICHEL", "ZINBF", "ZIBNB", 
                       "ZISICHEL"))) { # , "PSGIG"
    stop("family currently not supported")
  }
  mt <- terms(object)
  mf <- if (is.null(newdata)) {
    model.frame(object)
  }
  else {
    model.frame(mt, newdata, na.action = na.omit)
  }
   y <- model.response(mf)
   w <- object$weights
  if (is.null(w)) 
      w <- rep(1, NROW(y))
     mu <- predict(object, what = "mu", newdata = newdata, type = "response", 
                na.action = na.omit)
  sigma <- predict(object, what = "sigma", newdata = newdata, 
                   type = "response", na.action = na.omit)
     nu <- predict(object, what = "nu", newdata = newdata, type = "response", 
                na.action = na.omit)
    tau <- predict(object, what = "tau", newdata = newdata, type = "response", 
                   na.action = na.omit)
  if (family == "gaussian") {
         s <- sqrt(weighted.mean(residuals(object)^2, w))
    if (is.null(breaks)) 
    breaks <- "Sturges"
    breaks <- hist(y[w > 0], plot = FALSE, breaks = breaks)$breaks
    obsrvd <- as.vector(xtabs(w ~ cut(y, breaks, include.lowest = TRUE)))
    p <- matrix(NA, nrow = length(y), ncol = length(breaks) - 1L)
    for (i in 1L:ncol(p)) p[, i] <- pnorm(breaks[i + 1L], 
                        mean = mu, sd = s) - pnorm(breaks[i], mean = mu, 
                                                                sd = s)
    expctd <- colSums(p * w)
  }
  else if (family == "binomial") {
    if (NCOL(y) < 2L) 
      y <- cbind(y, 1L - y)
    size <- unique(rowSums(y))
    if (length(size) > 1L) 
      stop("rootogram only applicable to binomial distributions with same size")
    at <- 0L:size
    breaks <- -1L:size + 0.5
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = at)))
    p <- matrix(NA, length(mu), length(at))
    for (i in at) p[, i + 1L] <- dbinom(i, prob = mu, size = size)
    expctd <- colSums(p * w)
  }
  else {
    max0 <- if (is.null(max)) 
      max(1.5 * max(y[w > 0]), 20L)
    else max
    obsrvd <- as.vector(xtabs(w ~ factor(y, levels = 0L:max0)))
        at <- 0L:max0
         p <- matrix(NA, length(mu), length(at))
         
    if (family == "GEOM") {
           for (i in at) p[, i + 1L] <- dGEOM(i, mu = mu)
          }
    if (family == "GEOMo") {
           for (i in at) p[, i + 1L] <- dGEOMo(i, mu = mu)
         } 
    if (family == "LG") {
           for (i in at) p[, i + 1L] <- dLG(i, mu = mu)
         } 
    if (family == "PO") {
      for (i in at) p[, i + 1L] <- dpois(i, lambda = mu)
    }
    if (family == "YULE") {
           for (i in at) p[, i + 1L] <- dYULE(i, mu = mu)
         }     
    if (family == "ZIPF") {
           for (i in at) p[, i + 1L] <- dZIPF(i, mu = mu)
         }       
    if (family == "NBI") {
      for (i in at) p[, i + 1L] <- gamlss.dist::dNBI(i, mu = mu, sigma = sigma)
    }
    if (family == "NBII") {
      for (i in at) p[, i + 1L] <- gamlss.dist::dNBII(i, mu = mu, sigma = sigma)
      }
    if (family == "PIG") {
      for (i in at) p[, i + 1L] <- gamlss.dist::dPIG(i, mu = mu, sigma = sigma)
      }
    if (family == "WARING") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dWARING(i, mu = mu, sigma = sigma)
         }
    if (family == "ZALG") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZALG(i, mu = mu, sigma = sigma)
         }  
    if (family == "ZAP") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZAP(i, mu = mu, sigma = sigma)
         }   
    if (family == "ZAZIPF") {
           for (i in at) p[, i + 1L] <- gamlss.dist::ZAZIPF(i, mu = mu, sigma = sigma)
         }      
    if (family == "ZIP") {
      for (i in at) p[, i + 1L] <- gamlss.dist::dZIP(i, mu = mu, sigma = sigma)
         }
    if (family == "ZIP2") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZIP2(i, mu = mu, sigma = sigma)
         }
      if (family == "GPO") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dGPO(i, mu = mu, sigma = sigma)
         }
     if (family == "DPO") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dDPO(i, mu = mu, sigma = sigma)
         }        
      if (family == "BNB") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dBNB(i, mu = mu, sigma = sigma, nu = nu)
         }
      if (family == "DEL") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dDEL(i, mu = mu, sigma = sigma, nu = nu)
         }   
      if (family == "SICHEL") {
      for (i in at) p[, i + 1L] <- gamlss.dist::dSICHEL(i, mu = mu, sigma = sigma, nu = nu)
         }
      if (family == "SI") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dSICHEL(i, mu = mu, sigma = sigma, nu = nu)
         }
      if (family == "ZANBI") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZANBI(i, mu = mu, sigma = sigma, nu = nu)
         }   
      if (family == "ZAPIG") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZAPIG(i, mu = mu, sigma = sigma, nu = nu)
         }  
      if (family == "ZINBI") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZINBI(i, mu = mu, sigma = sigma, nu = nu)
         } 
      if (family == "ZIPIG") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZIPIG(i, mu = mu, sigma = sigma, nu = nu)
         }      
      if (family == "ZANBF") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZANBF(i, mu = mu, sigma = sigma, nu = nu, tau=tau)
         }    
      if (family == "ZABNB") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZABNB(i, mu = mu, sigma = sigma, nu = nu, tau=tau)
         }    
      if (family == "ZASICHEL") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZASICHEL(i, mu = mu, sigma = sigma, nu = nu, tau=tau)
         }    
      if (family == "ZINBF") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZINBFL(i, mu = mu, sigma = sigma, nu = nu, tau=tau)
         } 
      if (family == "ZIBNB") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZIBNB(i, mu = mu, sigma = sigma, nu = nu, tau=tau)
         } 
      if (family == "ZISICHEL") {
           for (i in at) p[, i + 1L] <- gamlss.dist::dZISICHEL(i, mu = mu, sigma = sigma, nu = nu, tau=tau)
         }   
    expctd <- colSums(p * w)
    if (is.null(max)) {
      max <- if (all(expctd >= 1L)) 
        max0
      else max(ceiling(mean(y)), min(which(expctd < 1L)) - 
                 1L)
      max <- min(max, length(expctd) - 1L)
    }
    breaks <- -1L:max + 0.5
    obsrvd <- obsrvd[1L:(max + 1L)]
    expctd <- expctd[1L:(max + 1L)]
  }
  if (is.null(xlab)) 
    xlab <- as.character(attr(mt, "variables"))[2L]
  if (is.null(main)) 
    main <- deparse(substitute(object))
  rootogram.default(obsrvd, expctd, breaks = breaks, xlab = xlab, 
                    main = main, width = if (family == "gaussian") 
                      1
                    else 0.9, ...)
}