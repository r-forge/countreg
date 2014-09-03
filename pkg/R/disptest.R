disptest <- function(object, 
  type = c("lrtNB2", "scoreNB2", "scoreNB2adj", "scoreNB1", "scoreNB1adj", "scoreKatz"), 
  trafo = NULL, alternative = c("greater", "two.sided", "less"))
{
## sanity checks  
  if(!inherits(object, "glm") || family(object)$family != "poisson")
    stop("only Poisson GLMs can be tested")  ## maybe later GLMs more generally?
  alternative <- match.arg(alternative)
  if(type != "scoreKatz") alternative <- "greater"  ## warning: quick fix
  type <- match.arg(type)

## preprocessing
  y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
  yhat <- fitted(object)
  n <- nobs(object)

## test statistics
  if(type == "lrtNB2") {
  m_nb <- glm.nb(object$formula, data = object$data)
  STAT <- 2 * (logLik(m_nb) - logLik(object))
  pval <- pchisq(STAT, df = 1, lower.tail = FALSE)/2 
  }
   
  if(type == "scoreNB2") {
  STAT <- sum((y - yhat)^2 - y) / sqrt(2 * sum(yhat^2))
  pval <- pnorm(STAT, lower.tail = FALSE)
  }
    
  if(type == "scoreNB2adj") {
  h <- hatvalues(object)
  STAT <- sum((y - yhat)^2 - y  + h * yhat) / sqrt(2 * sum(yhat^2)) 
  pval <- pnorm(STAT, lower.tail = FALSE)  
  }

  if(type == "scoreNB1") {
  STAT <- sum(( (y - yhat)^2 - y) / yhat ) / sqrt(2 * n)
  pval <- pnorm(STAT, lower.tail = FALSE)
  }  

  if(type == "scoreNB1adj") {
  h <- hatvalues(object)
  STAT <- sum(( (y - yhat)^2 - y + h * yhat) / yhat ) / sqrt(2 * n)
  pval <- pnorm(STAT, lower.tail = FALSE)
  }

  if(type == "scoreKatz") { 
  STAT <- sum(( (y - 1) * y  - yhat^2) / yhat ) / sqrt(2 * n)
  pval <- switch(alternative, "greater"   = pnorm(STAT, lower.tail = FALSE),
		                      "two.sided" = pnorm(abs(STAT), lower.tail = FALSE) * 2,
		                      "less"      = pnorm(STAT) )
  }
  
## collect output  
  rval <- list(
           statistic = c(z = STAT),
           p.value = pval,
#	       estimate = EST,
#	       null.value = NVAL,
	       alternative = alternative,
	       method = switch(alternative,
	                  "greater"   = "Overdispersion test",
			          "two.sided" = "Dispersion test",
			          "less"      = "Underdispersion test"),
	       data.name = deparse(substitute(object))
  )
  class(rval) <- "htest"
  return(rval)
}

## TODO:
## 2. also fix DCluster::test.nb.pois() and pscl::odTest()
## proposed interface:
##   poistest(object, object2 = NULL)
## where either a "negbin" and a "glm" object have to be
## supplied or only one of them, then update via either
##   cl <- object$call
##   cl[[1]] <- as.name("glm.nb")
##   cl$link <- object$family$link
##   cl$family <- NULL
## or
##   cl <- object$call
##   cl[[1]] <- as.name("glm")
##   cl$family <- call("poisson")
##   cl$family$link <- object$family$link
##   cl$link <- NULL
##   cl$init.theta <- NULL
## and evaluate the call "cl" appropriately.
