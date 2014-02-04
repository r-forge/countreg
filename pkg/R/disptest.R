disptest <- function(object, type = c("DL", "DLadj"), trafo = NULL, 
  alternative = c("greater", "two.sided", "less"))
{
## sanity checks  
  if(!inherits(object, "glm") || family(object)$family != "poisson")
    stop("only Poisson GLMs can be tested")  ## maybe later GLMs more generally?
  alternative <- match.arg(alternative)
  type <- match.arg(type)

## preprocessing
  y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
  yhat <- fitted(object)

## test statistics
  if(type == "DL") 
  STAT <- sum((y - yhat)^2 - y) / sqrt(2 * sum(yhat^2))
  
  if(type == "DLadj") {
  h <- hatvalues(object)
  STAT <- sum((y - yhat)^2 - y  + h * yhat) / sqrt(2 * sum(yhat^2)) 
  }

## collect output  
  rval <- list(
           statistic = c(z = STAT),
           p.value = switch(alternative,
                      "greater"   = pnorm(STAT, lower.tail = FALSE),
		              "two.sided" = pnorm(abs(STAT), lower.tail = FALSE) * 2,
		              "less"      = pnorm(STAT) ),
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
## 1. compare with DCluster::DeanB() and DCluster::DeanB2()
## and unify
##
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
