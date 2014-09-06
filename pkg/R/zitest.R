zitest <- function(object, type = c("scoreZIP")){
 
 ## sanity checks
 if(!inherits(object, "glm") || family(object)$family != "poisson")
    stop("only Poisson GLMs can be tested")
 type <- match.arg(type)

 n <- nobs(object)
 y <- if (is.null(object$y)) 
       model.response(model.frame(object))
      else object$y

 pzero <- exp(-fitted(object))
 s1    <- ( (y == 0) - pzero )/pzero
 s2    <- ( 1 - pzero )/pzero
 stat  <- sum(s1)^2 / ( sum(s2) - n * mean(y))

 rval <- list(statistic = c(S = stat), 
              p.value = pchisq(stat, df = 1, lower.tail = FALSE)/2, 
              alternative = NULL, 
              method = "Zero inflation test", 
              data.name = deparse(substitute(object))
              )
 
 class(rval) <- "htest"
 return(rval)
}
