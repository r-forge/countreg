dhpois <- function(x, lambda, hurdle, log = FALSE) {
  rval <- dpois(x, lambda = lambda, log = TRUE) -
    ppois(0, lambda = lambda, lower.tail = FALSE, log.p = TRUE) + 
    log(hurdle)
  x <- rep(x, length.out = length(rval))
  hurdle <- rep(hurdle, length.out = length(rval))
  rval <- ifelse(x == 0L, log(1 - hurdle), rval)
  if(log) rval else exp(rval)
}

phpois <- function(x, lambda, hurdle, log.p = FALSE) {
  rval <- log(ppois(x, lambda = lambda) - dpois(0L, lambda = lambda)) -
    ppois(0, lambda = lambda, lower.tail = FALSE, log.p = TRUE) + 
    log(hurdle)
  x <- rep(x, length.out = length(rval))
  hurdle <- rep(hurdle, length.out = length(rval))
  rval <- ifelse(x == 0L, 0, exp(rval)) + (1 - hurdle)
  if(log.p) log(rval) else rval
}

shpois <- function(x, lambda, hurdle, parameter = c("lambda", "hurdle"), drop = TRUE) {
  parameter <- sapply(parameter, function(x) match.arg(x, c("lambda", "hurdle")))
  s <- cbind(
    if("lambda" %in% parameter) x/lambda - 1 - exp(-lambda)/(1 - exp(-lambda)) else NULL,
    if("hurdle" %in% parameter) ifelse(x == 0L, -1, 0) else NULL
  )
  colnames(s) <- c("lambda", "hurdle")[c("lambda", "hurdle") %in% parameter]
  if(drop) drop(s) else s
}
