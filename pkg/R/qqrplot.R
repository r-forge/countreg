qqrplot <- function(object, aggregate = c("median", "mean", "none"), nsim = 1L,
  range = FALSE, diag = TRUE, xlim = NULL, ylim = NULL,
  main = "Q-Q residuals plot", xlab = "Theoretical quantiles", ylab = "Quantile residuals",
  ...)
{
  ## compute quantile residuals
  aggregate <- match.arg(aggregate)
  qres <- qresiduals(object, aggregate = aggregate, nsim = nsim)
  if(is.null(dim(qres))) qres <- matrix(qres, ncol = 1L)

  ## corresponding normal quantiles
  q2q <- function(y) qnorm(ppoints(length(y)))[order(order(y))]
  qnor <- apply(qres, 2L, q2q)
    
  ## default plotting ranges
  if(is.null(xlim)) xlim <- range(qnor)
  if(is.null(ylim)) ylim <- range(qres)

  ## set up coordinates
  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)

  ## polygon for range
  if(!identical(range, FALSE)) {
    if(isTRUE(range)) range <- "lightgray"
    rg <- qresiduals(object, aggregate = "range")
    y1 <- sort(rg[,1])
    y2 <- sort(rg[,2])
    x <- c(q2q(y1), rev(q2q(y2)))
    y <- c(y1, rev(y2))
    y[!is.finite(y)] <- 100 * sign(y[!is.finite(y)])
    x[!is.finite(x)] <- 100 * sign(x[!is.finite(x)])
    polygon(x, y, col = range, border = range)
    box()
  }

  ## add Q-Q plot(s)
  for(i in 1L:ncol(qres)) points(qnor[, i], qres[, i], ...)
  
  ## reference diagol
  if(!identical(diag, FALSE)) {
    if(isTRUE(diag)) diag <- "black"
    abline(0, 1, col = diag, lty = 2)
  }
  
  ## return coordinates invisibly
  invisible(list(normal = qnor, residuals = qres))
}
