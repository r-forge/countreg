qqrplot <- function(object, type = c("random", "quantile"),
  nsim = 1L, prob = 0.5, range = FALSE, diag = TRUE,
  col = "black", fill = "lightgray", xlim = NULL, ylim = NULL,
  main = "Q-Q residuals plot", xlab = "Theoretical quantiles", ylab = "Quantile residuals",
  ...)
{
  ## compute quantile residuals
  qres <- qresiduals(object, type = type, nsim = nsim, prob = prob)
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
    if(isTRUE(range)) range <- c(0.01, 0.99)
    rg <- qresiduals(object, type = "quantile", prob = range)
    y1 <- sort(rg[,1])
    y2 <- sort(rg[,2])
    x <- c(q2q(y1), rev(q2q(y2)))
    y <- c(y1, rev(y2))
    y[!is.finite(y)] <- 100 * sign(y[!is.finite(y)])
    x[!is.finite(x)] <- 100 * sign(x[!is.finite(x)])
    polygon(x, y, col = fill, border = fill)
    box()
  }

  ## add Q-Q plot(s)
  for(i in 1L:ncol(qres)) points(qnor[, i], qres[, i], col = col, ...)
  
  ## reference diagol
  if(!identical(diag, FALSE)) {
    if(isTRUE(diag)) diag <- "black"
    abline(0, 1, col = diag, lty = 2)
  }
  
  ## return coordinates invisibly
  invisible(list(normal = qnor, residuals = qres))
}
