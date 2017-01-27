pithist <- function(object, type = c("random", "proportional"), nsim = 1L,
  breaks = NULL, xlim = c(0, 1), ylim = NULL,
  xlab = "PIT", ylab = "Density", main = NULL,
  border = "black", fill = "lightgray", col = "#B61A51",
  lwd = 2, lty = 1, freq = FALSE, ...)
{
  ## either compute proportion exactly (to do...) or approximate by simulation
  type <- match.arg(type, c("random", "proportional"))
  if(type == "proportional") {
    stop("not yet implemented")
  } else {
    p <- pnorm(qresiduals.default(object, nsim = nsim))
  }

  ## breaks
  if(is.null(breaks)) breaks <- c(4, 10, 20, 25)[cut(NROW(p), c(0, 50, 5000, 1000000, Inf))]
  if(length(breaks) == 1L) breaks <- seq(xlim[1L], xlim[2L], length.out = breaks + 1L)

  ## labels
  if(is.null(main)) main <- deparse(substitute(object))

  ## histogram
  rval <- hist(p, breaks = breaks, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main,
    col = fill, border = border, freq = freq, ...)
  abline(h = 1, col = col, lty = lty, lwd = lwd)
  invisible(rval)
}
