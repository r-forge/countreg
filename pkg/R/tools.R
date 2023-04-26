## some helper functions

# d and p function for NB1
dnbinom1 <- function(x, mu, size, ...){
  dnbinom(x, size = mu*size, mu = mu, ...)
}

pnbinom1 <- function(x, mu, size, ...){
  pnbinom(x, size = mu*size, mu = mu, ...)
}



# linear predictor from covariates-matrix, coefficients and offset
linPred <- function(covars, par, offset) {
  as.vector(covars %*% par + offset)
}



# extend link functions with 2nd derivatives
make.link2 <- function(link) {
  linkobj <- make.link(link)
  linkobj$dmu.eta <- switch(link, # 2nd derivative of link function
                            "log" = function(x) exp(x),
                            "logit" = function(x) ( exp(x)*(1-exp(x)) )  /  ( (1+exp(x))^3 ),
                            NULL)
  return(linkobj)
}



## Hessians
# build hessian for distribution with one parameter
hessMat1 <- function(hm, gr, X, eta, dlinkX, ddlinkX, weights = 1) {
  # hm: hessian
  # gr: gradient
  # X: design matrix of first linear link (link for mean parameter)
  # eta: linear link
  # dlinkX: derivative of link function
  # ddlinkX: 2nd derivative of link function
  
  hessMu <- weights*hm*dlinkX(eta)^2 + gr*ddlinkX(eta)
  crossprod(hessMu*X, X)
}

# build hessian for distribution with two parameters
hessMat2 <- function(hm, gr, X, Z,
                    etaX, etaZ, fix = FALSE, dlinkX, ddlinkX, dlinkZ, ddlinkZ, weights = 1) {
  # hm: hessian
  # gr: gradient
  # X: design matrix of first linear link (link for mean parameter)
  # Z: design matrix (link for size parameter, e.g. theta in negative binomial), often Z is column of ones (constant theta)
  # eta: first linear link
  # etaTheta: second linear link
  # fix: is second linear link specified beforehand? (constant --> no derivatives)
  # dlinkX: derivative of first link function
  # ddlinkX: 2nd derivative of first link function
  # dlinkZ: derivative of second link function
  # ddlinkZ: second derivative of second link function
  #
  # assume following sorting of 2nd derivs in hm:
  # hm[,1]: 2nd deriv of first parameter
  # hm[,2]: 2nd deriv of second parameter
  # hm[,3]: cross-deriv of 2nd and 3rd param
  #
  # sorting of gradient gr:
  # gr[,1]: 1st deriv of first param
  # gr[,2]: 1st deriv of 2nd param
  
  hessX <- hessMat1(hm[,1], gr[,1], X, etaX, dlinkX, ddlinkX, weights = weights)
  if(!fix){
    hessZ <- hessMat1(hm = hm[,2], gr = gr[,2], X = Z, eta = etaZ, dlinkX = dlinkZ,
                      ddlinkX = ddlinkZ, weights = weights)
    hessXZ <- hessCross(crossDeriv = hm[,3], X = X, Z = Z, etaX, etaZ, dlinkX,
                        dlinkZ, weights = weights)
    hessZX <- t(hessXZ)
  }
  cbind(rbind(hessX, if(fix) NULL else hessZX), if(fix) NULL else rbind(hessXZ, hessZ))
}

# cross derivatives in hessians
hessCross <- function(crossDeriv, X, Z, etaX, etaZ, dlinkmuX, dlinkmuZ, weights = 1) {
  hc <- weights * crossDeriv * dlinkmuX(etaX) * dlinkmuZ(etaZ)
  crossprod(hc * X, Z)
}
