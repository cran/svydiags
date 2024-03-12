svycollinear <- function (mod, intcpt = TRUE, w, Vcov, sc = TRUE, svyglm.obj,
          rnd = 3, fuzz = 0.3)
{
  if (is.matrix(mod) & svyglm.obj==TRUE){
    stop("If mod is a matrix, e.g., from model.matrix, svyglm.obj must be set to FALSE. \n")
  }
  if (!(is.matrix(mod) || is.data.frame(mod))) {
    X <- as.data.frame(mod$model)
  }
  fac <- vector("logical", length = ncol(X))
  for (k in 1:ncol(X)) {
    fac[k] <- is.factor(X[, k])
  }
  if (any(fac)){
    stop("Factors need to be converted to columns of 0's and 1's. Use model.matrix(.) to create X matrix and input as mod=X.\n")
  }

  if (is.matrix(mod) || is.data.frame(mod)) {
    X <- as.matrix(mod)
  }
  else {
    X <- as.matrix(mod$model)
  }

  if (svyglm.obj) {
    X <- X[, -1]
    if (intcpt) {
      X <- cbind(rep(1, nrow(X)), X)
      colnames(X)[1] <- "Intercept"
    }
    X <- X[, 1:(ncol(X) - 1)]
  }
  rtwX <- sqrt(w) * X
  if (sc == TRUE) {
    X2 <- sqrt(apply(rtwX * rtwX, 2, sum))
    scwX <- rtwX
    scwX <- scwX %*% diag(1/X2)
  }
  UDV <- svd(scwX)
  condI <- max(UDV$d)/UDV$d
  G <- t(rtwX) %*% rtwX %*% Vcov
  UDV <- svd(rtwX)
  Q <- (UDV$v %*% diag((1/UDV$d)^2)) * t(t(UDV$v) %*% G)
  Qbar <- rowSums(Q)
  Pi <- t(Q) %*% diag(1/Qbar)
  out <- cbind(condI, Pi)
  colnames(out) <- c("Cond Index", colnames(X))
  out <- round(out, rnd)
  out[out <= fuzz] <- "."
  return(data.frame(out))
}
