svyvif <- function(X, w, V){
    if (
        dim(X)[1] != length(w) | dim(X)[1] != dim(V)[1]
    )
    stop("Dimensions of X, w, or V do not match. Check that input contains complete cases only. \n")


    rtw <- sqrt(w)
    rtwX <- rtw * X
    rtwVrtw <- t(rtw * t(rtw * V))
    xkbar <- apply(X, 2, function(col) weighted.mean(x=col, w=w))

    hat.N <- sum(w)
    nam <- dimnames(X)[[2]]

    vif.fit <- zeta <- varrho <- reg.vif <- rep(0, ncol(X))
    names(vif.fit)  <-  names(reg.vif) <- names(zeta) <- names(varrho) <- nam
    X.fit <- data.frame(1,X)    # Add intercept back to the X matrix

    for (k in 1:ncol(X)){
      X.reg <- rtw * X.fit[, -(k+1)]
      xkw <- rtwX[,k]
      xk <- X[,k]
      xk.bar <- sum(xk * w)/sum(w)

      fit <- glm.fit(x = X.reg, y = xkw)
      e.xk <- fit$residuals

      zeta[k] <- (t(e.xk) %*% rtwVrtw %*% e.xk)/(t(e.xk) %*% e.xk)
      varrho[k] <- (t(xkw) %*% xkw - hat.N * xk.bar^2)/
                    (t(xkw - rtw * xk.bar) %*% rtwVrtw %*% (xkw - rtw * xk.bar))
      reg.vif[k] <- (t(xkw) %*% xkw - hat.N * xk.bar^2)/(t(e.xk) %*% e.xk)
      vif.fit[k] <-  zeta[k] * varrho[k] * reg.vif[k]
    }

    zeta.varrho <- zeta * varrho
    v <- data.frame("svy.vif" = vif.fit,
                    "reg.vif" = reg.vif,
                    "zeta" = zeta,
                    "varrho" = varrho,
                    "zeta.x.varrho" = zeta.varrho)
    v
}
