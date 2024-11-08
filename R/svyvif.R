svyvif <- function(mobj, X, w, stvar=NULL, clvar=NULL){
  if (dim(X)[1] != length(w))
        stop("Dimensions of X and w do not match. Check that input contains complete cases only. \n")
  fam <- mobj$family$family
  lnk <- mobj$family$link
  hat.y <- mobj$fitted.values
  beta <- mobj$coefficients

  chk.fam <- fam %in% c("binomial","gaussian","poisson","quasibinomial","quasipoisson","Gamma",
                                       "inverse.gaussian")
  if (!chk.fam) stop("Only binomial, gaussian, poisson, quasibinomial, and quasipoisson families are supported. \n")

    V <- Vmat(mobj, stvar, clvar)
    nam <- dimnames(X)[[2]]
    vif.fit <- zeta <- varrho <- reg.vif <- Rsq <- rep(0, ncol(X))
    vif.m.fit <- varrho.m <- reg.vif.m <- Rsq.m <- rep(0, ncol(X))
    names(vif.fit)  <-  names(reg.vif) <- names(zeta) <- names(varrho) <- nam
    names(vif.m.fit)  <-  names(reg.vif.m) <- names(varrho.m) <- nam
    X.fit <- data.frame(1,X)    # Add intercept back to the X matrix
    hat.N <- sum(w)

  if (fam == "gaussian"){
    rtw <- sqrt(w)
    rtwX <- rtw * X
    rtwVrtw <- t(rtw * t(rtw * V))
    xkbar <- apply(X, 2, function(col) weighted.mean(x=col, w=w))

    for (k in 1:ncol(X)){
      X.reg <- rtw * X.fit[, -(k+1)]
      xkw <- rtwX[,k]
      xk <- X[,k]
      xkw.bar <- sum(xk * w)/sum(w)

      fit.xk <- glm.fit(x = X.reg, y = xkw)
      e.xkw <- fit.xk$residuals

      zeta[k] <- (t(e.xkw) %*% rtwVrtw %*% e.xkw)/(t(e.xkw) %*% e.xkw)

        # Intcpt adjusted
      varrho.m[k] <- (t(xkw) %*% xkw - hat.N * xkw.bar^2)/
                    (t(xkw - rtw * xkw.bar) %*% rtwVrtw %*% (xkw - rtw * xkw.bar))
      Rsq.m[k] <- 1 - (t(e.xkw) %*% e.xkw) / (t(xkw) %*% xkw - hat.N * xkw.bar^2)
      reg.vif.m[k] <- 1 / (1 - Rsq.m[k])
      vif.m.fit[k] <-  zeta[k] * varrho.m[k] * reg.vif.m[k]

        # No intcpt
      varrho[k] <- (t(xkw) %*% xkw)/ (t(xkw) %*% rtwVrtw) %*% xkw
      Rsq[k] <- 1 - (t(e.xkw) %*% e.xkw) / (t(xkw) %*% xkw)
      reg.vif[k] <- 1 / (1 - Rsq[k])
      vif.fit[k] <-  zeta[k] * varrho[k] * reg.vif[k]
    }
  }

  if (fam != "gaussian"){
    mu <- mobj$fitted.values  # works for all glm's
  if (fam == "binomial" || fam == "quasibinomial") v <- gam <- mu*(1-mu)
  if (fam == "poisson" || fam == "quasipoisson") v <- gam <- mu
  if (fam == "Gamma"){
      if (lnk != "inverse") stop("Only the inverse link is supported for the Gamma family.\n")
      if (lnk == "inverse") v <- gam <- hat.y^2
  }
  if (fam == "inverse.gaussian"){
      if (lnk != "1/mu^2") stop("Only the 1/mu^2 link is supported for the inverse.gaussian family.\n")
      if (lnk == "1/mu^2") {
        gam <- hat.y^3/4
        v <- hat.y^3
      }
  }

    Delta <- 1/v
    rtwg <- sqrt(w * gam)
    rtwgX <- rtwg * X

    XWgX <- t(rtwgX) %*% rtwgX
    VDel <- t(t(V * Delta) * Delta)       # equiv to diag(Delta) %*% V %*% diag(Delta)
    wgam <- w * gam
    hat.N <- sum(wgam)
    XWgamX <- t(rtwgX) %*% rtwgX
    rtwg.Vdel.rtwg <- t(t(VDel * rtwg) * rtwg)

    for (k in 1:ncol(X)){
      Xwg.reg <- rtwg * X.fit[, -(k+1)]
      xkwg <- rtwgX[,k]
      xk <- X[,k]
      X.reg <- X.fit[, -(k+1)]

      xkwg.bar <- sum(xk * wgam)/sum(wgam)

      fit.xkwg <- glm.fit(x = Xwg.reg, y = xkwg)
      e.xkwg <- fit.xkwg$residuals    # contains sqrt(w*gam)

      zeta[k] <- (t(e.xkwg) %*% rtwg.Vdel.rtwg %*% e.xkwg)/(t(e.xkwg) %*% e.xkwg)

      # Intcpt adjusted
      varrho.m[k] <- (t(xkwg) %*% xkwg - hat.N * xkwg.bar^2)/
        (t(xkwg - rtwg * xkwg.bar) %*% rtwg.Vdel.rtwg %*% (xkwg - rtwg * xkwg.bar))
      Rsq.m[k] <- 1 - (t(e.xkwg) %*% e.xkwg) / (t(xkwg) %*% xkwg - hat.N * xkwg.bar^2)
      reg.vif.m[k] <- 1 / (1 - Rsq.m[k])
      vif.m.fit[k] <-  zeta[k] * varrho.m[k] * reg.vif.m[k]

      # No intcpt
      varrho[k] <- (t(xkwg) %*% xkwg)/ (t(xkwg) %*% rtwg.Vdel.rtwg %*% xkwg)
      Rsq[k] <- 1 - (t(e.xkwg) %*% e.xkwg) / (t(xkwg) %*% xkwg)
      reg.vif[k] <- 1 / (1 - Rsq[k])
      vif.fit[k] <-  zeta[k] * varrho[k] * reg.vif[k]
    }
  }
    zeta.varrho.m <- zeta * varrho.m
    zeta.varrho <- zeta * varrho

    vif.m <- data.frame("svy.vif.m" = vif.m.fit,
                    "reg.vif.m" = reg.vif.m,
                    "zeta" = zeta,
                    "varrho.m" = varrho.m,
                    "zeta.x.varrho.m" = zeta.varrho.m,
                    "Rsq.m" = Rsq.m)
    vif <- data.frame("svy.vif" = vif.fit,
                   "reg.vif" = reg.vif,
                   "zeta" = zeta,
                   "varrho" = varrho,
                   "zeta.x.varrho" = zeta.varrho,
                   "Rsq" = Rsq)

  list("Intercept adjusted" = vif.m, "No intercept" = vif)
}
