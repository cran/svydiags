svycollinear <- function(
                    mobj,
                    X = NULL,
                    w,
                    sc = TRUE,
                    rnd = 3, fuzz = 0.05)
{
  fam <- mobj$family$family
  if (dim(X)[1] != length(w))
      stop("Dimensions of X and w do not match. Check that input contains complete cases only. \n")
  chk.fam <- mobj$family$family %in% c("binomial","gaussian","poisson","quasibinomial","quasipoisson",
                       "gamma","inverse.gaussian")
  if (!chk.fam) stop("Only binomial, gaussian, poisson, quasibinomial, quasipoisson, Gamma, and inverse.gaussian families are supported. \n")

  lnk <- mobj$family$link
  hat.y <- mobj$fitted.values
  beta <- mobj$coefficients

  if (fam == "gaussian"){
      if (lnk != "identity") stop("Only the identity link is supported for the gaussian family.\n")
      if (lnk == "identity") gam <- 1
  }
  if ((fam == "poisson") || (fam == "quasipoisson")){
      if (lnk != "log") stop("Only the log link is supported for the poisson and quasipoisson families.\n")
      if (lnk == "log") gam <- hat.y
  }
  if (fam == "Gamma"){
      if (lnk != "inverse") stop("Only the inverse link is supported for the gamma family.\n")
      if (lnk == "inverse") gam <- hat.y^2
  }
  if (fam == "inverse.gaussian"){
      if (lnk != "1/mu^2") stop("Only the 1/mu^2 link is supported for the inverse.gaussian family.\n")
      if (lnk == "1/mu^2") gam <- hat.y^3/4
  }

  if ((fam == "binomial") || (fam == "quasibinomial")){
      if ( !(lnk=="logit") & !(lnk=="probit") )
            stop("Only logit and probit links are supported for the binomial and quasibinomial families.\n")
      if (lnk == "logit"){
            gam <- hat.y * (1-hat.y)
      }
      if (lnk == "probit"){
            gam <- dnorm(X %*% beta)^2 / (hat.y * (1-hat.y))
      }
  }

   rtwgX <- as.vector(sqrt(w * gam)) * X
   if (sc == TRUE) {
         X2 <- sqrt(apply(rtwgX * rtwgX, 2, sum))
         scwgX <- rtwgX %*% diag(1/X2)
   }
   else {scwgX <- rtwgX}

   UDV <- svd(scwgX)
   condI <- max(UDV$d)/UDV$d

   if (fam == "gaussian"){
     v.beta <- vcov(mobj)
     G <- t(rtwgX) %*% rtwgX %*% v.beta     # G = t(\tilde{X}) %*% \tilde{X} %*% var(\hat{\beta}_{SW})
     Q <- (UDV$v %*% diag((1/UDV$d)^2)) * t(t(UDV$v) %*% G)
   }
   else {
     Q <- t(t(UDV$v * UDV$v)) / UDV$d
   }

   Qdot <- rowSums(Q)
   Pi <- t(Q) %*% diag(1/Qdot)
   out <- cbind(condI, Pi)
   colnames(out) <- c("Cond Index", colnames(X))
   out <- round(out, rnd)
   out[out <= fuzz] <- "."
   return(data.frame(out))
}
