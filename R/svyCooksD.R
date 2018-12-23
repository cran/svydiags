svyCooksD <- function(mobj, stvar=NULL, clvar=NULL, doplot=FALSE){
#browser()
    xdat <- mobj$data[, all.vars(mobj$formula)]
    X <- model.matrix(mobj$formula, xdat)
    w <- mobj$weights
    n <- length(w)
    p <- length(mobj$coefficients)

    rho <- 0
    mbar <- 1

    XWX <- (t(X) * w) %*% X
    Ainv <- ginv(XWX)
#    Ainv <- chol2inv(chol(XWX))
    h <- svyhat(mobj)

    num <- Ainv %*% t(X) * w * (1-h)
    vinv <- ginv(vcov(mobj))
#    vinv <- chol2inv(chol(vcov(mobj)))
    cook <- diag(t(num) %*% vinv %*% num)
    if (!is.null(clvar)){
        tmp <- svystdres(mobj, stvar, clvar)
        rho <- tmp$rhohat
        n <- tmp$n
        mbar <- tmp$mbar
    }
    adj <- n * mbar * (1 + rho*(mbar-1)) / p
    mcook <- sqrt(adj * cook)
    if (doplot){
        plot(1:length(mcook), mcook, xlab="sample element")
        abline(h = c(2,3), col="red", lwd=2)
    }
    mcook
}
