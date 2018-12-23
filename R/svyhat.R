svyhat <- function(mobj, doplot=FALSE){
    xdat <- mobj$data[, all.vars(mobj$formula)]
    X <- model.matrix(mobj$formula, xdat)
    w <- mobj$weights
    X <- X * sqrt(w)
    if (is.qr(X))
        n <- nrow(X$qr)
    else {
        n <- nrow(X)
        X <- qr(X)
    }
    h <- rowSums(qr.qy(X, diag(1, nrow = n, ncol = X$rank))^2)
    names(h) <- names(mobj$weights)
    if (doplot){
        plot(h, type= "h")
        abline(h = 3*mean(h), col="red")
    }
    h
}
