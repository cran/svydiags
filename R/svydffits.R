svydffits <- function(mobj, stvar=NULL, clvar=NULL, z=3){
    xdat <- mobj$data[, all.vars(mobj$formula)]
    X <- model.matrix(mobj$formula, xdat)
    w <- mobj$weights
    n <- length(w)
    p <- length(mobj$coefficients)
    e <- mobj$residuals
    cutoff <- z * sqrt(p/n)

    h <- svyhat(mobj)
    dffit <- e * h / (1-h)
    SEyhat <- SE(predict(mobj, se.fit=TRUE))
    dffits <- dffit / SEyhat

    if (!is.null(clvar)){
        tmp <- svystdres(mobj, stvar, clvar)
        cutoff <- z * sqrt(p / (tmp$n * tmp$mbar * (1 + (tmp$mbar-1)*tmp$rhohat)))
    }

    DFFITS <- list(Dffit = dffit,
                    Dffits = dffits,
                    cutoff = cutoff)
    DFFITS
}
