svystdres <- function (mobj, stvar = NULL, clvar = NULL, doplot = FALSE)
{
    p <- length(mobj$coefficients)
    e <- mobj$residuals
    mrows <- names(mobj$y)
#    mdata <- mobj$data[rownames(mobj$data) == mrows, ]
    mdata <- mobj$data[rownames(mobj$data) %in% mrows, ]
    sighat <- rhohat <- 0
    mbar <- 1
    if (is.null(clvar)) {
        if (!is.null(stvar)) {
            nh <- table(mdata[, stvar])
            n <- sum(nh)
            vst <- (nh - 1) * as.vector(by(e, mdata[, stvar],
                var))
            sighat <- sum(vst)/(n - p)
        }
        else {
            n <- length(e)
            sighat <- (n - 1) * var(e)/(n - p)
        }
    }

    if (!is.null(clvar)) {
        if (!is.null(stvar)) {
            new.psu <- paste(mdata[, stvar], mdata[, clvar],
                sep = ".")
       #     new.psu <- as.numeric(new.psu)
            new.str <- mdata[order(new.psu), stvar]
        }
        else {
            new.psu <- mdata[, clvar]
            new.psu <- as.numeric(new.psu)
        }
        n <- length(unique(new.psu))
        mhi <- table((new.psu))
        m <- sum(mhi)
        mbar <- mean(mhi)
        if (!is.null(stvar)) {
            e <- e[order(new.psu)]
           first <- do.call("rbind", list(by(mdata[,stvar],paste(mdata[,stvar],mdata[,clvar], sep="."), head,1)))
            nh <- table (first)
            eclvar <- by(e, INDICES = new.psu, var)
            eclvar[is.na(eclvar)] <- min(eclvar, na.rm = TRUE)
            Phat <- sum(eclvar)/n
            ebarhi <- by(e, new.psu, mean)
            ebarh <- by(e, new.str, mean)
            ebarh <- rep(ebarh, nh)
            Qhat <- sum(mhi * (ebarhi - ebarh)^2)/(n - 1)
            Dhat <- (m - sum(mhi^2)/m)/(n - 1)
        }
        else {
            e <- e[order(new.psu)]
            eclvar <- by(e, INDICES = new.psu, var)
            eclvar[is.na(eclvar)] <- min(eclvar, na.rm = TRUE)
            Phat <- sum(eclvar)/n
            ebari <- by(e, new.psu, mean)
            ebar <- mean(e)
            Qhat <- sum(mhi * (ebari - ebar)^2)/(n - 1)
            Dhat <- (m - sum(mhi^2)/m)/(n - 1)
        }
        sighat <- Phat + (Qhat - Phat)/Dhat
        rhohat <- (Qhat - Phat)/Dhat/sighat
    }
    stdresids <- e/sqrt(sighat)
    names(stdresids) <- names(e)
    if (doplot) {
        plot(1:length(e), stdresids, xlab = "sample element")
        abline(h = c(-3, 0, 3), col = "red", lwd = 2)
    }
    list(stdresids = stdresids, n = n, mbar = mbar, rtsighat = sqrt(sighat),
        rhohat = rhohat)
}
