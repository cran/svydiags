Vmat <- function(mobj, stvar = NULL, clvar = NULL){
    e <- residuals(mobj, type ="response")

    mrows <- names(mobj$y)
    mdata <- mobj$data[rownames(mobj$data) %in% mrows, ]
    if (is.null(clvar)) {
        if (!is.null(stvar)) {
            nh <- table(mdata[, stvar])
        }
        else {
            n <- length(e)
        }
    }

        # unstratified, unclustered
    if (is.null(clvar)) {
        if (is.null(stvar)) {
            Vstr <- diag(e * e)
        }
    }

        # stratified, unclustered
    if (is.null(clvar)) {
        if (!is.null(stvar)) {
            nh <- table(mdata[, stvar])
            u.str <- unique(mdata[,stvar])
            for (h in 1:length(u.str)){
                pick <- mdata[ ,stvar] == u.str[h]
                eh <- e[pick]
                Vh <- diag(eh * eh)
                if (h == 1) {Vstr <- Vh}
                else {Vstr <- bdiag(Vstr, Vh)}
            }
        }
    }
        # unstratified, clustered
    if (!is.null(clvar)) {
        if (is.null(stvar)) {
            u.psu <- unique(mdata[, clvar])
            n <- length(u.psu)
                for (i in 1:n){
                    pick <- mdata[, clvar] == u.psu[i]
                    ei <- e[pick]
                    Vi <- (ei - mean(ei)) %*% t(ei - mean(ei))
                    if (i == 1) {V <- Vi}
                    else {V <- bdiag(V, Vi)}
                }
            Vstr <- n/(n-1) * V
        }
    }

        # stratified, clustered
    if (!is.null(clvar)) {
        if (!is.null(stvar)) {
            new.psu <- paste(mdata[, stvar], mdata[, clvar],
                sep = ".")
            new.str <- mdata[order(new.psu), stvar]
            first <- do.call("rbind", list(by(mdata[,stvar],paste(mdata[,stvar],mdata[,clvar], sep="."),
                        head,1)))
            nh <- table(first)
        }
    }

    if (!is.null(clvar)) {
        if (!is.null(stvar)) {
            u.str <- unique(new.str)
            for (h in 1:length(u.str)){
                pick <- mdata[ ,stvar] == u.str[h]
                eh <- e[pick]
                str.psu <- unique(new.psu[pick])
                for (i in 1:length(str.psu)){
                    ehi <- e[new.psu == str.psu[i]]
                    Vhi <- (ehi - mean(ehi)) %*% t(ehi - mean(ehi))
                    if (i == 1) {Vh <- Vhi}
                    else {Vh <- bdiag(Vh, Vhi)}
                }
            Vstr.h <- nh[h]/(nh[h]-1) * Vh
            if (h == 1) {Vstr <- Vstr.h}
            else {Vstr <- bdiag(Vstr, Vstr.h)}
            }
        }
    }
    return(as.matrix(Vstr))
}
