#' Leftmost and rightmost isochrons
#'
#' Modifies the Galbraith and Laslett's minimum age module to fit
#' overdispersed (Th-U) isochron data.
#' @param x an IsoplotR data object or an \eqn{n \times 5} matrix of
#'     inverse isochron data.
#' @param left logical, switches between leftmost and rightmost
#'     isochrons
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the isochron age calculation.
#' @param inverse toggles between normal and inverse isochrons. See
#'     \code{\link{isochron}}.
#' @param ... unused optional arguments
#' @return a list with the intercept (\code{a}), slope (\code{b}) and
#'     other data that can be passed on to \code{scatterplot}.
#' @noRd
LRisochron <- function(x,...){ UseMethod("LRisochron",x) }
#' @noRd
#' @export
LRisochron.default <- function(x,left=FALSE,
                               hide=NULL,omit=NULL,
                               anchor=0,...){
    d2calc <- clear(x,hide,omit)
    init <- init_LRisochron(d2calc,left=left,anchor=anchor)
    fit <- optim(init,fn=LL_LRisochron,yd=d2calc,
                 anchor=anchor,hessian=TRUE)
    covmat <- matrix(0,4,4)
    if (anchor[1]==1 & length(anchor)>1){
        covmat[1:3,1:3] <- inverthess(fit$hessian)
        x0 <- exp(ifelse(left,1,-1)*fit$par[1])
        y0 <- anchor[2]
    } else if (anchor[1]==2 & length(anchor)>1){
        covmat[2:4,2:4] <- inverthess(fit$hessian)
        x0 <- anchor[2]
        y0 <- exp(fit$par[3])
    } else {
        covmat <- inverthess(fit$hessian)
        x0 <- exp(ifelse(left,1,-1)*fit$par[1])
        y0 <- exp(fit$par[4])
    }
    a <- y0
    b <- -y0/x0
    J <- matrix(0,nrow=2,ncol=4)
    J[1,4] <- y0    # d(a)/d(log[y0])
    J[2,1] <- y0/x0 # d(b)/d(log[x0])
    J[2,4] <- b     # d(b)/d(log[y0])
    E <- J %*% covmat %*% t(J)
    list(a=unname(c(a,sqrt(E[1,1]))),
         b=unname(c(b,sqrt(E[2,2]))),
         cov.ab=unname(E[1,2]),
         model=4,
         xyz=x,
         n=nrow(d2calc))
}
#' @noRd
#' @export
LRisochron.ArAr <- function(x,inverse=TRUE,left=FALSE,
                            hide=NULL,omit=NULL,anchor=0,...){
    if (anchor[1]==1){
        anchor[2] <- 1/iratio('Ar40Ar36')[1]
    } else if (anchor[1]==2 & length(anchor)>1){
        anchor[2] <- 1/age2ratio(tt=anchor[2],
                                 ratio='Ar40Ar39',
                                 J=x$J[1])[1]
    }
    yd <- data2york(x=x,inverse=TRUE)
    fit <- LRisochron(yd,left=left,hide=hide,
                      omit=omit,anchor=anchor)
    if (inverse){
        out <- fit
    } else {
        out <- invertfit(fit)
        out$xyz <- data2york(x=x,inverse=FALSE)
    }
    out
}
#' @noRd
#' @export
LRisochron.ThU <- function(x,inverse=TRUE,anchor=0,
                           hide=NULL,omit=NULL,left=FALSE,...){
    stop('Not implemented yet')
}

init_LRisochron <- function(yd,left=FALSE,anchor=0){
    if (anchor[1]==2 & length(anchor)>1){
        x0 <- anchor[2]
        dy <- max(yd[,3]) - min(yd[,3])
        dx <- x0 - min(yd[,1])
        y0 <- (max(yd[,3]) + min(yd[,1])*dy/dx)
        zs <- y02zs(y0=y0,yd=yd,left=left)
        sigma <- sqrt(mean((zs[,1]-log(x0))^2))
        out <- c(0,log(sigma),log(y0))
    } else if (anchor[1]==1 & length(anchor)>1){
        y0 <- anchor[2]
        zs <- y02zs(y0=y0,yd=yd,left=left)
        fit <- min_age_model(zs=zs)
        out <- fit$par
    } else {
        dy <- max(yd[,3]) - min(yd[,3])
        dx <- max(yd[,1]) - min(yd[,1])
        y0 <- (max(yd[,3]) + min(yd[,1])*dy/dx)
        zs <- y02zs(y0=y0,yd=yd,left=left)
        fit <- min_age_model(zs=zs)
        out <- c(fit$par,log(y0))
    }
    out
}

y02zs <- function(y0,yd,left=FALSE){
    DP <- yd[,'X']*y0/(y0-yd[,'Y'])
    J1 <- y0/(y0-yd[,'Y'])
    J2 <- DP/(y0-yd[,'Y'])
    E11 <- yd[,'sX']^2
    E22 <- yd[,'sY']^2
    E12 <- yd[,'rXY']*yd[,'sX']*yd[,'sY']
    sDP <- sqrt(errorprop1x2(J1=J1,J2=J2,E11=E11,E22=E22,E12=E12))
    pos <- DP>0
    sgn <- ifelse(left,1,-1)
    z <- sgn*rep(-Inf,length(DP))
    z[pos] <- sgn*log(DP[pos])
    s <- sDP/DP
    cbind(z,s)
}

LL_LRisochron <- function(p,yd,left=FALSE,anchor=0){
    if (length(anchor)>1 & length(p)==3){
        if (anchor[1]==1){        # y0
            p <- c(p,log(anchor[2]))
        } else if (anchor[1]==2){ # x0
            p <- c(ifelse(left,1,-1)*log(anchor[2]),p)
        }
    }
    gps <- min_age_model_mappar(p[1:3])
    y0 <- exp(p[4])
    zs <- y02zs(y0=y0,yd=yd,left=left)
    get_minage_L(c(gps,y0),zs)
}
