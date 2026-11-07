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
#' @param ... arguments to be passed on to \code{age2ratio}
#' @return a list with the intercept (\code{a}), slope (\code{b}) and
#'     other data that can be passed on to \code{scatterplot}.
#' @noRd
LRisochron <- function(x,...){ UseMethod("LRisochron",x) }
#' @noRd
#' @export
LRisochron.default <- function(x,left=FALSE,hide=NULL,omit=NULL,
                               anchor=0,inverse=FALSE,...){
    y0rat <- gety0rat(x)
    DPrat <- getDPrat(x)
    op <- iratio(y0rat)
    x2calc <- clear(x,hide=hide,omit=omit)
    init <- init_LRisochron(x2calc,left=left,anchor=anchor,y0rat=y0rat)
    fit <- optim(init$par,fn=LL_LRisochron,
                 x=x2calc,gam=init$gam,y0=init$y0,y0rat=y0rat,
                 left=left,hessian=TRUE)
    iratio(y0rat,x=op[1],e=op[2])
    E <- inverthess(fit$hessian)
    covmat <- matrix(0,4,4)
    J <- matrix(0,2,4)
    if (!is.null(init$gam)){
        a <- exp(fit$par[3])
        b <- age2ratio(tt=exp(init$gam),ratio=DPrat,...)[1]
        covmat[2:4,2:4] <- E
        J[1,4] <- a
    } else if (!is.null(init$y0)){
        tt <- exp(fit$par[1])
        a <- init$y0
        b <- age2ratio(tt=tt,ratio=DPrat,...)[1]
        covmat[1:3,1:3] <- E
        J[2,1] <- dDPdt(x=x,tt=tt,...)*tt
    } else {
        tt <- exp(fit$par[1])
        a <- exp(fit$par[4])
        b <- age2ratio(tt=tt,ratio=DPrat,...)[1]
        covmat <- E
        J[1,4] <- a
        J[2,1] <- dDPdt(x=x,tt=tt,...)*tt
    }
    E <- J %*% covmat %*% t(J)
    out <- list(a=c(a,sqrt(E[1,1])),
                b=c(b,sqrt(E[2,2])),
                cov.ab=E[1,2])
    if (inverse) out <- invertfit(out)
    out$xyz <- data2york(x,inverse=inverse)
    out$model <- 4
    out$n <- length(x2calc)
    out
}
#' @noRd
#' @export
LRisochron.ArAr <- function(x,left=FALSE,hide=NULL,omit=NULL,
                            anchor=0,inverse=FALSE,...){
    LRisochron.default(x,left=left,hide=hide,omit=omit,
                       anchor=anchor,inverse=inverse,J=x$J[1],...)
}
#' @noRd
#' @export
LRisochron.ThU <- function(x,left=FALSE,hide=NULL,omit=NULL,
                            anchor=0,inverse=FALSE,...){
    LRisochron.default(x,left=left,hide=hide,omit=omit,
                       anchor=anchor,inverse=inverse,...)
}

#' @noRd
dDPdt <- function(x,...){ UseMethod("dDPdt",x) }
#' @noRd
#' @export
dDPdt.default <- function(x,...){ stop("No default method implemented.") }
#' @noRd
#' @export
dDPdt.ArAr <- function(x,tt,J,...){
    L <- lambda('K40')[1]
    exp(L*tt)*L/J
}
#' @noRd
#' @export
dDPdt.ThU <- function(x,tt,...){
    L <- lambda('Th230')[1]
    exp(-L*tt)*L
}


init_LRisochron <- function(x,left=FALSE,anchor=0,y0rat){
    tst <- age(x,i2i=FALSE,exterr=FALSE)
    sgn <- ifelse(left,-1,1)
    zs <- cbind(sgn*log(tst[,1]),tst[,2]/tst[,1])
    fit <- min_age_model(zs=zs)
    y0 <- gam <- NULL
    gam <- NULL
    if (anchor[1]==1){
        p <- fit$par
        y0 <- iratio(y0rat)[1]
    } else if (anchor[1]==2 & length(anchor)>1){
        p <- c(fit$par[-1],log(iratio(y0rat)[1]))
        gam <- log(anchor[2])
    } else {
        p <- c(fit$par,log(iratio(y0rat)[1]))
    }
    list(par=p,y0=y0,gam=gam)
}

LL_LRisochron <- function(p,x,left=FALSE,gam=NULL,y0=NULL,y0rat){
    if (!is.null(gam)){
        p <- c(gam,p)
    }
    if (is.null(y0)){
        y0 <- exp(tail(p,1))
    }
    gps <- min_age_model_mappar(p[1:3])
    iratio(y0rat,y0)
    tst <- age(x,i2i=FALSE,exterr=FALSE)
    sgn <- ifelse(left,-1,1)
    if (any(tst[,1]<0)){
        out <- Inf
    } else {
        zs <- cbind(sgn*log(tst[,1]),tst[,2]/tst[,1])
        out <- get_minage_L(c(gps,gps[1]),zs)
    }
    out
}
