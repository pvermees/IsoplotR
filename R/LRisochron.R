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
#' @param anchor see \code{isochron}
#' @param inverse see \code{isochron}
#' @param inverse toggles between normal and inverse isochrons. See
#'     \code{\link{isochron}}.
#' @param ... arguments to be passed on to \code{age2ratio}
#' @return a list with the intercept (\code{a}), slope (\code{b}) and
#'     other data that can be passed on to \code{scatterplot}. If
#'     \code{x} has class \code{ThU}, also returns \code{age} and
#'     \code{y0}.
#' @noRd
LRisochron <- function(x,...){ UseMethod("LRisochron",x) }
#' @noRd
LRisochron.default <- function(x,left=FALSE,hide=NULL,omit=NULL,
                               anchor=0,inverse=FALSE,...){
    y0rat <- gety0rat(x)
    x2calc <- clear(x,hide=hide,omit=omit)
    init <- init_LRisochron(x2calc,left=left,
                            anchor=anchor,y0rat=y0rat)
    fit <- stats::optim(init$par,fn=LL_LRisochron,
                        x=x2calc,gam=init$gam,
                        y0=init$y0,y0rat=y0rat[1],
                        left=left,hessian=TRUE)
    out <- fit2LRisochron(x=x,fit=fit,gam=init$gam,
                          y0=init$y0,left=left,
                          inverse=inverse,...)
    out$xyz <- data2york(x,inverse=inverse)
    out
}
#' @noRd
LRisochron.ArAr <- function(x,left=FALSE,hide=NULL,omit=NULL,
                            anchor=0,inverse=FALSE,...){
    y0rat <- gety0rat(x)
    op <- iratio(y0rat)
    out <- LRisochron.default(x,left=left,hide=hide,omit=omit,
                              anchor=anchor,inverse=inverse,
                              J=x$J[1],...)
    iratio(y0rat,x=op[1],e=op[2])
    out
}
#' @noRd
LRisochron.ThU <- function(x,left=FALSE,hide=NULL,omit=NULL,
                           anchor=0,inverse=FALSE,...){
    U8Th2 <- x$U8Th2
    out <- LRisochron.default(x,left=left,hide=hide,omit=omit,
                              anchor=anchor,inverse=inverse,...)
    x$U8Th2 <- U8Th2
    out
}

#' Derivative of the daughter-parent ratio w.r.t time
#' @noRd
dDPdt <- function(x,...){ UseMethod("dDPdt",x) }
#' @noRd
dDPdt.default <- function(x,...){ stop("No default method implemented.") }
#' @param J the J-factor
#' @noRd
dDPdt.ArAr <- function(x,tt,J,...){
    L <- lambda('K40')[1]
    exp(L*tt)*L/J
}
#' @noRd
dDPdt.ThU <- function(x,tt,...){
    L <- lambda('Th230')[1]
    exp(-L*tt)*L
}

#' Initialise the LRisochron parameters
#' @noRd
init_LRisochron <- function(x,...){ UseMethod("init_LRisochron",x) }
#' @noRd
init_LRisochron.default <- function(x,left=FALSE,anchor=0,y0rat,...){
    if (anchor[1]==1){
        y0 <- iratio(y0rat)[1]
    } else {
        y0 <- NULL
    }
    gam <- NULL
    tst <- age(x)
    sgn <- ifelse(left,-1,1)
    zs <- cbind(sgn*log(tst[,1]),tst[,2]/tst[,1])
    fit <- min_age_model(zs=zs)
    if (anchor[1]==1){
        p <- fit$par
    } else if (anchor[1]==2 & length(anchor)>1){
        p <- c(fit$par[-1],log(iratio(y0rat)[1]))
        gam <- sgn*log(anchor[2])
    } else {
        p <- c(fit$par,log(iratio(y0rat)[1]))
    }
    list(par=p,y0=y0,gam=gam)
}
#' @noRd
init_LRisochron.ThU <- function(x,left=FALSE,anchor=0,y0rat,...){
    if (anchor[1]==1){
        y0 <- 1/x$U8Th2
        tst <- age(x,Th0i=2)
    } else {
        y0 <- NULL
        tst <- age(x)
    }
    gam <- NULL
    sgn <- ifelse(left,-1,1)
    zs <- cbind(sgn*log(tst[,1]),tst[,2]/tst[,1])
    fit <- min_age_model(zs=zs)
    if (anchor[1]==1){
        p <- fit$par
    } else if (anchor[1]==2 & length(anchor)>1){
        p <- c(fit$par[-1],-log(x$U8Th2))
        gam <- sgn*log(anchor[2])
    } else {
        p <- c(fit$par,0)
    }
    list(par=p,y0=y0,gam=gam)
}

LL_LRisochron <- function(p,x,left=FALSE,gam=NULL,
                          y0=NULL,y0rat){
    if (!is.null(gam)){
        p <- c(gam,p)
    }
    if (is.null(y0)){
        y0 <- exp(utils::tail(p,1))
    }
    sgn <- ifelse(left,-1,1)
    gps <- min_age_model_mappar(p[1:3])
    tst <- age4LRisochron(x,y0rat=y0rat,y0=y0)
    if (any(tst[,1]<0)){
        out <- Inf
    } else {
        zs <- cbind(sgn*log(tst[,1]),tst[,2]/tst[,1])
        out <- get_minage_L(c(gps,gps[1]),zs)
    }
    out
}

#' Calculate single grain ages given the non-radiogenic ratio \code{y0}
#' @noRd
age4LRisochron <- function(x,...){ UseMethod("age4LRisochron",x) }
#' @noRd
age4LRisochron.default <- function(x,y0rat,y0,...){
    iratio(y0rat,y0)
    age(x,i2i=FALSE)
}
#' @noRd
age4LRisochron.ThU <- function(x,y0,...){
    x$U8Th2 <- y0
    age(x,Th0i=2)
}

#' Error propagation for the LRisochron fit
#' @noRd
fit2LRisochron <- function(x,...){ UseMethod("fit2LRisochron",x) }
#' @noRd
fit2LRisochron.default <- function(x,fit,gam=NULL,y0=NULL,
                                   inverse=FALSE,left=FALSE,...){
    DPrat <- getDPrat(x)
    E <- inverthess(fit$hessian)
    covmat <- matrix(0,4,4)
    J <- matrix(0,2,4)
    sgn <- ifelse(left,-1,1)
    if (!is.null(gam)){
        a <- exp(fit$par[3])
        b <- age2ratio(tt=exp(sgn*gam),ratio=DPrat,...)[1]
        covmat[2:4,2:4] <- E
        J[1,4] <- a
    } else if (!is.null(y0)){
        tt <- exp(sgn*fit$par[1])
        a <- y0
        b <- age2ratio(tt=tt,ratio=DPrat,...)[1]
        covmat[1:3,1:3] <- E
        J[2,1] <- dDPdt(x=x,tt=tt,...)*tt
    } else {
        tt <- exp(sgn*fit$par[1])
        a <- exp(fit$par[4])
        b <- age2ratio(tt=tt,ratio=DPrat,...)[1]
        covmat <- E
        J[1,4] <- a
        J[2,1] <- dDPdt(x=x,tt=tt,...)*tt
    }
    fit2LRisochron_helper(x=x,covmat=covmat,J=J,
                          a=a,b=b,inverse=inverse)
}
#' @noRd
fit2LRisochron.ThU <- function(x,fit,gam=NULL,y0=NULL,
                               inverse=FALSE,left=FALSE,...){
    DPrat <- "Th230U238"
    E <- inverthess(fit$hessian)
    covmat <- matrix(0,4,4)
    J <- matrix(0,2,4)
    sgn <- ifelse(left,-1,1)
    if (!is.null(gam)){
        b <- age2ratio(tt=exp(sgn*gam),ratio=DPrat,...)[1]
        a <- exp(fit$par[3])*(1-b)
        covmat[2:4,2:4] <- E
        J[1,4] <- a
    } else if (!is.null(y0)){
        dtdg <- tt <- exp(sgn*fit$par[1])
        b <- age2ratio(tt=tt,ratio=DPrat,...)[1]
        a <- y0*(1-b)
        dbdt <- dDPdt(x=x,tt=tt,...)
        covmat[1:3,1:3] <- E
        J[1,1] <- -y0*dbdt*dtdg
        J[2,1] <- dbdt*dtdg
    } else {
        dtdg <- tt <- exp(sgn*fit$par[1])
        y0 <- exp(fit$par[4])
        b <- age2ratio(tt=tt,ratio=DPrat,...)[1]
        a <- y0*(1-b)
        dbdg <- dDPdt(x=x,tt=tt,...)
        covmat <- E
        J[1,4] <- a
        J[2,1] <- dbdg*dtdg
    }
    fit2LRisochron_helper(x=x,covmat=covmat,J=J,a=a,b=b,
                          inverse=inverse,type="d")
}
fit2LRisochron_helper <- function(x,covmat,J,a,b,
                                  inverse=FALSE,type="p"){
    E <- J %*% covmat %*% t(J)
    out <- list(a=c(a,sqrt(E[1,1])),
                b=c(b,sqrt(E[2,2])),
                cov.ab=E[1,2])
    if (inverse) out <- invertfit(out,type=type)
    out$model <- 4
    out$n <- length(x)
    out
}
