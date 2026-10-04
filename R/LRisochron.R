#' Leftmost and rightmost isochrons
#'
#' Modifies the Galbraith and Laslett's minimum age module to fit
#' overdispersed (Th-U) isochron data.
#' @param x an IsoplotR data object
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
LRisochron.default <- function(x,left=FALSE,
                               hide=NULL,omit=NULL,
                               a=NULL,b=NULL,...){
    x2calc <- clear(x,hide,omit)
    init <- init_LRisochron(x2calc,a=a,b=b,left=left,...)
    out <- stats::optim(par=init,fn=get_LRisochron_L,
                        yd=x2calc,a=a,b=b,left=left,
                        hessian=TRUE,...)
    covmat <- inverthess(out$hessian)
    np <- length(out$par)
    J <- matrix(0,nrow=2,ncol=np)
    if (is.null(a)){
        a <- exp(out$par[3])
        J[1,3] <- a
    }
    if (is.null(b)){
        b <- exp(out$par[np])
        J[2,np] <- b
    }
    E <- J %*% covmat %*% t(J)
    out$a <- c('a'=unname(a),'s[a]'=unname(sqrt(E[1,1])))
    out$b <- c('b'=unname(b),'s[b]'=unname(sqrt(E[2,2])))
    out$cov.ab <- unname(E[1,2])
    out$xyz <- x
    out$model <- 4
    out$n <- nrow(x2calc)
    out
}
#' @param anchor control parameters to fix the intercept age or
#'     non-radiogenic composition of the isochron fit. This can be a
#'     scalar or a vector.
#'
#' If \code{anchor[1]=0}: do not anchor the isochron.
#'
#' If \code{anchor[1]=1}: fix the intercept at the value stored in \code{x$U8Th2}.
#'
#' If \code{anchor[1]=2}: fix the age at the value stored in \code{anchor[2]}.
#' @noRd
LRisochron.ThU <- function(x,inverse=TRUE,anchor=0,
                           hide=NULL,omit=NULL,left=FALSE,...){
    if (x$format<3){
        stop("Rightmost isochrons are only available for ThU formats 3 and 4.")
    }
    yd <- data2york(x,inverse=FALSE)
    if (anchor[1]==1){
        x0 <- 1/x$U8Th2
        yd2calc <- yd
        yd2calc[,'X'] <- yd[,'X'] - x0
        out <- LRisochron.default(x=yd2calc,omit=omit,hide=hide,
                                  a=x0,left=left)
        out$a[1] <- out$a[1] - x0*out$b[1]
        J <- diag(2)
        J[1,2] <- -x0
        E <- diag(c(out$a[2],out$b[2]))^2
        E[1,2] <- E[2,1] <- out$cov.ab
        covmat <- J %*% E %*% t(J)
        out$a[2] <- sqrt(covmat[1,1])
        out$cov.ab <- covmat[1,2]
    } else if (anchor[1]==2 & length(anchor)>1){
        b <- age2ratio(tt=anchor[2],ratio='Th230U238')[1]
        out <- LRisochron.default(x=yd,b=b,omit=omit,hide=hide,
                                  left=left,ThU=TRUE)
    } else {
        out <- LRisochron.default(x=yd,omit=omit,hide=hide,
                                  left=left,ThU=TRUE)
    }
    if (inverse){
        out <- invertfit(out,type='d')
        out$xyz <- normal2inverse(yd,type='d')
    } else {
        out$xyz <- yd
    }
    out
}

init_lsi <- function(yd,a,b,left=FALSE,ThU=FALSE){
    if (ThU){
        x0 <- y0 <- a/(1-b)
    } else {
        x0 <- 0
        y0 <- a
    }
    zsz <- x0y02zs(yd=yd,x0=x0,y0=y0,left=left)
    log(stats::sd(zsz[,1]))
}

init_LRisochron <- function(yd,a=NULL,b=NULL,left=FALSE,ThU=FALSE){
    X <- yd[,1]
    Y <- yd[,3]
    lpi <- 0
    x0 <- 0
    if (is.null(a) && is.null(b)){
        h <- grDevices::chull(X,Y)
        nh <- length(h)
        vertices <- c(h,h[1])
        gr <- diff(Y[vertices])/diff(X[vertices])
        misfit <- Inf
        for (i in which(gr>0)){
            b <- gr[i]
            a <- Y[vertices[i]] - b * X[vertices[i]]
            if (a>0){
                lsi <- init_lsi(yd=yd,a=a,b=b,left=left,ThU=ThU)
                init <- c(lpi,lsi,log(a),log(b))
                fit <- stats::optim(par=init,fn=get_LRisochron_L,
                                    yd=yd,left=left,ThU=ThU)
                if (fit$value < misfit){
                    misfit <- fit$value
                    out <- fit$par
                }
            }
        }
    } else if (is.null(b)){
        fit <- stats::lm(I(Y - a) ~ X - 1)
        b <- unname(abs(fit$coefficients))
        lbi <- log(b)
        lsi <- init_lsi(yd=yd,a=a,b=b,left=left,ThU=ThU)
        out <- c(lpi,lsi,lbi)
    } else if (is.null(a)){
        fit <- stats::lm(Y ~ 1 + offset(b * X))
        a <- unname(abs(fit$coefficients))
        lsi <- init_lsi(yd=yd,a=a,b=b,left=left,ThU=ThU)
        out <- c(lpi,lsi,log(a))
    } else {
        lsi <- init_lsi(yd=yd,a=a,b=b,left=left,ThU=ThU)
        out <- c(lpi,lsi)
    }
    out
}

x0y02zs <- function(yd,x0=0,y0=0,left=FALSE){
    z <- (yd[,'Y']-y0)/(yd[,'X']-x0)
    sz <- sqrt(errorprop1x2(J1=-z/(yd[,'X']-x0),
                            J2=1/(yd[,'X']-x0),
                            E11=yd[,'sX']^2,
                            E22=yd[,'sY']^2,
                            E12=yd[,'rXY']*yd[,'sX']*yd[,'sY']))
    if (left){
        sz <- sz/z^2
        z <- 1/z
    }   
    cbind(z,sz)
}

get_LRisochron_L <- function(pars,yd,
                             a=NULL,b=NULL,
                             left=FALSE,ThU=FALSE){
    mappar <- function(pars,b=NULL,left=FALSE){
        prop <- logit(pars[1],inverse=TRUE)
        sig <- exp(pars[2])
        if (is.null(b)) b <- exp(utils::tail(pars,n=1))
        gam <- ifelse(left,1/b,b)
        mu <- gam
        c(gam,prop,sig,mu)
    }
    minage_pars <- mappar(pars,b=b,left=left)
    if (is.null(a)){
        a <- exp(pars[3])
    }
    if (is.null(b)){
        gam <- minage_pars[1]
        b <- ifelse(left,1/gam,gam)
    }
    if (ThU){
        x0 <- y0 <- a/(1-b)
    } else {
        x0 <- 0
        y0 <- a
    }
    zs <- x0y02zs(yd=yd,x0=x0,y0=y0,left=left)
    get_minage_L(pars=minage_pars,zs=zs)
}
