#' Leftmost and rightmost isochrons
#'
#' Modifies the Galbraith and Laslett's minimum age module to fit
#' overdispersed Pb-Pb and Th-U isochron data.
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
    init <- init_LRisochron(x2calc,a=a,b=b,left=left)
    out <- stats::optim(par=init,fn=get_LRisochron_L,
                        yd=x2calc,a=a,b=b,left=left,
                        hessian=TRUE)
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
#' If \code{anchor[1]=1}: fix the non-radiogenic composition at the
#' values stored in \code{settings('iratio',...)}, OR, if \code{x} has
#' class \code{ThU} and \code{x$format} = \code{3} or \code{4}, fix
#' the intercept at the value stored in \code{x$U8Th2}.
#'
#' If \code{anchor[1]=2}: fix the age at the value stored in \code{anchor[2]}.
#' @noRd
LRisochron.PbPb <- function(x,inverse=TRUE,anchor=0,hide=NULL,omit=NULL,...){
    yd <- data2york(x,inverse=TRUE)
    yd2calc <- clear(yd,hide,omit)
    if (anchor[1]<1){
        out <- LRisochron.default(yd2calc,left=TRUE)
    } else if (anchor[1]==1){
        Pb74 <- iratio('Pb207Pb204')[1]
        Pb64 <- iratio('Pb206Pb204')[1]
        y0 <- Pb74/Pb64
        offset <- 1/Pb64
        yd2calc[,'X'] <- yd2calc[,'X'] - offset
        stop("TODO")
    } else if (anchor[1]==2 & length(anchor)>1){
        y0 <- age2ratio(tt=anchor[2],ratio="Pb207Pb206")[1]
        stop("TODO")
    } else {
        stop("Invalid anchor.")
    }
    if (inverse){
        out$xyz <- yd
    } else {
        out <- invertfit(out,type='d')
        out$xyz <- normal2inverse(yd,type='d')
    }
    out
}
#' @noRd
LRisochron.ThU <- function(x,inverse=TRUE,anchor=0,hide=NULL,omit=NULL,...){
    if (x$format<3){
        stop("Rightmost isochrons are only available for ThU formats 3 and 4.")
    }
    yd <- data2york(x,inverse=FALSE)
    if (anchor[1]<1){
        out <- LRisochron.default(yd,left=FALSE,hide=hide,omit=omit)
    } else {
        yd2calc <- clear(yd,hide,omit)
        lpi <- 0
        y0i <- gi <- y0 <- b <- NULL # note: y0 = Th2U8i
        if (anchor[1]==1){
            y0i <- y0 <- 1/x$U8Th2
            g <- (yd2calc[,'Y']-y0i)/(yd2calc[,'X']-y0i)
            lsi <- log(stats::sd(g))
            gi <- mean(g)
            stop("TODO")
        } else if (anchor[1]==2 & length(anchor)>1){
            b <- age2ratio(tt=anchor[2],ratio='Th230U238')[1]
            leftmost <- which.min(yd[,'X'])
            x2 <- yd2calc[leftmost,'X']
            y2 <- yd2calc[leftmost,'Y']
            y0i <- y2 - b * x2
            g <- (yd2calc[,'Y']-y0i)/(yd2calc[,'X']-y0i)
            lsi <- log(stats::sd(g))
            stop("TODO")
        } else {
            stop("Invalid anchor")
        }
    }
    if (inverse){
        out <- invertfit(out,type='d')
        out$xyz <- normal2inverse(yd,type='d')
    } else {
        out$xyz <- yd
    }
    out
}

init_LRisochron <- function(yd,a=NULL,b=NULL,left=FALSE){
    X <- yd[,1]
    Y <- yd[,3]
    lpi <- 0
    if (is.null(a) && is.null(b)){
        anchor <- 0
        h <- chull(X,Y)
        nh <- length(h)
        vertices <- c(h,h[1])
        gr <- diff(Y[vertices])/diff(X[vertices])
        misfit <- Inf
        for (i in which(gr>0)){
            b <- gr[i]
            a <- Y[vertices[i]] - b * X[vertices[i]]
            if (a>0){
                if (left) lsi <- log(stats::sd((Y-a)/X))
                else lsi <- log(stats::sd(X/(Y-a)))
                init <- c(lpi,lsi,log(a),log(b))
                fit <- stats::optim(par=init,fn=get_LRisochron_L,yd=yd,left=left)
                if (fit$value < misfit){
                    misfit <- fit$value
                    out <- fit$par
                }
            }
        }
    } else if (is.null(b)){
        anchor <- 'a'
        fit <- stats::lm(I(Y - a) ~ X - 1)
        lbi <- log(abs(fit$coefficients[2]))
        lsi <- log(stats::sd((Y-a)/X))
        out <- c(lpi,lsi,lbi)
    } else if (is.null(a)){
        anchor <- 'b'
        fit <- stats::lm(Y ~ 1 + offset(b * X))
        ai <- fit$coefficients[1]
        lsi <- log(stats::sd((Y-ai)/X))
        out <- c(lpi,lsi,ai)
    } else {
        lsi <- log(stats::sd((Y-a)/X))
        out <- c(lpi,lsi)
    }
    out
}

get_LRisochron_L <- function(pars,yd,a=NULL,b=NULL,left=FALSE,...){
    mappar <- function(pars,b=NULL,left=FALSE){
        prop <- 1/(exp(pars[1])+1)
        sig <- exp(pars[2])
        if (is.null(b)) b <- exp(tail(pars,n=1))
        gam <- ifelse(left,1/b,b)
        mu <- gam
        c(gam,prop,sig,mu)
    }
    get_zs <- function(yd,a=0,left=FALSE){
        z <- (yd[,'Y']-a)/yd[,'X']
        sz <- sqrt(errorprop1x2(J1=-z/yd[,'X'],
                                J2=1/yd[,'X'],
                                E11=yd[,'sX']^2,
                                E22=yd[,'sY']^2,
                                E12=yd[,'rXY']*yd[,'sX']*yd[,'sY']))
        if (left){
            sz <- sz/z^2
            z <- 1/z
        }   
        cbind(z,sz)
    }
    minage_pars <- mappar(pars,b=b,left=left)
    if (is.null(a)) a <- exp(pars[3])
    zs <- get_zs(yd,a=a,left=left)
    get_minage_L(pars=minage_pars,zs=zs)
}
