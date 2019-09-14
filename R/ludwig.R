#' @title
#' Linear regression of U-Pb data with correlated errors, taking
#' into account decay constant uncertainties.
#'
#' @description
#' Implements the maximum likelihood algorithm for Total-Pb/U isochron
#' regression of Ludwig (1998)
#'
#' @details
#' The 3-dimensional regression algorithm of Ludwig and Titterington
#' (1994) was modified by Ludwig (1998) to fit so-called `Total Pb-U
#' isochrons'. These are constrained to a radiogenic endmember
#' composition that falls on the \code{\link{concordia}} line. In its
#' most sophisticated form, this algorithm does not only allow for
#' correlated errors between variables, but also between
#' aliquots. \code{IsoplotR} currently uses this algorithm to
#' propagate decay constant uncertainties in the total Pb-U isochron
#' ages.
#'
#' @param x an object of class \code{UPb}
#' 
#' @param alpha cutoff value for confidence intervals
#' 
#' @param ... optional
#' 
#' @return
#' \describe{
#'
#' \item{par}{a vector with the lower concordia intercept,
#' the common Pb ratios and the dispersion parameter.}
#'
#' \item{cov}{the covariance matrix of \code{par}}
#'
#' \item{df}{the degrees of freedom of the model fit (\eqn{3n-3},
#' where \eqn{n} is the number of aliquots).}
#'
#' \item{mswd}{the mean square of weighted deviates (a.k.a. reduced
#' Chi-square statistic) for the fit.}
#'
#' \item{p.value}{p-value of a Chi-square test for the linear fit}
#'
#' }
#'
#' @examples
#' f <- system.file("UPb4.csv",package="IsoplotR")
#' d <- read.data(f,method="U-Pb",format=4)
#' fit <- ludwig(d)
#' @references
#' Ludwig, K.R., 1998. On the treatment of concordant uranium-lead
#' ages. Geochimica et Cosmochimica Acta, 62(4), pp.665-676.
#'
#' Ludwig, K.R. and Titterington, D.M., 1994. Calculation of
#' \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#' Cosmochimica Acta, 58(22), pp.5031-5042.
#'
#' @export
ludwig <- function(x,...){ UseMethod("ludwig",x) }
#' @rdname ludwig
#' @export
ludwig.default <- function(x,...){
    stop( "No default method available (yet)." )
}
#' @param exterr propagate external sources of
#' uncertainty (i.e. decay constants)?
#' @param model one of three regression models:
#'
#' \code{1}: fit a discordia line through the data using the maximum
#' likelihood algorithm of Ludwig (1998), which assumes that the
#' scatter of the data is solely due to the analytical
#' uncertainties. In this case, \code{IsoplotR} will either calculate
#' an upper and lower intercept age (for Wetherill concordia), or a
#' lower intercept age and common
#' (\eqn{^{207}}Pb/\eqn{^{206}}Pb)\eqn{_\circ}-ratio intercept (for
#' Tera-Wasserburg). If \eqn{MSWD}>0, then the analytical
#' uncertainties are augmented by a factor \eqn{\sqrt{MSWD}}.
#'
#' \code{2}: fit a discordia line ignoring the analytical uncertainties
#'
#' \code{3}: fit a discordia line using a modified maximum likelihood
#' algorithm that includes accounts for any overdispersion by adding a
#' geological (co)variance term.
#' @param anchor control parameters to fix the intercept age or common
#'     Pb composition of the discordia fit. This is a two-element
#'     list.
#'
#' The first element is a boolean flag indicating whether the
#' discordia line should be anchored. If this is \code{FALSE}, then
#' the second item is ignored and both the common Pb composition and
#' age are estimated.
#'
#' If the first element is \code{TRUE} and the second element is
#' \code{NA}, then the common Pb composition is fixed at the values
#' stored in \code{settings('iratio',...)}.
#'
#' If the first element is \code{TRUE} and the second element is a
#' number, then the discordia line is forced to intersect the
#' concordia line at an age equal to that number.
#'
#' @seealso \code{\link{concordia}}, \code{\link{titterington}},
#'     \code{\link{isochron}}
#' @rdname ludwig
#' @export
ludwig.UPb <- function(x,exterr=FALSE,alpha=0.05,model=1,anchor=list(FALSE,NA)){
    out <- get.ta0b0w(x,exterr=exterr,model=model,anchor=anchor)
    if (model==3) out$wci <- profile_LL_ludwig(fit=out,x=x,alpha=alpha)
    out$n <- length(x)
    mswd <- mswd.lud(out$par,x=x,anchor=anchor)
    out <- c(out,mswd)
    if (x$format %in% c(1,2,3)) parnames <- c('t','76i','w')
    else if (x$format %in% c(4,5,6)) parnames <- c('t','64i','74i','w')
    else if (x$format %in% c(7,8)) parnames <- c('t','68i','78i','w')
    else stop("Illegal input format.")
    names(out$par) <- parnames
    rownames(out$cov) <- parnames
    colnames(out$cov) <- parnames
    out
}

mswd.lud <- function(ta0b0,x,anchor=list(FALSE,NA)){
    ns <- length(x)
    tt <- ta0b0[1]
    a0 <- ta0b0[2]
    out <- list()
    anchored <- anchor[[1]]
    tanchored <- is.numeric(anchor[[2]])
    if (x$format<4){
        b0 <- 0
        if (anchored) out$df <- ns-1
        else out$df <- ns-2
    } else {
        b0 <- ta0b0[3]
        if (anchored && tanchored) out$df <- 2*ns-1
        else if (tanchored) out$df <- 2*ns-2
        else out$df <- 2*ns-3
    }
    SS <- data2ludwig(x,tt=tt,a0=a0,b0=b0,w=0)$SS
    if (out$df>0){
        out$mswd <- as.vector(SS/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    } else {
        out$mswd <- 1
        out$p.value <- 1
    }
    out
}

get.ta0b0w <- function(x,exterr=FALSE,model=1,anchor=list(FALSE,NA),w=NA,...){
    out <- list()
    ta0b0 <- get.ta0b0.model2(x,anchor=anchor)
    if (model==2){
        out$par <- c(ta0b0,0)
        out$cov <- fisher.model2(x,fit=fit,anchor=anchor)
        out$cov <- matrix(0,length(ta0b0)+1,length(ta0b0)+1)
    } else {
        lower <- c(ta0b0/2,0)
        upper <- c(2*ta0b0,ta0b0[1])
        if (model==1){
            init <- c(ta0b0,0)  # no overdispersion
        } else {
            if (is.numeric(w)) init <- c(ta0b0,w)
            else init <- c(ta0b0,ta0b0[1]/100) # 1% overdispersion as a first guess
        }
        names(init)[3] <- 'w'
        fit <- optifix(parms=init,fn=LL.lud.UPb,gr=LL.lud.UPb.gr,
                       method="L-BFGS-B",x=x,exterr=exterr,
                       fixed=fixit(x,anchor=anchor,model=model,w=w),
                       lower=lower,upper=upper,
                       control=list(fnscale=-1),...)
        out <- list()
        out$LL <- fit$value
        out$par <- fit$par
        out$cov <- fisher.lud(x,fit=fit,anchor=anchor)
    }
    out$model <- model
    out$exterr <- exterr
    out
}

get.ta0b0.model2 <- function(x,anchor=list(FALSE,NA)){
    if (x$format<4){
        np <- 2
        init <- c(1,0)
        lower <- c(0,0)
        upper <- c(10000,100)
    } else {
        np <- 3
        init <- c(1,0,0)
        lower <- c(0,0,0)
        upper <- c(10000,100,100)
    }
    fit <- optifix(parms=init,fn=SS.model2,method="L-BFGS-B",
                   x=x,exterr=exterr,
                   fixed=fixit(x,anchor=anchor,model=model),
                   lower=lower,upper=upper,hessian=TRUE,...)
    out <- list()
    out$ta0b0 <- fit$par
    out$cov <- np*fit$minimum/(length(np)-2)*solve(fit$hessian)
    out
}

SS.model2 <- function(ta0b0,x){
    tt <- ta0b0[1]
    a0 <- ta0b0[2]
    out <- list()
    if (x$format<4){
        xy <- data2york(x,option=2)[,c('X','Y'),drop=FALSE]
        xr <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        yr <- age_to_Pb207Pb206_ratio(tt,st=0,d=x$d)[1]
        yp <- yr+(a0-yr)*(xr-xy[,'X'])/xr
        out <- sum((yp-xy[,'Y'])^2)
    } else {
        b0 <- ta0b0[2]
        ns <- length(x)
        xy <- matrix(0,ns,4)
        if (x$format<7){
            for (i in 1:ns){
                xy[i,] <- get.UPb.isochron.ratios.204(x,i)$x
            }
        } else {
            for (i in 1:ns){
                xy[i,] <- get.UPb.isochron.ratios.208(x,i,tt=tt)$x[1:4]
            }
        }
        x6 <- xy[,1] # U238Pb206
        y6 <- xy[,2] # Pb204Pb206 or Pb208cPb206
        x7 <- xy[,3] # U235Pb207
        y7 <- xy[,4] # Pb204Pb207 or Pb208cPb207
        r86 <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        r57 <- age_to_U235Pb207_ratio(tt,st=0,d=x$d)[1]
        out <- list()
        out$y0 <- c(0,0) # (6/48)i, (7/48)i
        out$y0[1] <- a0
        yp <- (x6-xr6)*a0/xr6
        SS06 <- sum((yp-y6)^2)
        out$y0[2] <- b0
        yp <- (x7-xr7)*b0/xr7
        SS07 <- sum((yp-y7)^2)
        out <- SS06 + SS07
    }
    out
}

get.ta0b0.model2 <- function(x,anchor=list(FALSE,NA)){
    if (x$format<4) ta0b0 <- get.ta0b0.model2.2D(x,anchor=anchor)
    else ta0b0 <- get.ta0b0.model2.3D(x,anchor=anchor)
    ta0b0
}
get.ta0b0.model2.2D <- function(x,anchor=list(FALSE,NA)){
    tlim <- c(1e-5,min(get.Pb206U238.age(x)[,1]))
    if (!anchor[[1]]){
        tt <- stats::optimise(SS.model2.2D,interval=tlim,x=x)$minimum
        fits <- model2fit.2D(tt,x=x)
        a0 <- fits$y0
    } else if (is.na(anchor[[2]])){
        a0 <- settings('iratio','Pb207Pb206')[1]
        tt <- stats::optimise(SS.model2.2D,interval=tlim,x=x,a0=a0)$minimum
    } else if (is.numeric(anchor[[2]])){
        tt <- anchor[[2]]
        fits <- model2fit.2D(tt,x=x)
        a0 <- fits$y0[1]
    }
    out <- c(tt,a0)
    names(out) <- c('t','76i')
    out
}
get.ta0b0.model2.3D <- function(x,anchor=list(FALSE,NA)){
    tlim <- c(1e-5,max(get.Pb206U238.age(x)[,1]))
    if (!anchor[[1]]){
        tt <- stats::optimise(SS.model2.3D,interval=tlim,x=x)$minimum
        fits <- model2fit.3D(tt,x=x)
        a0 <- fits$y0[1]
        b0 <- fits$y0[2]
    } else if (is.na(anchor[[2]])){
        if (x$format%in%c(4,5,6)){
            a0 <- settings('iratio','Pb206Pb204')[1]
            b0 <- settings('iratio','Pb207Pb204')[1]
        } else {
            a0 <- 1/settings('iratio','Pb208Pb206')[1]
            b0 <- 1/settings('iratio','Pb208Pb207')[1]
        }
        tt <- stats::optimise(SS.model2.3D,interval=tlim,x=x,a0=a0,b0=b0)$minimum
    } else if (is.numeric(anchor[[2]])){
        tt <- anchor[[2]]
        fits <- model2fit.3D(tt,x=x)
        a0 <- fits$y0[1]
        b0 <- fits$y0[2]
    }
    out <- c(tt,a0,b0)
    if (x$format%in%c(4,5,6))
        names(out) <- c('t','64i','74i')
    else
        names(out) <- c('t','68i','78i')
    out
}
SS.model2.2D <- function(tt,x,a0=NA){
    model2fit.2D(tt=tt,x=x,a0=a0)$SS
}
SS.model2.3D <- function(tt,x,a0=NA,b0=NA){
    model2fit.3D(tt=tt,x=x,a0=a0,b0=b0)$SS
}
model2fit.2D <- function(tt,x=x,a0=NA,b0=NA){
    ns <- length(x)
    xy <- data2york(x,option=2)[,c('X','Y'),drop=FALSE]
    xr <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
    yr <- age_to_Pb207Pb206_ratio(tt,st=0,d=x$d)[1]
    out <- list()
    out$y0 <- 0
    if (is.na(a0)){
        fit <- stats::lm(I(xy[,'Y']-yr) ~ I(xy[,'X']-xr) + 0)
        out$y0 <- yr - fit$coef[1]*xr
        out$SS <- sum(fit$residuals^2)
    } else {
        out$y0 <- a0
        yp <- yr+(a0-yr)*(xr-xy[,'X'])/xr
        out$SS <- sum((yp-xy[,'Y'])^2)
    }
    out
}
model2fit.3D <- function(tt,x=x,a0=NA,b0=NA){
    ns <- length(x)
    xy <- matrix(0,ns,4)
    for (i in 1:ns){
        if (x$format %in% c(4,5,6))
            xy[i,] <- get.UPb.isochron.ratios.204(x,i)$x
        else
            xy[i,] <- get.UPb.isochron.ratios.208(x,i,tt=tt)$x[1:4]
    }
    x6 <- xy[,1] # U238Pb206
    y6 <- xy[,2] # Pb204Pb206 or Pb208cPb206
    xr6 <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
    out <- list()
    out$y0 <- c(0,0) # (6/48)i, (7/48)i
    if (is.na(a0)){
        fit06 <- stats::lm(y6 ~ I(x6-xr6) + 0)
        out$y0[1] <- -1/(xr6*fit06$coef)
        SS06 <- sum(fit06$residuals^2)
    } else {
        out$y0[1] <- a0
        yp <- (x6-xr6)*a0/xr6
        SS06 <- sum((yp-y6)^2)
    }
    x7 <- xy[,3] # U235Pb207
    y7 <- xy[,4] # Pb204Pb207 or Pb208cPb207
    xr7 <- age_to_U235Pb207_ratio(tt,st=0,d=x$d)[1]
    if (is.na(b0)){
        fit07 <- stats::lm(y7 ~ I(x7-xr7) + 0)
        out$y0[2] <- -1/(xr7*fit07$coef)
        SS07 <- sum(fit07$residuals^2)
    } else {
        out$y0[2] <- b0
        yp <- (x7-xr7)*b0/xr7
        SS07 <- sum((yp-y7)^2)
    }
    out$SS <- SS06 + SS07
    out
}

LL.lud.UPb <- function(ta0b0w,x,exterr=FALSE,LL=TRUE){
    tt <- ta0b0w[1]
    a0 <- ta0b0w[2]
    if (x$format<4){
        b0 <- 0
        w <- ta0b0w[3]
    } else {
        b0 <- ta0b0w[3]
        w <- ta0b0w[4]        
    }
    l <- data2ludwig(x,tt=tt,a0=a0,b0=b0,w=w,exterr=exterr)
    if (LL) out <- l$LL
    else out <- l$SS
    out
}

LL.lud.UPb.gr <- function(ta0b0w,x,exterr=FALSE){
    tt <- ta0b0w[1]
    a0 <- ta0b0w[2]
    if (x$format<4){
        b0 <- 0
        w <- ta0b0w[3]
    } else {
        b0 <- ta0b0w[3]
        w <- ta0b0w[4]
    }
    data2ludwig(x,tt=tt,a0=a0,b0=b0,w=w,exterr=exterr,jacobian=TRUE)$jacobian
}

fisher.lud <- function(x,fit,exterr=TRUE,anchor=list(FALSE,NA),...){
    ns <- length(x)
    i1 <- 1:ns
    if (x$format %in% c(1,2,3)){
        fish <- data2ludwig_2D(x,tt=fit$par[1],a0=fit$par[2],w=fit$par[3],
                               exterr=exterr,hessian=TRUE)$hessian
        i2 <- (ns+1):(ns+3)
    } else if (x$format %in% c(4,5,6)){
        fish <- data2ludwig_3D(x,tt=fit$par[1],a0=fit$par[2],b0=fit$par[3],
                               w=fit$par[4],exterr=exterr,hessian=TRUE)$hessian
        i2 <- (ns+1):(ns+4)
    } else if (x$format %in% c(7,8)){
        fish <- data2ludwig_Th(x,tt=fit$par[1],a0=fit$par[2],b0=fit$par[3],
                               w=fit$par[4],exterr=exterr,hessian=TRUE)$hessian
        i2 <- (ns+1):(ns+4)
    }
    anchorfish(AA=fish[i1,i1],BB=fish[i1,i2],
               CC=fish[i2,i1],DD=fish[i2,i2],anchor=anchor)
}
anchorfish <- function(AA,BB,CC,DD,anchor=list(FALSE,NA)){
    npar <- nrow(DD)
    out <- matrix(0,npar,npar)
    if (anchor[[1]] && is.numeric(anchor[[2]])){
        bb <- subset(BB,select=-1)
        cc <- t(subset(t(CC),select=-1))
        dd <- DD[-1,-1]
        out[2:npar,2:npar] <- blockinverse(AA,bb,cc,dd)
    } else if (anchor[[1]]){
        bb <- subset(BB,select=1)
        cc <- t(subset(t(CC),select=1))
        dd <- DD[1,1]
        out[1,1] <- blockinverse(AA,bb,cc,dd)
    } else {
        out <- blockinverse(AA,BB,CC,DD)
    }
    out
}

data2ludwig <- function(x,tt,a0,b0=0,w=0,exterr=FALSE,
                        jacobian=FALSE,hessian=FALSE){
    if (x$format %in% c(1,2,3))
        out <- data2ludwig_2D(x,tt=tt,a0=a0,w=w,exterr=exterr,
                              jacobian=jacobian,hessian=hessian)
    else if (x$format %in% c(4,5,6))
        out <- data2ludwig_3D(x,tt=tt,a0=a0,b0=b0,w=w,exterr=exterr,
                              jacobian=jacobian,hessian=hessian)
    else if (x$format %in% c(7,8))
        out <- data2ludwig_Th(x,tt=tt,a0=a0,b0=b0,w=w,exterr=exterr,
                              jacobian=jacobian,hessian=hessian)
    else stop('Incorrect input format.')
    out
}
data2ludwig_2D <- function(x,tt,a0,w=0,exterr=FALSE,
                           jacobian=FALSE,hessian=FALSE){
    out <- list()
    U <- iratio('U238U235')[1]
    ns <- length(x)
    zeros <- rep(0,ns)
    E <- matrix(0,2*ns+6,2*ns+6)
    J <- matrix(0,2*ns,2*ns+6)
    J[1:(2*ns),1:(2*ns)] <- diag(2*ns)
    X <- zeros
    Y <- zeros
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    for (i in 1:ns){
        wd <- wetherill(x,i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        E[c(i,ns+i),c(i,ns+i)] <- wd$cov
        J[i,2*ns+2] <- -D$dPb207U235dl35     # dKdl35
        J[i,2*ns+4] <- -D$dPb207U235dl31     # dKdl31
        J[ns+i,2*ns+1] <- -D$dPb206U238dl38  # dLdl31
        J[ns+i,2*ns+3] <- -D$dPb206U238dl34  # dLdl34
        J[ns+i,2*ns+5] <- -D$dPb206U238dl30  # dLdl30
        J[ns+i,2*ns+6] <- -D$dPb206U238dl26  # dLdl26
    }
    E[2*ns+1,2*ns+1] <- lambda('U238')[2]^2
    E[2*ns+2,2*ns+2] <- lambda('U235')[2]^2
    E[2*ns+3,2*ns+3] <- lambda('U234')[2]^2
    E[2*ns+4,2*ns+4] <- lambda('Pa231')[2]^2
    E[2*ns+5,2*ns+5] <- lambda('Th230')[2]^2
    E[2*ns+6,2*ns+6] <- lambda('Ra226')[2]^2
    Ew <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D)
    ED <- J%*%E%*%t(J) + Ew
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    O <- blockinverse(AA=ED[i1,i1],BB=ED[i1,i2],
                      CC=ED[i2,i1],DD=ED[i2,i2],doall=TRUE)
    K0 <- X - D$Pb207U235 + a0*U*(D$Pb206U238 - Y)
    A <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*a0*U + K0%*%(O[i1,i2]+t(O[i2,i1])))
    B <- -(a0*U*(O[i1,i1]+t(O[i1,i1]))*a0*U + (O[i2,i2]+t(O[i2,i2])) +
           a0*U*(O[i1,i2]+t(O[i1,i2])) + (O[i2,i1]+t(O[i2,i1]))*a0*U)
    L <- as.vector(solve(B,A))
    c0 <- Y - D$Pb206U238 - L
    K <- X - D$Pb207U235 - a0*U*c0
    KL <- c(K,L)
    out$SS <- KL%*%O%*%KL
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- -(2*ns*log(2*pi) + detED + out$SS)/2
    if (jacobian | hessian){
        JKL <- matrix(0,2*ns,ns+2) # derivatives of KL w.r.t. c0, t, and a0
        colnames(JKL) <- c(paste0('c0[',i1,']'),'t','a0')
        rownames(JKL) <- c(paste0('K[',i1,']'),paste0('L[',i1,']'))
        diag(JKL[i1,i1]) <- -U*a0       # dKdc0
        diag(JKL[i2,i1]) <- -1          # dLdc0
        JKL[i1,'t']  <- -D$dPb207U235dt # dKdt
        JKL[i2,'t']  <- -D$dPb206U238dt # dLdt
        JKL[i1,'a0'] <- -U*c0           # dKda0
        out$jacobian <- rep(0,3) # derivatives of LL w.r.t t, a0 and w
        names(out$jacobian) <- c('t','a0','w')
        out$jacobian['t'] <- -KL%*%O%*%JKL[,'t']
        out$jacobian['a0'] <- -KL%*%O%*%JKL[,'a0']
        dEDdw <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D,deriv=1)
        dlnDetEDdw <- trace(O%*%dEDdw)
        dOdw <- -O%*%dEDdw%*%O
        out$jacobian['w'] <- -(dlnDetEDdw + KL%*%dOdw%*%KL)/2
    }
    if (hessian){
        out$hessian <- matrix(0,ns+3,ns+3)
        hesnames <- c(paste0('c0[',1:ns,']'),'t','a0','w')
        rownames(out$hessian) <- hesnames
        colnames(out$hessian) <- hesnames
        d2KLdt2 <- rep(0,2*ns)
        d2KLdt2[i1] <- -D$d2Pb207U235dt2 # d2Kdt2
        d2KLdt2[i2] <- -D$d2Pb206U238dt2 # d2Ldt2
        d2KLdc0da0 <- matrix(0,2*ns,ns)
        diag(d2KLdc0da0[i1,i1]) <- -U    # d2Kdc0da0
        out$hessian[i1,i1] <- t(JKL[,i1])%*%O%*%JKL[,i1]              # d2dc02
        out$hessian['t',i1] <- t(JKL[,'t'])%*%O%*%JKL[,i1]            # d2dc0dt
        out$hessian['a0',i1] <-
            t(JKL[,'a0'])%*%O%*%JKL[,i1] + KL%*%O%*%d2KLdc0da0        # d2dc0da0
        out$hessian['t','t'] <-
            t(JKL[,'t'])%*%O%*%JKL[,'t'] + KL%*%O%*%d2KLdt2           # d2dt2
        out$hessian['a0','a0'] <- t(JKL[,'a0'])%*%O%*%JKL[,'a0']      # d2da02
        d2EDdw2 <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D,deriv=2)
        d2lnDetEDdw2 <- trace(O%*%d2EDdw2 - trace(O%*%dEDdw%*%O%*%dEDdw))
        d2Odw2 <- dOdw%*%dEDdw%*%O + O%*%d2EDdw2%*%O + O%*%dEDdw%*%dOdw 
        out$hessian['w','w'] <- -(d2lnDetEDdw2 + KL%*%d2Odw2%*%KL )/2 # d2dw2
        out$hessian['w',i1] <- KL%*%dOdw%*%JKL[,i1]                   # d2dwdc0
        out$hessian['w','t'] <- KL%*%dOdw%*%JKL[,'t']                 # d2dwdt
        out$hessian['w','a0'] <- KL%*%dOdw%*%JKL[,'a0']               # d2dwda0
        out$hessian[i1,'t'] <- out$hessian['t',i1]                    # d2dc0dt
        out$hessian[i1,'a0'] <- out$hessian['a0',i1]                  # d2dc0da0
        out$hessian[i1,'w'] <- out$hessian['w',i1]                    # d2dc0dw
        out$hessian['t','w'] <- out$hessian['w','t']                  # d2dtdw
        out$hessian['a0','w'] <- out$hessian['w','a0']                # d2da0dw
    }
    out
}
# rederived from Ludwig (1998):
data2ludwig_3D <- function(x,tt,a0,b0,w=0,exterr=FALSE,
                           jacobian=FALSE,hessian=FALSE){
    out <- list()
    U <- iratio('U238U235')[1]
    ns <- length(x)
    zeros <- rep(0,ns)
    E <- matrix(0,3*ns+6,3*ns+6)
    J <- matrix(0,3*ns,3*ns+6)
    J[1:(3*ns),1:(3*ns)] <- diag(3*ns)
    X <- zeros
    Y <- zeros
    Z <- zeros
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    x75 <- D$Pb207U235
    x68 <- D$Pb206U238
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        Z[i] <- wd$x['Pb204U238']
        E[c(i,ns+i,2*ns+i),c(i,ns+i,2*ns+i)] <- wd$cov
        J[i,3*ns+2] <- -D$dPb207U235dl35     #dKdl35
        J[i,3*ns+4] <- -D$dPb207U235dl31     #dKdl31
        J[ns+i,3*ns+1] <- -D$dPb206U238dl38  #dLdl38
        J[ns+i,3*ns+3] <- -D$dPb206U238dl34  #dLdl34
        J[ns+i,3*ns+5] <- -D$dPb206U238dl30  #dLdl30
        J[ns+i,3*ns+6] <- -D$dPb206U238dl26  #dLdl26
    }
    E[3*ns+1,3*ns+1] <- lambda('U238')[2]^2
    E[3*ns+2,3*ns+2] <- lambda('U235')[2]^2
    E[3*ns+3,3*ns+3] <- (lambda('U234')[2]*1000)^2
    E[3*ns+4,3*ns+4] <- (lambda('Pa231')[2]*1000)^2
    E[3*ns+5,3*ns+5] <- (lambda('Th230')[2]*1000)^2
    E[3*ns+6,3*ns+6] <- (lambda('Ra226')[2]*1000)^2
    Ew <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D)
    ED <- J%*%E%*%t(J) + Ew
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    i3 <- (2*ns+1):(3*ns)
    O <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                         DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                         GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    o11 <- O[i1,i1]; o12 <- O[i1,i2]; o13 <- O[i1,i3]
    o21 <- O[i2,i1]; o22 <- O[i2,i2]; o23 <- O[i2,i3]
    o31 <- O[i3,i1]; o32 <- O[i3,i2]; o33 <- O[i3,i3]
    K0 <- X - U*b0*Z - D$Pb207U235
    L0 <- Y - a0*Z - D$Pb206U238
    V <- t(K0%*%(o11+t(o11))*U*b0 + L0%*%(o12+t(o21))*U*b0 + K0%*%(o12+t(o21))*a0 +
           L0%*%(o22+t(o22))*a0 + K0%*%(o13+t(o31)) + L0%*%(o23+t(o32)))
    W <- -(U*b0*(o11+t(o11))*U*b0 + U*b0*(o12+t(o12))*a0 + U*b0*(o13+t(o13)) +
           a0*(o21+t(o21))*U*b0 + a0*(o22+t(o22))*a0 + a0*(o23+t(o23)) +
           (o31+t(o31))*U*b0 + (o32+t(o32))*a0 + (o33+t(o33)))
    M <- as.vector(solve(W,V))
    c0 <- as.vector(Z - M)
    K <- as.vector(X - U*b0*c0 - x75)
    L <- as.vector(Y - a0*c0 - x68)
    KLM <- c(K,L,M)
    out$SS <- KLM%*%O%*%KLM
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- -(3*log(2*pi)/2 + detED/2 + out$SS)/2
    if (jacobian | hessian){
        JKLM <- matrix(0,3*ns,ns+3) # derivatives of KLM w.r.t. c0, t, a0 and b0
        colnames(JKLM) <- c(paste0('c0[',i1,']'),'t','a0','b0')
        diag(JKLM[i1,i1]) <- -U*b0     # dKdc0
        diag(JKLM[i2,i1]) <- -a0       # dLdc0
        diag(JKLM[i3,i1]) <- -1        # dMdc0
        JD[i1,'t']  <- -D$dPb207U235dt # dKdt
        JD[i2,'t']  <- -D$dPb206U238dt # dLdt
        JD[i2,'a0'] <- -c0             # dLda0
        JD[i1,'b0'] <- -U*c0           # dKdb0
        out$jacobian <- rep(0,4) # derivatives of LL w.r.t t, a0, b0 and w
        names(out$jacobian) <- c('t','a0','b0','w')
        out$jacobian['t'] <- -KLM%*%O%*%JKLM[,'t']
        out$jacobian['a0'] <- -KLM%*%O%*%JKLM[,'a0']
        out$jacobian['b0'] <- -KLM%*%O%*%JKLM[,'b0']
        dEDdw <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D,deriv=1)
        dlnDetEDdw <- trace(O%*%dEDdw)
        dOdw <- -O%*%dEDdw%*%O
        out$jacobian['w'] <- -(dlnDetEDdw + KLM%*%dOdw%*%KLM)/2
    }
    if (hessian){
        d2Kdt2 <- rep(-D$d2Pb207U235dt2,ns)
        d2Ldt2 <- rep(-D$d2Pb206U238dt2,ns)
        d2Ldc0da0 <- rep(-1,ns)
        d2Kdc0db0 <- rep(-U,ns)
    }
    out
}
data2ludwig_Th <- function(x,tt,a0,b0,w=0,exterr=FALSE,
                           jacobian=FALSE,hessian=FALSE){
    out <- list()
    l2 <- lambda('Th232')[1]
    U <- iratio('U238U235')[1]
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    x75 <- D$Pb207U235
    x68 <- D$Pb206U238
    x82 <- exp(l2*tt)-1
    ns <- length(x)
    zeros <- rep(0,ns)
    X <- zeros
    Y <- zeros
    Z <- zeros
    W <- zeros
    K0 <- zeros
    L0 <- zeros
    E <- matrix(0,4*ns+7,4*ns+7)
    J <- matrix(0,3*ns,4*ns+7)
    J[1:(3*ns),1:(3*ns)] <- diag(3*ns)
    for (i in 1:ns){
        wd <- wetherill(x,i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        Z[i] <- wd$x['Pb208Th232']
        W[i] <- wd$x['Th232U238']
        K0[i] <- X[i] - (Z[i]-x82)*b0*U*W[i] - x75
        L0[i] <- Y[i] - (Z[i]-x82)*a0*W[i] - x68
        E[(0:3)*ns+i,(0:3)*ns+i] <- wd$cov
        J[i,4*ns+2] <- -D$dPb207U235dl35      # dKdl35
        J[i,4*ns+5] <- -D$dPb207U235dl31      # dKdl31
        J[ns+i,4*ns+1] <- -D$dPb206U238dl38   # dLdl38
        J[ns+i,4*ns+3] <- -D$dPb206U238dl34   # dLdl34
        J[ns+i,4*ns+6] <- -D$dPb206U238dl30   # dLdl30
        J[ns+i,4*ns+7] <- -D$dPb206U238dl26   # dLdl26
        J[2*ns+i,4*ns+4] <- -tt*exp(l2*tt)    # dMdl32
    }
    E[4*ns+1,4*ns+1] <- lambda('U238')[2]^2
    E[4*ns+2,4*ns+2] <- lambda('U235')[2]^2
    E[4*ns+3,4*ns+3] <- (lambda('U234')[2]*1000)^2
    E[4*ns+4,4*ns+4] <- lambda('Th232')[2]^2
    E[4*ns+5,4*ns+5] <- (lambda('Pa231')[2]*1000)^2
    E[4*ns+6,4*ns+6] <- (lambda('Th230')[2]*1000)^2
    E[4*ns+7,4*ns+7] <- (lambda('Ra226')[2]*1000)^2
    Ew <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D)
    ED <- J%*%E%*%t(J) + Ew
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    i3 <- (2*ns+1):(3*ns)
    O <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                         DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                         GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    o11 <- O[i1,i1]; o12 <- O[i1,i2]; o13 <- O[i1,i3]
    o21 <- O[i2,i1]; o22 <- O[i2,i2]; o23 <- O[i2,i3]
    o31 <- O[i3,i1]; o32 <- O[i3,i2]; o33 <- O[i3,i3]
    Wd <- diag(W)
    AA <- (Wd%*%o11%*%Wd)*(U*b0)^2 + (Wd%*%o22%*%Wd)*a0^2 + o33 +
        U*a0*b0*Wd%*%(o12+o21)%*%Wd +
        U*b0*(Wd%*%o13+o31%*%Wd) + a0*(Wd%*%o23+o32%*%Wd)
    BT <- t(U*b0*K0%*%o11%*%Wd + a0*L0%*%o22%*%Wd +
            a0*K0%*%o12%*%Wd + U*b0*L0%*%o21%*%Wd +
            K0%*%o13 + L0%*%o23)
    CC <- U*b0*Wd%*%o11%*%K0 + a0*Wd%*%o22%*%L0 +
        a0*Wd%*%o12%*%K0 + U*b0*Wd%*%o21%*%L0 +
        o13%*%K0 + o23%*%L0
    M <- as.vector(solve(-(AA+t(AA)),(BT+CC)))
    c0 <- as.vector(Z - M - x82)
    K <- as.vector(X - c0*b0*U*W - x75)
    L <- as.vector(Y - c0*a0*W - x68)
    KLM <- c(K,L,M)
    out$SS <- KLM%*%O%*%KLM
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- 3*log(2*pi)/2 - detED/2 - out$SS/2
    if (jacobian | hessian){
        JD <- matrix(0,3*ns,4)
        jacnames <- c('t','a0','b0','c0','w')
        colnames(JD) <- jacnames[1:4]
        JD[i1,'t']  <- -D$dPb207U235dt # dKdt
        JD[i2,'t']  <- -D$dPb206U238dt # dLdt
        JD[i3,'t']  <- -exp(l2*tt)*l2  # dMdt
        JD[i2,'a0'] <- -c0*W           # dLda0
        JD[i1,'b0'] <- -c0*U*W         # dKdb0
        JD[i1,'c0'] <- -b0*U*W         # dKdc0
        JD[i2,'c0'] <- -a0*W           # dLdc0
        JD[i3,'c0'] <- -1              # dMdc0
        out$jacobian <- rep(0,5)
        names(out$jacobian) <- jacnames
        for (pn in jacnames[1:4]){
            out$jacobian[pn] <- -KLM%*%O%*%JD[,pn]
        }
        dEDdw <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D,deriv=1)
        out$jacobian['w'] <- (KLM%*%O%*%dEDdw%*%O%*%KLM - trace(O%*%dEDdw))/2
    }
    if (hessian){
        d2Kdt2 <- rep(-D$d2Pb207U235dt2,ns)
        d2Ldt2 <- rep(-D$d2Pb206U238dt2,ns)
        d2Mdt2 <- rep(-exp(l2*tt)*l2^2,ns)
        d2Ldc0da0 <- -W
        d2Kdc0db0 <- -U*W
    }
    out
}

get.Ew <- function(w=0,format=1,ns=1,tt=0,D=mclean(),deriv=0){
    if (format<4) ndim <- 2
    else ndim <- 3
    J <- matrix(0,ndim,1)
    J[1,1] <- -D$dPb207U235dt # dKdt
    J[2,1] <- -D$dPb206U238dt # dLdt
    if (format>6) J[3,1] <- -D$dPb208Th232dt # dMdt
    if (deriv==1) {
        dEdw <- 2*w
        Ew <- J%*%dEdw%*%t(J)
    } else if (deriv==2) {
        d2Edw2 <- 2
        Ew <- J%*%d2Edw2%*%t(J)
    } else {
        E <- w^2
        Ew <- J%*%E%*%t(J)
    }
    out <- matrix(0,ndim*ns,ndim*ns)
    diag(out[1:ns,1:ns]) <- Ew[1,1]
    diag(out[(ns+1):(2*ns),(ns+1):(2*ns)]) <- Ew[2,2]
    diag(out[1:ns,(ns+1):(2*ns)]) <- Ew[1,2]
    diag(out[(ns+1):(2*ns),1:ns]) <- Ew[2,1]
    if (format>4){
        diag(out[(ns+1):(2*ns),(2*ns+1):(3*ns)]) <- Ew[3,3]
        diag(out[1:ns,(2*ns+1):(3*ns)]) <- Ew[1,3]
        diag(out[(2*ns+1):(3*ns),1:ns]) <- Ew[3,1]
        diag(out[(ns+1):(2*ns),(2*ns+1):(3*ns)]) <- Ew[2,3]
        diag(out[(2*ns+1):(3*ns),(ns+1):(2*ns)]) <- Ew[3,2]
    }
    out
}

fixit <- function(x,anchor=list(FALSE,NA),model=1,w=NA){
    if (x$format<4) np <- 3
    else np <- 4
    out <- rep(FALSE,np)
    if (model==1 | is.numeric(w)) out[np] <- TRUE # fix w
    if (anchor[[1]]){ # anchor t or a0(,b0)
        if (is.numeric(anchor[[2]])) out[1] <- TRUE # fix t
        else out[2:(np-1)] <- TRUE # fix a0(,b0)
    }
    out
}
