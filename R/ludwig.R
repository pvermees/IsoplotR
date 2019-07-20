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
#' ages. Future versions of the program will generalise this approach
#' to other chronometers as well.
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
#' \item{par}{a two-element vector with the lower concordia intercept
#' and initial \eqn{^{207}}Pb/\eqn{^{206}}Pb-ratio.}
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
#' \item{w}{the overdispersion, i.e., a three-element vector with the
#' estimated standard deviation of the (assumedly) Normal distribution
#' that underlies the true isochron; and the lower and upper
#' half-widths of its \eqn{100(1-\alpha)\%} confidence interval (only
#' relevant if \code{model = 3}).}
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
#' uncertainty (e.g., decay constants)?
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
#' \itemize{
#'
#' \item The first element is a boolean flag indicating whether the
#' discordia line should be anchored. If this is \code{FALSE}, then
#' the second item is ignored and both the common Pb composition and
#' age are estimated.
#'
#' \item If the first element is \code{TRUE} and the second element is
#' \code{NA}, then the common Pb composition is fixed at the values
#' stored in \code{settings('iratio',...)}.
#'
#' item If the first element is \code{TRUE} and the second element is a
#' number, then the discordia line is forced to intersect the
#' concordia line at an age equal to that number.
#' }
#'
#' @seealso \code{\link{concordia}}, \code{\link{titterington}},
#'     \code{\link{isochron}}
#' @rdname ludwig
#' @export
ludwig.UPb <- function(x,exterr=FALSE,alpha=0.05,model=1,anchor=list(FALSE,NA),...){
    fit <- get.ta0b0(x,exterr=exterr,model=model,anchor=anchor)
    out <- fit[c('par','w','model')]
    out$cov <- fisher.lud(x,fit=fit,anchor=anchor)
    out$n <- length(x)
    mswd <- mswd.lud(fit$par,x=x,anchor=anchor)
    out <- c(out,mswd)
    if (model==3){
        out$w <- c(fit$w,profile_LL_discordia_disp(fit,x=x,alpha=alpha))
        names(out$w) <- c('s','ll','ul')
    }
    if (x$format %in% c(1,2,3)) parnames <- c('t','76i')
    else if (x$format %in% c(4,5,6)) parnames <- c('t','64i','74i')
    else if (x$format %in% c(7,8)) parnames <- c('t','68i','78i')
    else stop("Illegal input format.")
    names(out$par) <- parnames
    rownames(out$cov) <- parnames
    colnames(out$cov) <- parnames
    out
}

mswd.lud <- function(ta0b0,x,anchor=list(FALSE,NA)){
    ns <- length(x)
    # Mistake in Ludwig (1998)? Multiply the following by 2?
    SS <- LL.lud.UPb(ta0b0,x=x,exterr=FALSE,w=0,LL=FALSE)
    out <- list()
    anchored <- anchor[[1]]
    tanchored <- is.numeric(anchor[[2]])
    if (x$format<4){
        if (anchored) out$df <- ns-1
        else out$df <- ns-2
    } else if (x$format>3){
        if (anchored && tanchored) out$df <- 2*ns-1
        else if (tanchored) out$df <- 2*ns-2
        else out$df <- 2*ns-3
    } else {
        stop('Incorrect input format')
    }
    if (out$df>0){
        out$mswd <- as.vector(SS/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    } else {
        out$mswd <- 1
        out$p.value <- 1
    }
    out
}

get.ta0b0 <- function(x,exterr=FALSE,model=1,anchor=list(FALSE,NA)){
    init <- get.ta0b0.model2(x,anchor=anchor)
    if (model==1)
        out <- get.ta0b0.model1(x,init=init,exterr=exterr,anchor=anchor)
    else if (model==2)
        out <- list(par=init,w=0)
    else if (model==3)
        out <- get.ta0b0.model3(x,init=init,exterr=exterr,anchor=anchor)
    out$model <- model
    out$exterr <- exterr
    out
}

get.ta0b0.model1 <- function(x,init,exterr=FALSE,anchor=list(FALSE,NA)){
    out <- fit_ludwig_discordia(x,init=init,w=0,exterr=exterr,anchor=anchor)
    out$w <- 0
    out
}
get.ta0b0.model2 <- function(x,anchor=list(FALSE,NA)){
    if (x$format %in% c(1,2,3))
        out <- get.ta0b0.model2.2D(x,anchor=anchor)
    else if (x$format %in% c(4,5,6,7,8))
        out <- get.ta0b0.model2.3D(x,anchor=anchor)
    else
        stop("Illegal input format.")
    out
}
get.ta0b0.model3 <- function(x,init,exterr=FALSE,anchor=list(FALSE,NA)){
    out <- list(w=0,par=init)
#    for (i in 1:5){ # loop for more accurate but slower and more unstable results
    out <- fit_ludwig_discordia(x,init=out$par,w=out$w,exterr=exterr,anchor=anchor)
    out$w <- stats::optimise(LL.lud.disp,interval=c(0,1),x=x,ta0b0=out$par,
                             exterr=exterr,anchor=anchor,maximum=TRUE)$maximum
#    }
    out
}

get.ta0b0.model2.2D <- function(x,anchor=list(FALSE,NA)){
    xy <- data2york(x,option=2)[,c('X','Y'),drop=FALSE]
    if (!anchor[[1]]) {
        xyfit <- stats::lm(xy[,'Y'] ~ xy[,'X'])
        intercept <- xyfit$coef[1]
        slope <- xyfit$coef[2]
        ta0b0 <- concordia.intersection.ab(intercept,slope,wetherill=FALSE,d=x$d)
    } else if (is.na(anchor[[2]])){
        intercept <- settings('iratio','Pb207Pb206')[1]
        xyfit <- stats::lm(I(xy[,'Y']-intercept) ~ 0 + xy[,'X'])
        slope <- xyfit$coef
        ta0b0 <- concordia.intersection.ab(intercept,slope,wetherill=FALSE,d=x$d)
    } else if (is.numeric(anchor[[2]])){
        TW <- age_to_terawasserburg_ratios(anchor[[2]],st=0,exterr=FALSE,d=x$d)
        xyfit <- stats::lm(I(xy[,'Y']-TW$x['Pb207Pb206']) ~
                           0 + I(xy[,'X']-TW$x['U238Pb206']))
        slope <- xyfit$coef
        intercept <- TW$x['Pb207Pb206'] - slope*TW$x['U238Pb206']
        ta0b0 <- c(anchor[[2]],intercept)
    }
    ta0b0
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
            a0 <- settings('iratio','Pb206Pb208')[1]
            b0 <- settings('iratio','Pb207Pb208')[1]
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
SS.model2.3D <- function(tt,x,a0=NA,b0=NA){
    fits <- model2fit.3D(tt=tt,x=x,a0=a0,b0=b0)
    fits$SS
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

fit_ludwig_discordia <- function(x,init,w=0,exterr=FALSE,anchor=list(FALSE,NA),...){
    optifix(parms=init,fn=LL.lud.UPb,gr=LL.lud.UPb.gr,method="L-BFGS-B",x=x,w=w,
            exterr=exterr,fixed=fixit(x,anchor),lower=c(0,0,0),upper=c(10000,100,100))
}

LL.lud.disp <- function(w,x,ta0b0,exterr=FALSE,anchor=list(FALSE,NA)){
    # these two lines produce slightly more accurate but much slower results:
    # fit <- fit_ludwig_discordia(x,init=ta0b0,w=w,exterr=exterr,anchor=anchor)
    # LL.lud.UPb(ta0b0=fit$par,x=x,exterr=exterr,w=w,LL=TRUE)
    LL.lud.UPb(ta0b0=ta0b0,x=x,exterr=exterr,w=w,LL=TRUE)
}
LL.lud.UPb <- function(ta0b0,x,exterr=FALSE,w=0,LL=FALSE){
    if (x$format %in% c(1,2,3)){
        return(LL.lud.2D(ta0b0,x=x,exterr=exterr,w=w,LL=LL))
    } else if (x$format %in% c(4,5,6)){
        return(LL.lud.3D(ta0b0,x=x,exterr=exterr,w=w,LL=LL))
    } else if (x$format %in% c(7,8)){
        return(LL.lud.Th(ta0b0,x=x,exterr=exterr,w=w,LL=LL))
    } else {
        stop('Incorrect input format.')
    }
}
LL.lud.2D <- function(ta0,x,exterr=FALSE,w=0,LL=FALSE){
    l <- data2ludwig(x,tt=ta0[1],a0=ta0[2],exterr=exterr,w=w)
    R <- rbind(l$rx,l$ry)
    out <- t(R) %*% l$omega %*% R
    if (LL){
        k <- length(R)
        detE <- determinant(2*pi*l$omegainv,logarithm=TRUE)$modulus
        out <- -0.5*(out + k*log(2*pi) + detE)
    }
    out
}
LL.lud.3D <- function(ta0b0,x,exterr=FALSE,w=0,LL=FALSE){
    l <- data2ludwig(x,tt=ta0b0[1],a0=ta0b0[2],b0=ta0b0[3],exterr=exterr,w=w)
    phi <- l$phi
    R <- l$R
    r <- l$r
    out <- 0
    ns <- length(l$R)
    omega <- l$omega
    for (i in 1:ns){
        i1 <- i
        i2 <- i + ns
        i3 <- i + 2*ns
        for (j in 1:ns){
            j1 <- j
            j2 <- j + ns
            j3 <- j + 2*ns
            out <- out + R[i]*R[j]*omega[i1,j1] +
                r[i]*r[j]*omega[i2,j2] + phi[i]*phi[j]*omega[i3,j3] +
                2*( R[i]*r[j]*omega[i1,j2] +
                    R[i]*phi[j]*omega[i1,j3] + r[i]*phi[j]*omega[i2,j3] )
        }
    }
    if (LL){
        k <- 3*ns
        detE <- determinant(l$omegainv,logarithm=TRUE)$modulus
        out <- -0.5*(out + k*log(2*pi) + detE)
    }
    out
}
LL.lud.Th <- function(ta0b0,x,exterr=FALSE,w=0,LL=FALSE){
    l <- data2ludwig(x,tt=ta0b0[1],a0=ta0b0[2],b0=ta0b0[3],exterr=exterr,w=w)
    KLM <- matrix(c(l$K,l$L,l$M),nrow=1)
    out <- KLM %*% l$omega %*% t(KLM)
    if (LL){
        k <- length(KLM)
        detE <- determinant(l$omegainv,logarithm=TRUE)$modulus
        out <- -0.5*(out + k*log(2*pi) + detE)
    }
    out
}

LL.lud.UPb.gr <- function(ta0b0,x,exterr=FALSE,w=0){
    if (x$format %in% c(1,2,3)){
        return(LL.lud.2D.gr(ta0b0,x,exterr=exterr,w=w))
    } else if (x$format %in% c(4,5,6)){
        return(LL.lud.3D.gr(ta0b0,x=x,exterr=exterr,w=w))
    } else if (x$format %in% c(7,8)){
        return(LL.lud.Th.gr(ta0b0,x=x,exterr=exterr,w=w))
    } else {
        stop('Incorrect input format.')
    }
}
LL.lud.2D.gr <- function(ta0,x,exterr=FALSE,w=0){
    tt <- ta0[1]
    a0 <- ta0[2]
    l5 <- settings('lambda','U235')[1]
    l8 <- settings('lambda','U238')[1]
    U <- settings('iratio','U238U235')[1]
    l <- data2ludwig(x,tt=tt,a0=a0,exterr=exterr,w=w)
    rx <- l$rx
    X <- get.U238Pb206.ratios(x)[,1]
    ry <- l$ry # ry = Y - a - bx
    ns <- length(x)
    dRdt <- matrix(0,2*ns,1)
    dRda0 <- matrix(0,2*ns,1)
    D <- mclean(tt,d=x$d)
    # b = D$Pb207U235/U - a0*D$Pb206U238
    dbdt <- D$dPb207U235dt/U - a0*D$dPb206U238dt
    drydt <- -dbdt*(X-rx) # x = X - rX
    dryda0 <- D$Pb206U238 - 1
    dRdt[(ns+1):(2*ns),1] <- drydt
    dRda0[(ns+1):(2*ns),1] <- dryda0
    dSdt <- 2*t(dRdt) %*% l$omega %*% dRdt
    dSda0 <- 2*t(dRda0) %*% l$omega %*% dRda0
    c(dSdt,dSda0)
}
LL.lud.3D.gr <- function(ta0b0,x,exterr=FALSE,w=0){
    tt <- ta0b0[1]
    a0 <- ta0b0[2]
    b0 <- ta0b0[3]
    l2 <- settings('lambda','Th232')[1]
    l5 <- settings('lambda','U235')[1]
    l8 <- settings('lambda','U238')[1]
    U <- settings('iratio','U238U235')[1]
    l <- data2ludwig(x,tt=tt,a0=0,b0=b0,exterr=exterr,w=w)
    phi <- l$phi
    R <- l$R
    r <- l$r
    omega <- l$omega
    Z <- get.Pb204U238.ratios(x)[,1]
    z <- Z - phi
    dRdt <- -l5*exp(l5*tt)
    drdt <- -l8*exp(l8*tt)
    drda0 <- -z
    dRdb0 <- -U*z
    dSdt <- 0
    dSda0 <- 0
    dSdb0 <- 0
    ns <- length(x)
    for (i in 1:ns){
        i1 <- i
        i2 <- i + ns
        i3 <- i + 2*ns
        for (j in 1:ns){
            j1 <- j
            j2 <- j + ns
            j3 <- j + 2*ns
            dSdt <- dSdt + dRdt*R[j]*omega[i1,j1] + R[i]*dRdt*omega[i1,j1] +
                drdt*r[j]*omega[i2,j2] + r[i]*drdt*omega[i2,j2] +
                2*( dRdt*r[j]*omega[i1,j2] + R[i]*drdt*omega[i1,j2] +
                    dRdt*phi[j]*omega[i1,j3] + drdt*phi[j]*omega[i2,j3] )
            dSda0 <- dSda0 + drda0[i]*r[j]*omega[i2,j2] + r[i]*drda0[j]*omega[i2,j2] +
                2*( R[i]*drda0[j]*omega[i1,j2] + drda0[i]*phi[j]*omega[i2,j3] )
            dSdb0 <- dSdb0 + dRdb0[i]*R[j]*omega[i1,j1] + R[i]*dRdb0[j]*omega[i1,j1] +
                2*( dRdb0[i]*r[j]*omega[i1,j2] + dRdb0[i]*phi[j]*omega[i1,j3] )
        }
    }
    c(dSdt,dSda0,dSdb0)
}
LL.lud.Th.gr <- function(ta0b0,x,exterr=FALSE,w=0){
    tt <- ta0b0[1]
    a0 <- ta0b0[2]
    b0 <- ta0b0[3]
    l <- data2ludwig_Th(x,tt=tt,a0=a0,b0=b0,w=w,exterr=exterr)
    ns <- length(l$K)
    l5 <- settings('lambda','U235')[1]
    l8 <- settings('lambda','U238')[1]
    l2 <- settings('lambda','Th232')[1]
    U <- settings('iratio','U238U235')[1]
    K <- matrix(l$K,nrow=ns)
    L <- matrix(l$L,nrow=ns)
    M <- matrix(l$M,nrow=ns)
    Z <- get.Pb208Th232.ratios(x)[,1,drop=FALSE]
    c0 <- Z - exp(l2*tt) + 1 - M
    W <- x$x[,'Th232U238']
    dKdt <- matrix(rep(-l5*exp(l5*tt),ns),nrow=ns)
    dLdt <- matrix(rep(-l8*exp(l8*tt),ns),nrow=ns)
    dMdt <- matrix(rep(-l2*exp(l2*tt),ns),nrow=ns)
    dKdb0 <- -U*c0*W
    dLda0 <- -c0*W
    dSdt <- 2*t(K)%*%l$omega11%*%dKdt +
        t(K)%*%l$omega12%*%dLdt + t(dLdt)%*%t(l$omega12)%*%dKdt +
        t(K)%*%l$omega13%*%dMdt + t(dMdt)%*%t(l$omega13)%*%dKdt +
        t(L)%*%l$omega21%*%dKdt + t(dKdt)%*%t(l$omega21)%*%dLdt +
        2*t(L)%*%l$omega22%*%dLdt +
        t(L)%*%l$omega23%*%dMdt + t(dMdt)%*%t(l$omega23)%*%dLdt +
        t(M)%*%l$omega31%*%dKdt + t(dKdt)%*%t(l$omega31)%*%dMdt +
        t(M)%*%l$omega32%*%dLdt + t(dLdt)%*%t(l$omega32)%*%dMdt +
        2*t(M)%*%l$omega33%*%dMdt
    dSda0 <- 2*t(L)%*%l$omega22%*%dLda0 +
        t(K)%*%l$omega12%*%dLda0 + t(K)%*%t(l$omega21)%*%dLda0 +
        t(M)%*%l$omega32%*%dLda0 + t(M)%*%t(l$omega23)%*%dLda0        
    dSdb0 <- 2*t(K)%*%l$omega11%*%dKdb0 +
        t(L)%*%l$omega21%*%dKdb0 + t(L)%*%t(l$omega12)%*%dKdb0 + 
        t(M)%*%l$omega31%*%dKdb0 + t(M)%*%t(l$omega13)%*%dKdb0
    c(dSdt,dSda0,dSdb0)/2
}

fisher.lud <- function(x,...){ UseMethod("fisher.lud",x) }
fisher.lud.default <- function(x,...){
    stop( "No default method available (yet)." )
}
fisher.lud.UPb <- function(x,fit,exterr=TRUE,anchor=list(FALSE,NA),...){
    ns <- length(x)
    if (x$format %in% c(1,2,3)){
        fish <- fisher_lud_2D(x,fit)
        AA <- fish[1:ns,1:ns]
        BB <- fish[1:ns,(ns+1):(ns+2)]
        CC <- fish[(ns+1):(ns+2),1:ns]
        DD <- fish[(ns+1):(ns+2),(ns+1):(ns+2)]
    } else if (x$format %in% c(4,5,6)){
        fish <- fisher_lud_3D(x,fit)
        AA <- fish[1:ns,1:ns]
        BB <- fish[1:ns,(ns+1):(ns+3)]
        CC <- fish[(ns+1):(ns+3),1:ns]
        DD <- fish[(ns+1):(ns+3),(ns+1):(ns+3)]
    } else if (x$format %in% c(7,8)){
        fish <- fisher_lud_Th(x,fit)
        AA <- fish[1:ns,1:ns]
        BB <- fish[1:ns,(ns+1):(ns+3)]
        CC <- fish[(ns+1):(ns+3),1:ns]
        DD <- fish[(ns+1):(ns+3),(ns+1):(ns+3)]
    } else {
        stop("Illegal data format.")
    }
    anchorfish(AA,BB,CC,DD,anchor=anchor)
}
anchorfish <- function(AA,BB,CC,DD,anchor){
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
fisher_lud_2D <- function(x,fit){
    tt <- fit$par[1]
    a0 <- fit$par[2]
    l <- data2ludwig_2D(x,tt=tt,a0=a0,w=0,exterr=fit$exterr)
    l5 <- settings('lambda','U235')[1]
    l8 <- settings('lambda','U238')[1]
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    ones <- matrix(1,ns,1)
    D <- mclean(tt,d=x$d)
    BB <- D$Pb207U235/U - a0*D$Pb206U238
    CC <- -t(l$X)%*%l$omega12 + t(l$x)%*%l$omega12 + t(l$x)%*%l$omega21 -
        t(l$Y)%*%l$omega22 + t(ones)%*%l$omega22*a0 + 2*t(l$x)%*%l$omega22*BB
    DD <- -t(l$X)%*%l$omega12%*%l$x + t(l$x)%*%l$omega12%*%l$x -
        t(l$Y)%*%l$omega22%*%l$x + t(ones)%*%l$omega22%*%l$x*a0 -
        t(l$x)%*%l$omega21%*%l$X + t(l$x)%*%l$omega21%*%l$x -
        t(l$x)%*%l$omega22%*%l$Y + t(l$x)%*%l$omega22%*%ones*a0 +
        2*t(l$x)%*%l$omega22%*%l$x*BB
    dBdt <- D$dPb207U235dt/U - a0*D$dPb206U238dt
    dBda0 <- -D$Pb206U238
    d2Bda0dt <- -D$dPb206U238dt
    d2Bdtda0 <- d2Bda0dt
    d2Bdt2 <- D$d2Pb207U235dt2/U - a0*D$d2Pb206U238dt2
    d2Bda02 <- 0
    dDda0 <- t(ones)%*%l$omega22%*%l$x +
        t(l$x)%*%l$omega22%*%ones + 2*t(l$x)%*%l$omega22%*%l$x*dBda0
    ns <- length(x)
    out <- matrix(0,ns+2,ns+2)
    out[1:ns,1:ns] <- # d2S/dx2
        2*(l$omega11 + l$omega21*BB + l$omega12*BB + l$omega22*BB^2)
    out[ns+2,1:ns] <-   # d2S/dxda0
        2*(CC*dBda0 + t(ones)%*%l$omega21 + t(ones)%*%l$omega22*BB)
    out[ns+2,ns+2] <- # d2S/da02
        2*(t(ones)%*%l$omega22%*%ones + t(ones)%*%l$omega22%*%l$x*dBda0 +
           t(l$x)%*%l$omega22%*%ones*dBda0 + t(l$x)%*%l$omega22%*%l$x*(dBda0^2))
    out[ns+1,1:ns] <- 2*CC*dBdt # d2S/dxdt
    out[ns+1,ns+1] <- # d2S/dt2
        DD*d2Bdt2 + 2*t(l$x)%*%l$omega22%*%l$x*(dBdt^2)
    out[ns+1,ns+2] <- # d2S/dt0da0
        DD*d2Bdtda0 + dDda0*dBdt
    out[1:ns,ns+1] <- t(out[ns+1,1:ns]) # d2S/dtdx
    out[1:ns,ns+2] <- t(out[ns+2,1:ns]) # d2S/da0dx
    out[ns+2,ns+1] <- out[ns+1,ns+2]
    out/2
}
fisher_lud_3D <- function(x,fit){
    tt <- fit$par[1]
    a0 <- fit$par[2]
    b0 <- fit$par[3]
    l <- data2ludwig_3D(x,tt=tt,a0=a0,b0=b0,w=0,exterr=fit$exterr)
    z <- l$z
    omega <- l$omega
    ns <- length(z)
    l5 <- settings('lambda','U235')
    l8 <- settings('lambda','U238')
    U <- settings('iratio','U238U235')[1]
    Q235 <- l5[1]*exp(l5[1]*tt)
    Q238 <- l8[1]*exp(l8[1]*tt)
    d2L.dt2 <- 0
    d2L.da02 <- 0
    d2L.db02 <- 0
    d2L.da0dt <- 0
    d2L.db0dt <- 0
    d2L.da0db0 <- 0
    out <- matrix(0,ns+3,ns+3)
    for (i in 1:ns){
        i1 <- i
        i2 <- i + ns
        i3 <- i + 2*ns
        d2L.dzdt <- 0
        d2L.dzda0 <- 0
        d2L.dzdb0 <- 0
        d2L.dzdz <- 0
        for (j in 1:ns){
            j1 <- j
            j2 <- j + ns
            j3 <- j + 2*ns
            d2L.dt2 <- d2L.dt2 + omega[i1,j1]*Q235^2 +
                omega[i2,j2]*Q238^2 + 2*Q235*Q238*omega[i1,j2]
            d2L.da02 <- d2L.da02 + z[i]*z[j]*omega[i2,j2]
            d2L.db02 <- d2L.db02 + z[i]*z[j]*omega[i1,j1]*U^2
            d2L.da0dt <- d2L.da0dt +
                Q238*(z[i]+z[j])*omega[i2,j2]/2 +
                Q235*z[j]*omega[i1,j2]
            d2L.db0dt <- d2L.db0dt + U*(
                Q235*(z[i]+z[j])*omega[i1,j1]/2 +
                Q238*z[i]*omega[i1,j2] )
            d2L.da0db0 <- d2L.da0db0 +
                U * z[i]*z[j]*omega[i1,j2]
            d2L.dzdt <- d2L.dzdt +
                Q235*(U*b0*omega[i1,j1] + a0*omega[i2,j1] + omega[i3,j1]) +
                Q238*(U*b0*omega[i1,j2] + a0*omega[i2,j2] + omega[i3,j2])
            d2L.dzda0 <- d2L.dzda0 +
                z[j]*(a0*omega[i2,j2] + U*b0*omega[i1,j2] + omega[i3,j2])
            d2L.dzdb0 <- d2L.dzdb0 +
                U*z[j]*(a0*omega[i2,j1] + U*b0*omega[i1,j1] + omega[i3,j1])
            d2L.dzdz <- (omega[i1,j1] + omega[i1,j3] +
                         omega[i3,j1])*(U*b0)^2 +
                         omega[i2,j2]*a0^2 + omega[i3,j3] +
                         a0*(omega[i2,j3]+omega[i3,j2]) +
                         a0*U*b0*(omega[i1,j2]+omega[i2,j1])
            out[i,j] <- d2L.dzdz # dij
            out[j,i] <- d2L.dzdz # dij
        }
        d2L.dz2 <- omega[i1,i1]*(U*b0)^2 + omega[i2,i2]*a0^2 +
                   omega[i3,i3] + 2*(a0*U*b0*omega[i1,i2] +
                   U*b0*omega[i1,i3] + a0*omega[i3,i3])
        out[i,i] <- d2L.dz2 # dii
        out[i,ns+1] <- d2L.dzdt # ni1
        out[i,ns+2] <- d2L.dzda0 # ni2
        out[i,ns+3] <- d2L.dzdb0 # ni3
        out[ns+1,i] <- d2L.dzdt # ni1
        out[ns+2,i] <- d2L.dzda0 # ni2
        out[ns+3,i] <- d2L.dzdb0 # ni3
    }
    out[ns+1,ns+1] <- d2L.dt2 # m11
    out[ns+2,ns+2] <- d2L.da02 # m22
    out[ns+3,ns+3] <- d2L.db02 # m33
    out[ns+1,ns+2] <- d2L.da0dt # m12
    out[ns+1,ns+3] <- d2L.db0dt # m13
    out[ns+2,ns+1] <- d2L.da0dt # m12
    out[ns+3,ns+1] <- d2L.db0dt # m13
    out[ns+2,ns+3] <- d2L.da0db0 # m23
    out[ns+3,ns+2] <- d2L.da0db0 # m23
    out
}
fisher_lud_Th <- function(x,fit){
    tt <- fit$par[1]
    a0 <- fit$par[2]
    b0 <- fit$par[3]
    l <- data2ludwig_Th(x,tt=tt,a0=a0,b0=b0,w=0,exterr=fit$exterr)
    ns <- length(l$K)
    l5 <- settings('lambda','U235')[1]
    l8 <- settings('lambda','U238')[1]
    l2 <- settings('lambda','Th232')[1]
    U <- settings('iratio','U238U235')[1]
    K <- matrix(l$K,nrow=ns)
    L <- matrix(l$L,nrow=ns)
    M <- matrix(l$M,nrow=ns)
    Z <- get.Pb208Th232.ratios(x)[,1,drop=FALSE]
    c0 <- Z - exp(l2*tt) + 1 - M
    W <- x$x[,'Th232U238']
    dKdt <- matrix(rep(-l5*exp(l5*tt),ns),nrow=ns) # TODO add disequilibrium correction, do the same for the other two Fisher functions!
    dLdt <- matrix(rep(-l8*exp(l8*tt),ns),nrow=ns)
    dMdt <- matrix(rep(-l2*exp(l2*tt),ns),nrow=ns)
    d2Kdt2 <- l5*dKdt
    d2Ldt2 <- l8*dLdt
    d2Mdt2 <- l2*dMdt
    dKdb0 <- -U*c0*W
    dKdc0 <- diag(-U*b0*W)
    d2Kdc0db0 <- diag(-U*W)
    dLda0 <- -c0*W
    dLdc0 <- diag(-a0*W)
    d2Ldc0da0 <- diag(-W)
    dMdc0 <- diag(rep(-1,ns))
    d2Sdc0dt <- 2*t(dKdc0)%*%l$omega11%*%dKdt +
        t(dLdc0)%*%t(l$omega12)%*%dKdt + t(dKdc0)%*%l$omega12%*%dLdt +
        t(dMdc0)%*%t(l$omega13)%*%dKdt + t(dKdc0)%*%l$omega13%*%dMdt +
        t(dKdc0)%*%t(l$omega21)%*%dLdt + t(dLdc0)%*%l$omega21%*%dKdt +
        2*t(dLdc0)%*%t(l$omega22)%*%dLdt +
        t(dMdc0)%*%t(l$omega23)%*%dLdt + t(dLdc0)%*%l$omega23%*%dMdt +
        t(dKdc0)%*%t(l$omega31)%*%dMdt + t(dMdc0)%*%l$omega31%*%dKdt +
        t(dLdc0)%*%t(l$omega32)%*%dMdt + t(dMdc0)%*%l$omega32%*%dLdt +
        2*t(dMdc0)%*%t(l$omega33)%*%dMdt
    d2Sdc0da0 <- t(dKdc0)%*%l$omega12%*%dLda0 +
        t(dKdc0)%*%t(l$omega21)%*%dLda0 + 
        2*t(dLdc0)%*%t(l$omega22)%*%dLda0 +
        t(dMdc0)%*%t(l$omega23)%*%dLda0 +
        t(dMdc0)%*%l$omega32%*%dLda0
    d2Sdc0db0 <- 2*t(dKdc0)%*%l$omega11%*%dKdb0 +
        t(dLdc0)%*%t(l$omega12)%*%dKdb0 + 
        t(dMdc0)%*%t(l$omega13)%*%dKdb0 +
        t(dLdc0)%*%l$omega21%*%dKdb0 +
        t(dMdc0)%*%l$omega31%*%dKdb0
    d2Sdc02 <- 2*t(dKdc0)%*%l$omega11%*%dKdc0 +
        t(dLdc0)%*%t(l$omega12)%*%dKdc0 + t(dKdc0)%*%l$omega12%*%dLdc0 +
        t(dMdc0)%*%t(l$omega13)%*%dKdc0 + t(dKdc0)%*%l$omega13%*%dMdc0 +
        t(dKdc0)%*%t(l$omega21)%*%dLdc0 + t(dLdc0)%*%l$omega21%*%dKdc0 +
        2*t(dLdc0)%*%l$omega22%*%dLdc0 +
        t(dMdc0)%*%t(l$omega23)%*%dLdc0 + t(dLdc0)%*%l$omega23%*%dMdc0 +
        t(dKdc0)%*%t(l$omega31)%*%dMdc0 + t(dMdc0)%*%l$omega31%*%dKdc0 +
        t(dLdc0)%*%t(l$omega32)%*%dMdc0 + t(dMdc0)%*%l$omega32%*%dLdc0 +
        2*t(dMdc0)%*%l$omega33%*%dMdc0
    d2Sdt2 <- 2*t(K)%*%l$omega11%*%d2Kdt2 + 2*t(dKdt)%*%l$omega11%*%dKdt +
        t(K)%*%l$omega12%*%d2Ldt2 + t(dLdt)%*%t(l$omega12)%*%dKdt +
        t(dLdt)%*%t(l$omega12)%*%d2Kdt2 + t(dKdt)%*%l$omega12%*%dLdt +
        t(K)%*%l$omega13%*%d2Mdt2 + t(dMdt)%*%t(l$omega13)%*%dKdt +
        t(dMdt)%*%t(l$omega13)%*%d2Kdt2 + t(dKdt)%*%l$omega13%*%dMdt +
        t(L)%*%l$omega21%*%d2Kdt2 + t(dKdt)%*%t(l$omega21)%*%dLdt +
        t(dKdt)%*%t(l$omega21)%*%d2Ldt2 + t(dLdt)%*%l$omega21%*%dKdt +
        2*t(L)%*%l$omega22%*%d2Ldt2 + 2*t(dLdt)%*%l$omega22%*%dLdt +
        t(L)%*%l$omega23%*%d2Mdt2 + t(dMdt)%*%t(l$omega23)%*%dLdt +
        t(dMdt)%*%t(l$omega23)%*%d2Ldt2 + t(dLdt)%*%l$omega23%*%dMdt +
        t(M)%*%l$omega31%*%d2Kdt2 + t(dKdt)%*%t(l$omega31)%*%dMdt +
        t(dKdt)%*%t(l$omega31)%*%d2Mdt2 + t(dMdt)%*%l$omega31%*%dKdt +
        t(M)%*%l$omega32%*%d2Ldt2 + t(dLdt)%*%t(l$omega32)%*%dMdt +
        t(dLdt)%*%t(l$omega32)%*%d2Mdt2 + t(dMdt)%*%l$omega32%*%dLdt +
        2*t(M)%*%l$omega33%*%d2Mdt2 + 2*t(dMdt)%*%l$omega33%*%dMdt
    d2Sdtda0 <- t(dKdt)%*%t(l$omega21)%*%dLda0 + 
        2*t(dLdt)%*%l$omega22%*%dLda0 +
        t(dMdt)%*%t(l$omega23)%*%dLda0 
    d2Sdtdb0 <- 2*t(dKdt)%*%l$omega11%*%dKdb0 +
        t(dLdt)%*%t(l$omega12)%*%dKdb0 + 
        t(dMdt)%*%t(l$omega13)%*%dKdb0
    d2Sda02 <- 2*t(dLda0)%*%l$omega22%*%dLda0
    d2Sdb02 <- 2*t(dKdb0)%*%l$omega11%*%dKdb0
    d2Sda0db0 <- t(dLda0)%*%l$omega12%*%dKdb0 + t(dLda0)%*%l$omega21%*%dKdb0
    out <- matrix(0,ns+3,ns+3)
    out[1:ns,1:ns] <- d2Sdc02
    out[1:ns,ns+1] <- d2Sdc0dt
    out[1:ns,ns+2] <- d2Sdc0da0
    out[1:ns,ns+3] <- d2Sdc0db0
    out[ns+1,ns+1] <- d2Sdt2
    out[ns+2,ns+2] <- d2Sda02
    out[ns+3,ns+3] <- d2Sdb02
    out[ns+1,ns+2] <- d2Sdtda0
    out[ns+1,ns+3] <- d2Sdtdb0
    out[ns+2,ns+3] <- d2Sda0db0
    out[ns+1,1:ns] <- out[1:ns,ns+1]
    out[ns+2,1:ns] <- out[1:ns,ns+2]
    out[ns+3,1:ns] <- out[1:ns,ns+3]
    out[ns+1,ns+2] <- out[ns+2,ns+1]
    out[ns+1,ns+3] <- out[ns+3,ns+1]
    out[ns+2,ns+3] <- out[ns+3,ns+2]
    out/2
}

data2ludwig <- function(x,...){ UseMethod("data2ludwig",x) }
data2ludwig.default <- function(x,...){ stop('default function undefined') }
data2ludwig.UPb <- function(x,tt,a0,b0=0,g0=rep(0,length(x)),exterr=FALSE,w=0,...){
    if (x$format %in% c(1,2,3))
        out <- data2ludwig_2D(x,tt=tt,a0=a0,w=w,exterr=exterr)
    else if (x$format %in% c(4,5,6))
        out <- data2ludwig_3D(x,tt=tt,a0=a0,b0=b0,w=w,exterr=exterr)
    else if (x$format %in% c(7,8))
        out <- data2ludwig_Th(x,tt=tt,a0=a0,b0=b0,w=w,exterr=exterr)
    else stop('Incorrect input format.')
    out
}
data2ludwig_2D <- function(x,tt,a0,w=0,exterr=FALSE){
    l5 <- settings('lambda','U235')
    l8 <- settings('lambda','U238')
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    E <- matrix(0,2*ns+2,2*ns+2)
    J <- matrix(0,2*ns,2*ns+2)
    X <- matrix(0,ns,1)
    Y <- matrix(0,ns,1)
    for (i in 1:ns){
        D <- mclean(tt=tt,d=x$d[i],exterr=exterr)
        B <-  D$Pb207U235/U - a0*D$Pb206U238
        XY <- tera.wasserburg(x,i)
        X[i] <- XY$x['U238Pb206']
        Y[i] <- XY$x['Pb207Pb206']
        E[(2*i-1):(2*i),(2*i-1):(2*i)] <- XY$cov[1:2,1:2]
        E[2*i,2*i] <- E[2*i,2*i] + (a0*w)^2
        J[i,2*i-1] <- 1  # drx/dX
        J[ns+i,2*i] <- 1 # dry/dY
        dBdl5 <- D$dPb207U235dl35/U
        dBdl8 <- -a0*D$dPb206U238dl38
        J[ns+i,2*ns+1] <- -dBdl5*X[i]
        J[ns+i,2*ns+2] <- -dBdl8*X[i]
    }
    E[2*ns+1,2*ns+1] <- l5[2]^2
    E[2*ns+2,2*ns+2] <- l8[2]^2
    out <- list()
    OI <- J %*% E %*% t(J)
    O <- blockinverse(AA=OI[1:ns,1:ns],
                      BB=OI[1:ns,(ns+1):(2*ns)],
                      CC=OI[(ns+1):(2*ns),1:ns],
                      DD=OI[(ns+1):(2*ns),(ns+1):(2*ns)],
                      doall=TRUE)
    out$omegainv <- OI
    out$omega <- O
    out$omega11 <- O[1:ns,1:ns]
    out$omega12 <- O[1:ns,(ns+1):(2*ns)]
    out$omega21 <- O[(ns+1):(2*ns),1:ns]
    out$omega22 <- O[(ns+1):(2*ns),(ns+1):(2*ns)]
    left <- out$omega11 + B*(out$omega12+out$omega21) + out$omega22*B^2
    right <- out$omega11%*%X + B*out$omega12%*%X +
        out$omega21%*%(Y-a0) + B*out$omega22%*%(Y-a0)
    out$x <- solve(left,right)
    out$X <- X
    out$Y <- Y
    out$rx <- X - out$x
    out$ry <- Y - a0 - B*out$x
    out
}
data2ludwig_3D <- function(x,tt,a0,b0,w=0,exterr=FALSE){
    l5 <- settings('lambda','U235')
    l8 <- settings('lambda','U238')
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    E <- matrix(0,3*ns+2,3*ns+2)
    J <- matrix(0,3*ns,3*ns+2)
    R <- rep(0,ns)
    r <- rep(0,ns)
    Z <- rep(0,ns)
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        D <- mclean(tt=tt,d=x$d[i],exterr=exterr)
        Z[i] <- wd$x['Pb204U238']
        R[i] <- wd$x['Pb207U235'] - U*b0*Z[i] - D$Pb207U235
        r[i] <- wd$x['Pb206U238'] - a0*Z[i] - D$Pb206U238
        Ew <- get.Ew(w=w,Z=Z[i],a0=a0,b0=b0,U=U)
        E[(3*i-2):(3*i),(3*i-2):(3*i)] <- wd$cov + Ew
        J[i,3*i-2] <- 1                     # dRdX
        J[ns+i,3*i-1] <- 1                  # drdY
        J[2*ns+i,3*i] <- 1                  # dphidZ
        J[i,3*ns+1] <- -D$dPb207U235dl35    # dRdl5
        J[ns+i,3*ns+2] <- -D$dPb206U238dl38 # drdl8
    }
    E[3*ns+1,3*ns+1] <- l5[2]^2
    E[3*ns+2,3*ns+2] <- l8[2]^2
    omega <- solve(J %*% E %*% t(J))
    # rearrange sum of squares:
    V <- matrix(0,ns,ns)
    W <- rep(0,ns)
    for (i in 1:ns){
        i1 <- i
        i2 <- i + ns
        i3 <- i + 2*ns
        for (j in 1:ns){
            j1 <- j
            j2 <- j + ns
            j3 <- j + 2*ns
            W[i] <- W[i]-R[j]*(U*b0*omega[i1,j1]+a0*omega[i2,j1]+omega[i3,j1])-
                    r[j]*(U*b0*omega[i1,j2]+a0*omega[i2,j2]+omega[i3,j2])
            # Ludwig (1998): V[i,j]<-U*b0*omega[i1,j3]+a0*omega[i2,j3]+omega[i3,j3]
            V[i,j] <- omega[i1,j1]*(U*b0)^2+omega[i2,j2]*a0^2+omega[i3,j3] +
                      2*(U*a0*b0*omega[i1,j2]+U*b0*omega[i1,j3]+a0*omega[i2,j3])
        }
    }
    phi <- solve(V,W)
    z <- Z - phi
    out <- list(R=R,r=r,phi=phi,z=z,omega=omega,omegainv=E)    
}
data2ludwig_Th <- function(x,tt,a0,b0,w=0,exterr=FALSE){
    l2 <- settings('lambda','Th232')
    l5 <- settings('lambda','U235')
    l8 <- settings('lambda','U238')
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    K <- rep(0,ns)
    L <- rep(0,ns)
    M <- rep(0,ns)
    c0 <- rep(0,ns)
    W <- x$x[,'Th232U238']
    E <- matrix(0,4*ns+3,4*ns+3)
    J <- matrix(0,3*ns,4*ns+3)
    for (i in 1:ns){
        D <- mclean(tt,d=x$d[i],exterr=exterr)
        wd <- wetherill(x,i)
        ii <- 4*i-3
        c0[i] <- wd$x['Pb208Th232'] - exp(l2[1]*tt) + 1
        K[i] <- wd$x['Pb207U235'] - c0[i]*U*b0*W[i] - D$Pb207U235
        L[i] <- wd$x['Pb206U238'] - c0[i]*a0*W[i] - D$Pb206U238
        # M[i] <- wd$x['Pb208Th232'] - exp(l2[1]*tt) + 1 - c0[i]
        E[ii:(ii+3),ii:(ii+3)] <- wd$cov
        J[i,ii] <- 1 # dKi/dPb7U5
        J[i,ii+3] <- -c0[i]*U*b0 # dKi/dTh2U8
        J[ns+i,ii+1] <- 1 # dLi/dPb6U8
        J[ns+i,ii+3] <- -c0[i]*a0 # dLi/dTh2U8
        J[2*ns+i,ii+2] <- 1 # dMi/dPb8Th2
        J[i,4*ns+1] <- -D$dPb207U235dl35 # dKi/dl5
        J[ns+i,4*ns+2] <- -D$dPb206U238dl38 # dLi/dl8
        J[2*ns+i,4*ns+3] <- - tt*exp(l2[1]*tt) # dMi/dl2
    }
    E[4*ns+1,4*ns+1] <- l5[2]^2
    E[4*ns+2,4*ns+2] <- l8[2]^2
    E[4*ns+3,4*ns+3] <- l2[2]^2
    E2 <- J %*% E %*% t(J)
    if (w>0) E2 <- E2 + get.Ew_Th(w=w,W=x$x[,'Th232U238'],a0=a0,b0=b0,c0=c0)
    omega <- blockinverse3x3(AA=E2[1:ns,1:ns],
                             BB=E2[1:ns,(ns+1):(2*ns)],
                             CC=E2[1:ns,(2*ns+1):(3*ns)],
                             DD=E2[(ns+1):(2*ns),1:ns],
                             EE=E2[(ns+1):(2*ns),(ns+1):(2*ns)],
                             FF=E2[(ns+1):(2*ns),(2*ns+1):(3*ns)],
                             GG=E2[(2*ns+1):(3*ns),1:ns],
                             HH=E2[(2*ns+1):(3*ns),(ns+1):(2*ns)],
                             II=E2[(2*ns+1):(3*ns),(2*ns+1):(3*ns)])
    omega11 <- omega[1:ns,1:ns]
    omega12 <- omega[1:ns,(ns+1):(2*ns)]
    omega13 <- omega[1:ns,(2*ns+1):(3*ns)]
    omega21 <- omega[(ns+1):(2*ns),1:ns]
    omega22 <- omega[(ns+1):(2*ns),(ns+1):(2*ns)]
    omega23 <- omega[(ns+1):(2*ns),(2*ns+1):(3*ns)]
    omega31 <- omega[(2*ns+1):(3*ns),1:ns]
    omega32 <- omega[(2*ns+1):(3*ns),(ns+1):(2*ns)]
    omega33 <- omega[(2*ns+1):(3*ns),(2*ns+1):(3*ns)]
    AA <- W*right(omega11,W)*(U*b0)^2 + W*right(omega12,W)*U*a0*b0 +
        W*omega13*U*b0 + W*right(omega21,W)*U*a0*b0 + w*omega23*a0 +
        right(omega31,W)*U*b0 + right(omega32,W)*a0 + omega33
    BB <- U*b0*W*omega11%*%K + U*b0*W*omega12%*%L +
        a0*W*omega21%*%K + omega31%*%K + omega32%*%L
    CC <- U*b0*t(K)%*%right(omega11,W) + a0*t(K)%*%right(omega12,W) +
        t(K)%*%omega13 + U*b0*t(L)%*%right(omega21,W) + t(L)%*%omega21
    M <- -solve(t(AA), BB + t(CC))/2
    list(K=K,L=L,M=M,omega11=omega11,omega12=omega12,omega13=omega13,
         omega21=omega21,omega22=omega22,omega23=omega23,omega31=omega31,
         omega32=omega32,omega33=omega33,omega=omega,omegainv=E2)
}

get.Ew <- function(w,Z,a0,b0,U){
    E <- diag(c(a0,b0)*w)^2
    J <- matrix(0,3,2)
    J[1,2] <- -U*Z # dRda0
    J[2,1] <- -Z   # dRdb0
    J %*% E %*% t(J)
}

get.Ew_Th <- function(w,W,a0,b0,c0){
    ns <- length(W)
    U <- settings('iratio','U238U235')[1]
    E <- diag(c(a0,b0)*w)^2
    J <- matrix(0,3*ns,2)
    J[1:ns,1] <- c0*U*W         # dKda0
    J[(ns+1):(2*ns),2] <- c0*W  # dLdb0
    J %*% E %*% t(J)
}

fixit <- function(x,anchor=list(FALSE,NA)){
    if (x$format %in% c(1,2,3)){
        if (!anchor[[1]]) out <- rep(FALSE,2)
        else if (is.numeric(anchor[[2]])) out <- c(TRUE,FALSE)
        else out <- c(FALSE,TRUE)
    } else if (x$format %in% c(4,5,6,7,8)){
        if (!anchor[[1]]) out <- rep(FALSE,3)
        else if (is.numeric(anchor[[2]])) out <- c(TRUE,FALSE,FALSE)
        else out <- c(FALSE,TRUE,TRUE)
    } else {
        stop('Illegal input format.')
    }
    out
}
