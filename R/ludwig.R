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
ludwig.UPb <- function(x,exterr=FALSE,alpha=0.05,model=1,anchor=list(FALSE,NA),...){
    fit <- get.ta0b0(x,exterr=exterr,model=model,anchor=anchor,...)
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

get.ta0b0 <- function(x,exterr=FALSE,model=1,anchor=list(FALSE,NA),...){
    init <- get.ta0b0.model2(x,anchor=anchor)
    if (model==1)
        out <- get.ta0b0.model1(x,init=init,exterr=exterr,anchor=anchor,...)
    else if (model==2)
        out <- list(par=init,w=0)
    else if (model==3)
        out <- get.ta0b0.model3(x,init=init,exterr=exterr,anchor=anchor,...)
    out$model <- model
    out$exterr <- exterr
    out
}
get.ta0b0.model1 <- function(x,init,exterr=FALSE,anchor=list(FALSE,NA),...){
    out <- fit_ludwig_discordia(x,init=init,w=0,exterr=exterr,anchor=anchor,...)
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
get.ta0b0.model3 <- function(x,init,exterr=FALSE,anchor=list(FALSE,NA),...){
    out <- list(w=0,par=init)
#    for (i in 1:5){ # loop for more accurate but slower and more unstable results
    out <- fit_ludwig_discordia(x,init=out$par,w=out$w,
                                exterr=exterr,anchor=anchor,...)
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
    optifix(parms=init,fn=LL.lud.UPb,gr=LL.lud.UPb.gr,method="L-BFGS-B",
            x=x,w=w,exterr=exterr,fixed=fixit(x,anchor),
            lower=c(0,0,0),upper=c(1000,100,100),...)
}

LL.lud.disp <- function(w,x,ta0b0,exterr=FALSE,anchor=list(FALSE,NA)){
    # these two lines produce slightly more accurate but much slower results:
    # fit <- fit_ludwig_discordia(x,init=ta0b0,w=w,exterr=exterr,anchor=anchor)
    # LL.lud.UPb(ta0b0=fit$par,x=x,exterr=exterr,w=w,LL=TRUE)
    LL.lud.UPb(ta0b0=ta0b0,x=x,exterr=exterr,w=w,LL=TRUE)
}
LL.lud.UPb <- function(ta0b0,x,exterr=FALSE,w=0,LL=FALSE){
    if (x$format < 4){
        l <- data2ludwig(x,tt=ta0b0[1],a0=ta0b0[2],exterr=exterr,w=w)
        D <- c(l$K,l$L)
    } else {
        l <- data2ludwig(x,tt=ta0b0[1],a0=ta0b0[2],b0=ta0b0[3],exterr=exterr,w=w)
        D <- c(l$K,l$L,l$M)
    }
    out <- D %*% l$omega %*% D
    if (LL){
        k <- length(D)
        detE <- determinant(2*pi*l$omegainv,logarithm=TRUE)$modulus
        out <- -0.5*(out + k*log(2*pi) + detE)
    }
    out    
}

LL.lud.UPb.gr <- function(ta0b0,x,exterr=FALSE,w=0){
    tt <- ta0b0[1]
    a0 <- ta0b0[2]
    if (x$format<4){
        l <- data2ludwig(x,tt=tt,a0=a0,exterr=exterr,w=w)
        dSdt <- dSdx_2D(l=l,tt=tt,a0=a0,x='t')
        dSda0 <- dSdx_2D(l=l,tt=tt,a0=a0,x='a0')
        out <- c(dSdt,dSda0)
    } else {
        b0 <- ta0b0[3]
        l <- data2ludwig(x,tt=tt,a0=a0,b0=b0,exterr=exterr,w=w)
        dSdt <- dSdx_3D(l=l,tt=tt,a0=a0,b0=b0,x='t')
        dSda0 <- dSdx_3D(l=l,tt=tt,a0=a0,b0=b0,x='a0')
        dSdb0 <- dSdx_3D(l=l,tt=tt,a0=a0,b0=b0,x='b0')
        out <- c(dSdt,dSda0,dSdb0)
    }
    out
}

fisher.lud <- function(x,fit,exterr=TRUE,anchor=list(FALSE,NA),...){
    ns <- length(x)
    i1 <- 1:ns
    if (x$format < 4){
        fish <- fisher_lud_2D(x,fit)
        i2 <- (ns+1):(ns+2)
    } else {
        fish <- fisher_lud_3D(x,fit)
        i2 <- (ns+1):(ns+3)
    }
    anchorfish(AA=fish[i1,i1],BB=fish[i1,i2],
               CC=fish[i2,i1],DD=fish[i2,i2],anchor=anchor)
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
    ns <- length(x)
    out <- matrix(0,ns+2,ns+2)
    for (i in 1:ns){
        out[i,ns+1] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,x='c0',y='t',i=i)
        out[i,ns+2] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,x='c0',y='a0',i=i)
        for (j in 1:ns){
            out[i,j] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,x='c0',y='c0',i=i,j=j)
        }
    }
    out[ns+1,ns+1] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,x='t',y='t')
    out[ns+2,ns+2] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,x='a0',y='a0')
    out[ns+1,ns+2] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,x='t',y='a0')
    out[ns+2,ns+1] <- out[ns+1,ns+2]
    out[ns+1,1:ns] <- out[1:ns,ns+1]
    out[ns+2,1:ns] <- out[1:ns,ns+2]
    out[(ns+2),(ns+1)] <- out[(ns+1),(ns+2)]
    out
}
fisher_lud_3D <- function(x,fit){
    tt <- fit$par[1]
    a0 <- fit$par[2]
    b0 <- fit$par[3]
    if (x$format>6){
        l <- data2ludwig_Th(x,tt=tt,a0=a0,b0=b0,w=0,exterr=fit$exterr)
    } else {
        l <- data2ludwig_3D(x,tt=tt,a0=a0,b0=b0,w=0,exterr=fit$exterr)
    }
    ns <- length(x)
    out <- matrix(0,ns+3,ns+3)
    for (i in 1:ns){
        out[i,ns+1] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='c0',y='t',i=i)
        out[i,ns+2] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='c0',y='a0',i=i)
        out[i,ns+3] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='c0',y='b0',i=i)        
        for (j in 1:ns){
            out[i,j] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='c0',y='c0',i=i,j=j)
        }
    }
    out[ns+1,ns+1] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='t',y='t')
    out[ns+2,ns+2] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='a0',y='a0')
    out[ns+3,ns+3] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='b0',y='b0')
    out[ns+1,ns+2] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='t',y='a0')
    out[ns+1,ns+3] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='t',y='b0')
    out[ns+2,ns+3] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,x='a0',y='b0')
    out[ns+1,1:ns] <- out[1:ns,ns+1]
    out[ns+2,1:ns] <- out[1:ns,ns+2]
    out[ns+3,1:ns] <- out[1:ns,ns+3]
    out[(ns+2),(ns+1)] <- out[(ns+1),(ns+2)]
    out[(ns+3),(ns+1)] <- out[(ns+1),(ns+3)]
    out[(ns+3),(ns+2)] <- out[(ns+2),(ns+3)]
    out
}

dSdx_2D <- function(l,tt=0,a0=0,x='a0'){
    if (identical(x,'t')){
        dKdx <- l$dKdt
        dLdx <- l$dLdt
    } else if (identical(x,'a0')){
        dKdx <- l$dKda0
        dLdx <- 0*l$K
    }
    return(l$K%*%l$omega11%*%dKdx + dKdx%*%l$omega11%*%l$K +
           l$K%*%l$omega12%*%dLdx + dKdx%*%l$omega12%*%l$L +
           l$L%*%l$omega21%*%dKdx + dLdx%*%l$omega21%*%l$K +
           l$L%*%l$omega22%*%dLdx + dLdx%*%l$omega22%*%l$L)
}
dSdx_3D <- function(l,tt=0,a0=0,b0=0,x='a0'){
    zeros <- 0*l$K
    if (identical(x,'t')){
        dKdx <- l$dKdt
        dLdx <- l$dLdt
        dMdx <- l$dMdt
    } else if (identical(x,'a0')){
        dKdx <- zeros
        dLdx <- l$dLda0
        dMdx <- zeros
    } else if (identical(x,'b0')){
        dKdx <- l$dKdb0
        dLdx <- zeros
        dMdx <- zeros
    }
    return(l$K%*%l$omega11%*%dKdx + dKdx%*%l$omega11%*%l$K +
           l$K%*%l$omega12%*%dLdx + dKdx%*%l$omega12%*%l$L +
           l$K%*%l$omega13%*%dMdx + dKdx%*%l$omega13%*%l$M +
           l$L%*%l$omega21%*%dKdx + dLdx%*%l$omega21%*%l$K +
           l$L%*%l$omega22%*%dLdx + dLdx%*%l$omega22%*%l$L +
           l$L%*%l$omega23%*%dMdx + dLdx%*%l$omega23%*%l$M +
           l$M%*%l$omega31%*%dKdx + dMdx%*%l$omega31%*%l$K +
           l$M%*%l$omega32%*%dLdx + dMdx%*%l$omega32%*%l$L +
           l$M%*%l$omega33%*%dMdx + dMdx%*%l$omega33%*%l$M)
}
d2Sdxdy_2D <- function(l,tt=0,a0=0,x='a0',y='c0',i=1,j=1){
    ns <- length(l$K)
    zeros <- rep(0,ns)
    dKdx <- zeros; dLdx <- zeros
    dKdy <- zeros; dLdy <- zeros
    d2Kdxdy <- zeros; d2Ldxdy <- zeros;
    if (identical(x,'t')){
        dKdx <- l$dKdt
        dLdx <- l$dLdt
        if (identical(y,'t')) d2Kdxdy <- l$d2Kdt2
    } else if (identical(x,'a0')){
        dKdx <- l$dKda0
    } else if (identical(x,'c0')){
        dKdx[i] <- l$dKdc0[i]
        dLdx[i] <- l$dLdc0[i]
        if (identical(y,'a0')) d2Kdxdy[i] <- l$d2Kdc0da0[i]
    } else {
        stop('Invalid option.')
    }
    if (identical(y,'t')){
        dKdy <- l$dKdt
        dLdy <- l$dLdt
    } else if (identical(y,'a0')){
        dKdy <- l$dKda0
    } else if (identical(y,'c0')){
        dKdy[j] <- l$dKdc0[j]
        dLdy[j] <- l$dLdc0[j]
        if (identical(x,'a0')) d2Kdxdy[j] <- l$d2Kdc0da0[j]
    } else {
        stop('Invalid option.')
    }
    return(l$K%*%l$omega11%*%d2Kdxdy + dKdy%*%l$omega11%*%dKdx +
           dKdx%*%l$omega11%*%dKdy + d2Kdxdy%*%l$omega11%*%l$K +
           l$K%*%l$omega12%*%d2Ldxdy + dKdy%*%l$omega12%*%dLdx +
           dKdx%*%l$omega12%*%dLdy + d2Kdxdy%*%l$omega12%*%l$L +
           l$L%*%l$omega21%*%d2Kdxdy + dLdy%*%l$omega21%*%dKdx +
           dLdx%*%l$omega21%*%dKdy + d2Ldxdy%*%l$omega21%*%l$K +
           l$L%*%l$omega22%*%d2Ldxdy + dLdy%*%l$omega22%*%dLdx +
           dLdx%*%l$omega22%*%dLdy + d2Ldxdy%*%l$omega22%*%l$L)
}
d2Sdxdy_3D <- function(l,tt=0,a0=0,b0=0,x='a0',y='c0',i=1,j=1){
    ns <- length(l$K)
    zeros <- rep(0,ns)
    dKdx <- zeros; dLdx <- zeros; dMdx <- zeros;
    dKdy <- zeros; dLdy <- zeros; dMdy <- zeros;
    d2Kdxdy <- zeros; d2Ldxdy <- zeros; d2Mdxdy <- zeros;
    if (identical(x,'t')){
        dKdx <- l$dKdt
        dLdx <- l$dLdt
        dMdx <- l$dMdt
        if (identical(y,'t')){
            d2Kdxdy <- l$d2Kdt2
            d2Ldxdy <- l$d2Ldt2
            d2Mdxdy <- l$d2Mdt2
        }
    } else if (identical(x,'a0')){
        dLdx <- l$dLda0
    } else if (identical(x,'b0')){
        dKdx <- l$dKdb0
    } else if (identical(x,'c0')){
        dKdx[i] <- l$dKdc0[i]
        dLdx[i] <- l$dLdc0[i]
        dMdx[i] <- l$dMdc0[i]
        if (identical(y,'a0')) d2Ldxdy[i] <- l$d2Ldc0da0[i]
        else if (identical(y,'b0')) d2Kdxdy[i] <- l$d2Kdc0db0[i]
    } else {
        stop('Invalid option.')
    }
    if (identical(y,'t')){
        dKdy <- l$dKdt
        dLdy <- l$dLdt
        dMdy <- l$dMdt
    } else if (identical(y,'a0')){
        dLdy <- l$dLda0
    } else if (identical(y,'b0')){
        dKdy <- l$dKdb0
    } else if (identical(y,'c0')){
        dKdy[j] <- l$dKdc0[j]
        dLdy[j] <- l$dLdc0[j]
        dMdy[j] <- l$dMdc0[j]
        if (identical(x,'a0')) d2Ldxdy[j] <- l$d2Ldc0da0[j]
        else if (identical(x,'b0')) d2Kdxdy[j] <- l$d2Kdc0db0[j]
    } else {
        stop('Invalid option.')
    }
    return(l$K%*%l$omega11%*%d2Kdxdy + dKdy%*%l$omega11%*%dKdx +
           dKdx%*%l$omega11%*%dKdy + d2Kdxdy%*%l$omega11%*%l$K +
           l$K%*%l$omega12%*%d2Ldxdy + dKdy%*%l$omega12%*%dLdx +
           dKdx%*%l$omega12%*%dLdy + d2Kdxdy%*%l$omega12%*%l$L +
           l$K%*%l$omega13%*%d2Mdxdy + dKdy%*%l$omega13%*%dMdx +
           dKdx%*%l$omega13%*%dMdy + d2Kdxdy%*%l$omega13%*%l$M +
           l$L%*%l$omega21%*%d2Kdxdy + dLdy%*%l$omega21%*%dKdx +
           dLdx%*%l$omega21%*%dKdy + d2Ldxdy%*%l$omega21%*%l$K +
           l$L%*%l$omega22%*%d2Ldxdy + dLdy%*%l$omega22%*%dLdx +
           dLdx%*%l$omega22%*%dLdy + d2Ldxdy%*%l$omega22%*%l$L +
           l$L%*%l$omega23%*%d2Mdxdy + dLdy%*%l$omega23%*%dMdx +
           dLdx%*%l$omega23%*%dMdy + d2Ldxdy%*%l$omega23%*%l$M +
           l$M%*%l$omega31%*%d2Kdxdy + dMdy%*%l$omega31%*%dKdx +
           dMdx%*%l$omega31%*%dKdy + d2Mdxdy%*%l$omega31%*%l$K +
           l$M%*%l$omega32%*%d2Ldxdy + dMdy%*%l$omega32%*%dLdx +
           dMdx%*%l$omega32%*%dLdy + d2Mdxdy%*%l$omega32%*%l$L +
           l$M%*%l$omega33%*%d2Mdxdy + dMdy%*%l$omega33%*%dMdx +
           dMdx%*%l$omega33%*%dMdy + d2Mdxdy%*%l$omega33%*%l$M)
}

data2ludwig <- function(x,tt,a0,b0=0,g0=rep(0,length(x)),exterr=FALSE,w=0,...){
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
    OI <- J%*%E%*%t(J) + get.Ew_2D(w=w,Y=Y,a0=a0)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    O <- blockinverse(AA=OI[i1,i1],BB=OI[i1,i2],
                      CC=OI[i2,i1],DD=OI[i2,i2],doall=TRUE)
    K0 <- X - D$Pb207U235 - a0*U*Y + a0*U*D$Pb206U238
    A <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*a0*U + K0%*%(O[i1,i2]+t(O[i2,i1])))
    B <- -(a0*U*(O[i1,i1]+t(O[i1,i1]))*a0*U + (O[i2,i2]+t(O[i2,i2])) +
           a0*U*(O[i1,i2]+t(O[i1,i2])) + (O[i2,i1]+t(O[i2,i1]))*a0*U)
    out <- list(omega11=O[i1,i1],omega12=O[i1,i2],omega21=O[i2,i1],
                omega22=O[i2,i2],omega=O,omegainv=OI)
    out$L <- as.vector(solve(B,A))
    out$c0 <- Y - D$Pb206U238 - out$L
    out$K <- X - D$Pb207U235 - a0*U*out$c0
    out$dKdt <- rep(-D$dPb207U235dt,ns)
    out$dLdt <- rep(-D$dPb206U238dt,ns)
    out$dKda0 <- -U*out$c0
    out$dKdc0 <- rep(-U*a0,ns)
    out$dLdc0 <- rep(-1,ns)
    out$d2Kdt2 <- rep(-D$d2Pb207U235dt2,ns)
    out$d2Ldt2 <- rep(-D$d2Pb206U238dt2,ns)
    out$d2Kdc0da0 <- rep(-U,ns)
    out
}
# rederived from Ludwig (1998):
data2ludwig_3D <- function(x,tt,a0,b0,w=0,exterr=FALSE){
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
    ED <- J%*%E%*%t(J) + get.Ew_3D(w=w,Z=Z,a0=a0,b0=b0)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    i3 <- (2*ns+1):(3*ns)
    omega <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                             DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                             GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    o11 <- omega[i1,i1]; o12 <- omega[i1,i2]; o13 <- omega[i1,i3]
    o21 <- omega[i2,i1]; o22 <- omega[i2,i2]; o23 <- omega[i2,i3]
    o31 <- omega[i3,i1]; o32 <- omega[i3,i2]; o33 <- omega[i3,i3]
    K0 <- X - U*b0*Z - D$Pb207U235
    L0 <- Y - a0*Z - D$Pb206U238
    V <- t(K0%*%(o11+t(o11))*U*b0 + L0%*%(o12+t(o21))*U*b0 + K0%*%(o12+t(o21))*a0 +
           L0%*%(o22+t(o22))*a0 + K0%*%(o13+t(o31)) + L0%*%(o23+t(o32)))
    W <- -(U*b0*(o11+t(o11))*U*b0 + U*b0*(o12+t(o12))*a0 + U*b0*(o13+t(o13)) +
           a0*(o21+t(o21))*U*b0 + a0*(o22+t(o22))*a0 + a0*(o23+t(o23)) +
           (o31+t(o31))*U*b0 + (o32+t(o32))*a0 + (o33+t(o33)))
    out <- list(omega=omega,omegainv=ED,omega11=o11,omega12=o12,
                omega13=o13,omega21=o21,omega22=o22,omega23=o23,
                omega31=o31,omega32=o32,omega33=o33)
    out$M <- as.vector(solve(W,V))
    out$c0 <- as.vector(Z - out$M)
    out$K <- as.vector(X - U*b0*out$c0 - x75)
    out$L <- as.vector(Y - a0*out$c0 - x68)
    out$dKdt <- rep(-D$dPb207U235dt,ns)
    out$dLdt <- rep(-D$dPb206U238dt,ns)
    out$dMdt <- zeros
    out$dKdb0 <- as.vector(-U*out$c0)
    out$dLda0 <- as.vector(-out$c0)
    out$dKdc0 <- rep(-U*b0,ns)
    out$dLdc0 <- rep(-a0,ns)
    out$dMdc0 <- rep(-1,ns)
    out$d2Kdt2 <- rep(-D$d2Pb207U235dt2,ns)
    out$d2Ldt2 <- rep(-D$d2Pb206U238dt2,ns)
    out$d2Mdt2 <- zeros
    out$d2Ldc0da0 <- rep(-1,ns)
    out$d2Kdc0db0 <- rep(-U,ns)
    out
}
data2ludwig_Th <- function(x,tt,a0,b0,w=0,exterr=FALSE){
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
    ED <- J%*%E%*%t(J) + get.Ew_Th(w=w,W=W,Z=Z,x82=x82,a0=10,b0=b0)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    i3 <- (2*ns+1):(3*ns)
    omega <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                             DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                             GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    o11 <- omega[i1,i1]; o12 <- omega[i1,i2]; o13 <- omega[i1,i3]
    o21 <- omega[i2,i1]; o22 <- omega[i2,i2]; o23 <- omega[i2,i3]
    o31 <- omega[i3,i1]; o32 <- omega[i3,i2]; o33 <- omega[i3,i3]
    A <- t(K0%*%(o11+t(o11))*U*b0*W + K0%*%(o12+t(o21))*a0*W +
           K0%*%(o13+t(o31)) + L0%*%(o21+t(o12))*U*b0*W +
           L0%*%(o22+t(o22))*a0*W + L0%*%(o23+t(o32)))
    B <- -(U*b0*W*(o11+t(o11))*U*b0*W + a0*W*(o22+t(o22))*a0*W + (o33+t(o33)) +
           U*b0*W*(o12+t(o21))*a0*W + a0*W*(o12+t(o21))*U*b0*W +
           U*b0*W*(o13+t(o31)) + (o13+t(o31))*U*b0*W +
           a0*W*(o23+t(o32)) + (o23+t(o32))*a0*W )
    out <- list(omega11=o11,omega12=o12,omega13=o13,
                omega21=o21,omega22=o22,omega23=o23,omega31=o31,
                omega32=o32,omega33=o33,omega=omega,omegainv=ED)
    out$M <- as.vector(solve(B,A))
    out$c0 <- as.vector(Z - out$M - x82)
    out$K <- as.vector(X - out$c0*b0*U*W - x75)
    out$L <- as.vector(Y - out$c0*a0*W - x68)
    out$dKdt <- rep(-D$dPb207U235dt,ns)
    out$dLdt <- rep(-D$dPb206U238dt,ns)
    out$dMdt <- rep(-exp(l2*tt)*l2,ns)
    out$dLda0 <- -out$c0*W
    out$dKdb0 <- -out$c0*U*W
    out$dKdc0 <- -b0*U*W
    out$dLdc0 <- -a0*W
    out$dMdc0 <- rep(-1,ns)
    out$d2Kdt2 <- rep(-D$d2Pb207U235dt2,ns)
    out$d2Ldt2 <- rep(-D$d2Pb206U238dt2,ns)
    out$d2Mdt2 <- rep(-exp(l2*tt)*l2^2,ns)
    out$d2Ldc0da0 <- -W
    out$d2Kdc0db0 <- -U*W
    out
}

get.Ew_2D <- function(w,Y,a0){
    if (w>0){
        ns <- length(Y)
        U <- iratio('U238U235')[1]
        sa0 <- (a0*w)^2
        J <- matrix(0,2*ns,1)
        J[1:ns,1] <- -U*Y               # dK0da0
        out <- J %*% sa0 %*% t(J)
    } else {
        out <- 0
    }
    out
}
get.Ew_3D <- function(w,Z,a0,b0){
    if (w>0){
        ns <- length(Z)
        U <- iratio('U238U235')[1]
        E <- diag(c(a0,b0)*w)^2
        J <- matrix(0,3*ns,2)
        J[1:ns,2] <- -U*Z               # dK0db0
        J[(ns+1):(2*ns),1] <- -Z        # dL0da0
        out <- J %*% E %*% t(J)
    } else {
        out <- 0
    }
    out
}
get.Ew_Th <- function(w,W,Z,x82,a0,b0){
    if (w>0){
        ns <- length(W)
        U <- iratio('U238U235')[1]
        E <- diag(c(a0,b0)*w)^2
        J <- matrix(0,3*ns,2)
        J[1:ns,1] <- (x82-Z)*U*W        # dK0db0
        J[(ns+1):(2*ns),2] <- (x82-Z)*W # dL0da0
        out <- J %*% E %*% t(J)
    } else {
        out <- 0
    }
    out
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
