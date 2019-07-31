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
            exterr=exterr,fixed=fixit(x,anchor),lower=c(0,0,0),upper=c(1000,100,100))
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
    D <- c(l$K,l$L)
    out <- D %*% l$omega %*% D
    if (LL){
        k <- length(D)
        detE <- determinant(2*pi*l$omegainv,logarithm=TRUE)$modulus
        out <- -0.5*(out + k*log(2*pi) + detE)
    }
    out
}
LL.lud.3D <- function(ta0b0,x,exterr=FALSE,w=0,LL=FALSE){
    l <- data2ludwig(x,tt=ta0b0[1],a0=ta0b0[2],b0=ta0b0[3],exterr=exterr,w=w)
    D <- c(l$K,l$L,l$M)
    out <- D %*% l$omega %*% D
    if (LL){
        k <- 3*ns
        detE <- determinant(l$omegainv,logarithm=TRUE)$modulus
        out <- -0.5*(out + k*log(2*pi) + detE)
    }
    out
}
LL.lud.Th <- function(ta0b0,x,exterr=FALSE,w=0,LL=FALSE){
    l <- data2ludwig(x,tt=ta0b0[1],a0=ta0b0[2],b0=ta0b0[3],exterr=exterr,w=w)
    D <- c(l$K,l$L,l$M)
    out <- D %*% l$omega %*% D
    if (LL){
        k <- length(D)
        detE <- determinant(l$omegainv,logarithm=TRUE)$modulus
        out <- -0.5*(out + k*log(2*pi) + detE)
    }
    out
}

LL.lud.UPb.gr <- function(ta0b0,x,exterr=FALSE,w=0){
    if (x$format<4){
        return(LL.lud.2D.gr(ta0b0,x,exterr=exterr,w=w))
    } else {
        return(LL.lud.3D.gr(ta0b0,x=x,exterr=exterr,w=w))
    }
}
LL.lud.2D.gr <- function(ta0,x,exterr=FALSE,w=0){
    tt <- ta0[1]
    a0 <- ta0[2]
    U <- settings('iratio','U238U235')[1]
    l <- data2ludwig(x,tt=tt,a0=a0,exterr=exterr,w=w)
    X <- get.U238Pb206.ratios(x)[,1]
    ns <- length(x)
    dRdt <- matrix(0,2*ns,1)
    dRda0 <- matrix(0,2*ns,1)
    D <- mclean(tt=tt,d=x$d)
    # b = D$Pb207U235/U - a0*D$Pb206U238
    dbdt <- D$dPb207U235dt/U - a0*D$dPb206U238dt
    drydt <- -dbdt*(X-l$R) # x = X - R
    dryda0 <- D$Pb206U238*(X-l$R) - 1
    dRdt[(ns+1):(2*ns),1] <- drydt
    dRda0[(ns+1):(2*ns),1] <- dryda0
    Rr <- rbind(R,r)
    dSdt <- 2 * t(R) %*% l$omega %*% dRdt
    dSda0 <- 2 * t(R) %*% l$omega %*% dRda0
    c(dSdt,dSda0)
}
LL.lud.3D.gr <- function(ta0b0,x,exterr=FALSE,w=0){
    tt <- ta0b0[1]
    a0 <- ta0b0[2]
    b0 <- ta0b0[3]
    U <- iratio('U238U235')[1]
    l <- data2ludwig(x,tt=tt,a0=a0,b0=b0,exterr=exterr,w=w)
    D <- mclean(tt=tt,d=x$d)
    dSdt <- dSdx_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='t')
    dSda0 <- dSdx_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='a0')
    dSdb0 <- dSdx_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='b0')
    c(dSdt,dSda0,dSdb0)
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
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    D <- mclean(tt,d=x$d)
    out <- matrix(0,ns+2,ns+2)
    for (i in 1:ns){
        out[i,ns+1] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,D=D,x='y',y='t',i=i)
        out[i,ns+2] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,D=D,x='y',y='a0',i=i)
        for (j in 1:ns){
            out[i,j] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,D=D,x='y',y='y',i=i,j=j)
        }
    }
    out[ns+1,ns+1] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,D=D,x='t',y='t')
    out[ns+2,ns+2] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,D=D,x='a0',y='a0')
    out[ns+1,ns+2] <- d2Sdxdy_2D(l=l,tt=tt,a0=a0,D=D,x='t',y='a0')
    out[ns+1,1:ns] <- out[1:ns,ns+1]
    out[ns+2,1:ns] <- out[1:ns,ns+2]
    out[(ns+2),(ns+1)] <- out[(ns+1),(ns+2)]
    out/2    
}
fisher_lud_2D.old <- function(x,fit){
    tt <- fit$par[1]
    a0 <- fit$par[2]
    l <- data2ludwig_2D(x,tt=tt,a0=a0,w=0,exterr=fit$exterr)
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
    U <- iratio('U238U235')[1]
    D <- mclean(tt=tt,d=x$d)
    ns <- length(x)
    out <- matrix(0,ns+3,ns+3)
    for (i in 1:ns){
        out[i,ns+1] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='z',y='t',i=i)
        out[i,ns+2] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='z',y='a0',i=i)
        out[i,ns+3] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='z',y='b0',i=i)        
        for (j in 1:ns){
            out[i,j] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='z',y='z',i=i,j=j)
        }
    }
    out[ns+1,ns+1] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='t',y='t')
    out[ns+2,ns+2] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='a0',y='a0')
    out[ns+3,ns+3] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='b0',y='b0')
    out[ns+1,ns+2] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='t',y='a0')
    out[ns+1,ns+3] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='t',y='b0')
    out[ns+2,ns+3] <- d2Sdxdy_3D(l=l,tt=tt,a0=a0,b0=b0,U=U,D=D,x='a0',y='b0')
    out[ns+1,1:ns] <- out[1:ns,ns+1]
    out[ns+2,1:ns] <- out[1:ns,ns+2]
    out[ns+3,1:ns] <- out[1:ns,ns+3]
    out[(ns+2),(ns+1)] <- out[(ns+1),(ns+2)]
    out[(ns+3),(ns+1)] <- out[(ns+1),(ns+3)]
    out[(ns+3),(ns+2)] <- out[(ns+2),(ns+3)]
    out/2
}
fisher_lud_Th <- function(x,fit){
    tt <- fit$par[1]
    a0 <- fit$par[2]
    b0 <- fit$par[3]
    l <- data2ludwig_Th(x,tt=tt,a0=a0,b0=b0,w=0,exterr=fit$exterr)
    D <- mclean(tt=tt,d=x$d)
    ns <- length(l$K)
    out <- matrix(0,ns+3,ns+3)
    for (i in 1:ns){
        out[i,ns+1] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='c0',y='t',i=i)
        out[i,ns+2] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='c0',y='a0',i=i)
        out[i,ns+3] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='c0',y='b0',i=i)        
        for (j in 1:ns){
            out[i,j] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='c0',y='c0',i=i,j=j)
        }
    }
    out[ns+1,ns+1] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='t',y='t')
    out[ns+2,ns+2] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='a0',y='a0')
    out[ns+3,ns+3] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='b0',y='b0')
    out[ns+1,ns+2] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='t',y='a0')
    out[ns+1,ns+3] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='t',y='b0')
    out[ns+2,ns+3] <- d2Sdxdy_Th(l=l,tt=tt,a0=a0,b0=b0,D=D,x='a0',y='b0')
    out[ns+1,1:ns] <- out[1:ns,ns+1]
    out[ns+2,1:ns] <- out[1:ns,ns+2]
    out[ns+3,1:ns] <- out[1:ns,ns+3]
    out[(ns+2),(ns+1)] <- out[(ns+1),(ns+2)]
    out[(ns+3),(ns+1)] <- out[(ns+1),(ns+3)]
    out[(ns+3),(ns+2)] <- out[(ns+2),(ns+3)]
    out/2
}

S.lud <- function(l){
    return(l$K%*%l$omega11%*%l$K + l$K%*%l$omega12%*%l$L + l$K%*%l$omega13%*%l$M +
           l$L%*%l$omega21%*%l$K + l$L%*%l$omega22%*%l$L + l$L%*%l$omega23%*%l$M +
           l$M%*%l$omega31%*%l$K + l$M%*%l$omega32%*%l$L + l$M%*%l$omega33%*%l$M)
}
dSdx_3D <- function(l,tt=0,a0=0,b0=0,U=iratio('U238U235')[1],
                    D=mclean(tt=tt,d=diseq()),x='a0'){
    ns <- length(l$z)
    zeros <- rep(0,ns)
    dKdx <- zeros; dLdx <- zeros; dMdx <- zeros
    if (identical(x,'t')){
        dKdx <- l$dKdt
        dLdx <- l$dLdt
        dMdx <- l$dMdt
    } else if (identical(x,'a0')){
        dKdx <- l$dKda0
        dLdx <- l$dLda0
        dMdx <- l$dMda0
    } else if (identical(x,'b0')){
        dKdx <- l$dKdb0
        dLdx <- l$dLdb0
        dMdx <- l$dMdb0
    }
    dSdx(dKdx,dLdx,dMdx,l)
}
dSdx <- function(dKdx,dLdx,dMdx,l){
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
d2Sdxdy_2D <- function(l,tt=0,a0=0,U=iratio('U238U235')[1],
                       D=mclean(tt=tt,d=diseq()),x='a0',y='y',i=1,j=1){
    ns <- length(l$y)
    zeros <- rep(0,ns)
    dKdx <- zeros; dLdx <- zeros
    dKdy <- zeros; dLdy <- zeros
    d2Kdxdy <- zeros; d2Ldxdy <- zeros;
    if (identical(x,'t')){
        dKdx <- rep(-D$dPb207U235dt,ns)
        if (identical(y,'t')) d2Kdxdy <- rep(-D$d2Pb207U235dt2,ns)
    } else if (identical(x,'a0')){
        dKdx <- -U*l$y
    } else if (identical(x,'y')){
        dKdx[i] <- -U*a0
        dLdx[i] <- -1
        if (identical(y,'a0')) d2Kdxdy[j] <- -U
    }
    if (identical(y,'t')){
        dKdy <- rep(-D$dPb207U235dt,ns)
        if (identical(x,'t')){
            d2Kdxdy <- rep(-D$d2Pb207U235dt2,ns)
        }
    } else if (identical(y,'a0')){
        dKdy <- -U*l$y
    } else if (identical(y,'y')){
        dKdy[i] <- -U*a0
        dLdy[i] <- -1
        if (identical(x,'a0')) d2Kdxdy[j] <- -U
    }
    d2Sdxdy(K=l$K,dKdx=dKdx,dKdy=dKdy,d2Kdxdy=d2Kdxdy,
            L=l$L,dLdx=dLdx,dLdy=dLdy,d2Ldxdy=d2Ldxdy,
            M=zeros,dMdx=zeros,dMdy=zeros,d2Mdxdy=zeros,l=l)
}
d2Sdxdy_3D <- function(l,tt=0,a0=0,b0=0,U=iratio('U238U235')[1],
                       D=mclean(tt=tt,d=diseq()),x='a0',y='z',i=1,j=1){
    ns <- length(l$z)
    zeros <- rep(0,ns)
    dKdx <- zeros; dLdx <- zeros; dMdx <- zeros;
    dKdy <- zeros; dLdy <- zeros; dMdy <- zeros;
    d2Kdxdy <- zeros; d2Ldxdy <- zeros; d2Mdxdy <- zeros;
    if (identical(x,'t')){
        dKdx <- rep(-D$dPb207U235dt,ns)
        dLdx <- rep(-D$dPb206U238dt,ns)
        if (identical(y,'t')){
            d2Kdxdy <- rep(-D$d2Pb207U235dt2,ns)
            d2Ldxdy <- rep(-D$d2Pb206U238dt2,ns)
        }
    } else if (identical(x,'a0')){
        dLdx <- -l$z
    } else if (identical(x,'b0')){
        dKdx <- -U*l$z
    } else if (identical(x,'z')){
        dKdx[i] <- -U*b0
        dLdx[i] <- -a0
        dMdx[i] <- -1
        if (identical(y,'a0')) d2Ldxdy[j] <- -1
        else if (identical(y,'b0')) d2Kdxdy[j] <- -U
    }
    if (identical(y,'t')){
        dKdy <- rep(-D$dPb207U235dt,ns)
        dLdy <- rep(-D$dPb206U238dt,ns)
        if (identical(x,'t')){
            d2Kdxdy <- rep(-D$d2Pb207U235dt2,ns)
            d2Ldxdy <- rep(-D$d2Pb206U238dt2,ns)
        }
    } else if (identical(y,'a0')){
        dLdy <- -l$z
    } else if (identical(y,'b0')){
        dKdy <- -U*l$z
    } else if (identical(y,'z')){
        dKdy[j] <- -U*b0
        dLdy[j] <- -a0
        dMdy[j] <- -1
        if (identical(x,'a0')) d2Ldxdy[i] <- -1
        else if (identical(x,'b0')) d2Kdxdy[i] <- -U
    }
    d2Sdxdy(l$K,dKdx,dKdy,d2Kdxdy,
            l$L,dLdx,dLdy,d2Ldxdy,
            l$M,dMdx,dMdy,d2Mdxdy,l=l)
}
d2Sdxdy_Th <- function(l,tt=0,a0=0,b0=0,D=mclean(tt=tt,d=diseq()),x='a0',y='c0',i=1,j=1){
    U <- iratio('U238U235')[1]
    l32 <- lambda('Th232')[1]
    ns <- length(l$K)
    dKdx <- rep(0,ns); dLdx <- rep(0,ns); dMdx <- rep(0,ns);
    dKdy <- rep(0,ns); dLdy <- rep(0,ns); dMdy <- rep(0,ns);
    d2Kdxdy <- rep(0,ns); d2Ldxdy <- rep(0,ns); d2Mdxdy <- rep(0,ns);
    if (identical(x,'t')){
        dKdx <- rep(-D$dPb207U235dt,ns)
        dLdx <- rep(-D$dPb206U238dt,ns)
        dMdx <- rep(-exp(l32*tt)*l32,ns)
        if (identical(y,'t')){
            d2Kdxdy <- rep(-D$d2Pb207U235dt2,ns)
            d2Ldxdy <- rep(-D$d2Pb206U238dt2,ns)
            d2Mdxdy <- rep(-exp(l32*tt)*l32^2,ns)
        }
    } else if (identical(x,'a0')){
        dLdx <- -l$c0*l$W
    } else if (identical(x,'b0')){
        dKdx <- -U*l$c0*l$W
    } else if (identical(x,'c0')){
        dKdx[i] <- -U*b0*l$W[i]
        dLdx[i] <- -a0*l$W[i]
        dMdx[i] <- -1
        if (identical(y,'a0')) d2Ldxdy[j] <- -l$W[j]
        else if (identical(y,'b0')) d2Kdxdy[j] <- -U*l$W[j]
    }
    if (identical(y,'t')){
        dKdy <- rep(-D$dPb207U235dt,ns)
        dLdy <- rep(-D$dPb206U238dt,ns)
        dMdy <- rep(-exp(l32*tt)*l32,ns)
        if (identical(x,'t')){
            d2Kdxdy <- rep(-D$d2Pb207U235dt2,ns)
            d2Ldxdy <- rep(-D$d2Pb206U238dt2,ns)
            d2Mdxdy <- rep(-exp(l32*tt)*l32^2,ns)
        }
    } else if (identical(y,'a0')){
        dLdy <- -l$c0*l$W
    } else if (identical(y,'b0')){
        dKdy <- -U*l$c0*l$W
    } else if (identical(y,'c0')){
        dKdy[j] <- -U*b0*l$W[i]
        dLdy[j] <- -a0*l$W[i]
        dMdy[j] <- -1
        if (identical(x,'a0')) d2Ldxdy[i] <- -l$W[i]
        else if (identical(x,'b0')) d2Kdxdy[i] <- -U*l$W[i]
    }
    d2Sdxdy(l$K,dKdx,dKdy,d2Kdxdy,
            l$L,dLdx,dLdy,d2Ldxdy,
            l$M,dMdx,dMdy,d2Mdxdy,l)
}
d2Sdxdy <- function(K,dKdx,dKdy,d2Kdxdy,
                    L,dLdx,dLdy,d2Ldxdy,
                    M,dMdx,dMdy,d2Mdxdy,l){
    return(K%*%l$omega11%*%d2Kdxdy + dKdy%*%l$omega11%*%dKdx +
           dKdx%*%l$omega11%*%dKdy + d2Kdxdy%*%l$omega11%*%K +
           K%*%l$omega12%*%d2Ldxdy + dKdy%*%l$omega12%*%dLdx +
           dKdx%*%l$omega12%*%dLdy + d2Kdxdy%*%l$omega12%*%L +
           K%*%l$omega13%*%d2Mdxdy + dKdy%*%l$omega13%*%dMdx +
           dKdx%*%l$omega13%*%dMdy + d2Kdxdy%*%l$omega13%*%M +
           L%*%l$omega21%*%d2Kdxdy + dLdy%*%l$omega21%*%dKdx +
           dLdx%*%l$omega21%*%dKdy + d2Ldxdy%*%l$omega21%*%K +
           L%*%l$omega22%*%d2Ldxdy + dLdy%*%l$omega22%*%dLdx +
           dLdx%*%l$omega22%*%dLdy + d2Ldxdy%*%l$omega22%*%L +
           L%*%l$omega23%*%d2Mdxdy + dLdy%*%l$omega23%*%dMdx +
           dLdx%*%l$omega23%*%dMdy + d2Ldxdy%*%l$omega23%*%M +
           M%*%l$omega31%*%d2Kdxdy + dMdy%*%l$omega31%*%dKdx +
           dMdx%*%l$omega31%*%dKdy + d2Mdxdy%*%l$omega31%*%K +
           M%*%l$omega32%*%d2Ldxdy + dMdy%*%l$omega32%*%dLdx +
           dMdx%*%l$omega32%*%dLdy + d2Mdxdy%*%l$omega32%*%L +
           M%*%l$omega33%*%d2Mdxdy + dMdy%*%l$omega33%*%dMdx +
           dMdx%*%l$omega33%*%dMdy + d2Mdxdy%*%l$omega33%*%M)
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
    U <- iratio('U238U235')[1]
    ns <- length(x)
    E <- matrix(0,2*ns+6,2*ns+6)
    J <- matrix(0,2*ns,2*ns+6)
    J[1:(2*ns),1:(2*ns)] <- diag(2*ns)
    X <- rep(0,ns)
    Y <- rep(0,ns)
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    for (i in 1:ns){
        wd <- wetherill(x,i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        E[c(i,ns+i),c(i,ns+i)] <- wd$cov
        J[i,ns+i] <- -a0*U                   # dKdY
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
    out$y <- Y - D$Pb206U238 - out$L
    out$K <- X - D$Pb207U235 - a0*U*out$y
    out
}
# rederived from Ludwig (1998):
data2ludwig_3D <- function(x,tt,a0,b0,w=0,exterr=FALSE){
    U <- iratio('U238U235')[1]
    ns <- length(x)
    E <- matrix(0,3*ns+6,3*ns+6)
    J <- matrix(0,3*ns,3*ns+6)
    J[1:(3*ns),1:(3*ns)] <- diag(3*ns)
    X <- rep(0,ns)
    Y <- rep(0,ns)
    Z <- rep(0,ns)
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
    Winv <- solve(W)
    out <- list(omega=omega,omegainv=ED,omega11=o11,omega12=o12,
                omega13=o13,omega21=o21,omega22=o22,omega23=o23,
                omega31=o31,omega32=o32,omega33=o33)
    out$M <- as.vector(Winv%*%V)
    out$z <- as.vector(Z - out$M)
    out$K <- as.vector(X - U*b0*out$z - x75)
    out$L <- as.vector(Y - a0*out$z - x68)
    out$dMdt <- rep(0,ns)
    out$dKdt <- rep(-D$dPb207U235dt,ns)
    out$dLdt <- rep(-D$dPb206U238dt,ns)
    out$dMda0 <- rep(0,ns)
    out$dKda0 <- rep(0,ns)
    out$dLda0 <- as.vector(-out$z)
    out$dMdb0 <- rep(0,ns)
    out$dKdb0 <- as.vector(-U*out$z)
    out$dLdb0 <- rep(0,ns)
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
    X <- rep(0,ns)
    Y <- rep(0,ns)
    Z <- rep(0,ns)
    W <- rep(0,ns)
    K0 <- rep(0,ns)
    L0 <- rep(0,ns)
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
        E[c(i,ns+i,2*ns+i,3*ns+i),c(i,ns+i,2*ns+i,3*ns+i)] <- wd$cov
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
    ED <- J%*%E%*%t(J) + get.Ew_Th(w=w,W=x$x[,'Th232U238'],a0=a0,b0=b0,c0=c0)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    i3 <- (2*ns+1):(3*ns)
    omega <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                             DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                             GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    o11 <- omega[i1,i1]; o12 <- omega[i1,i2]; o13 <- omega[i1,i3]
    o21 <- omega[i2,i1]; o22 <- omega[i2,i2]; o23 <- omega[i2,i3]
    o31 <- omega[i3,i1]; o32 <- omega[i3,i2]; o33 <- omega[i3,i3]
    A <- t(K0%*%(o11+t(o11))*U*b0*W + K0%*%(o12+t(o21))*a0*W + K0%*%(o13+t(o31)) +
           L0%*%(o21+t(o12))*U*b0*W + L0%*%(o22+t(o22))*a0*W + L0%*%(o23+t(o32)))
    B <- -(U*b0*W*(o11+t(o11))*U*b0*W + a0*W*(o22+t(o22))*a0*W + (o33+t(o33)) +
           U*b0*W*(o12+t(o21))*a0*W + a0*W*(o12+t(o21))*U*b0*W + U*b0*W*(o13+t(o31)) +
           (o13+t(o31))*U*b0*W + a0*W*(o23+t(o32)) )
    out <- list(omega11=o11,omega12=o12,omega13=o13,
                omega21=o21,omega22=o22,omega23=o23,omega31=o31,
                omega32=o32,omega33=o33,omega=omega,omegainv=ED)
    out$M <- as.vector(solve(B,A))
    out$c0 <- as.vector(Z - out$M - x82)
    out$K <- as.vector(X - out$c0*b0*U*W - x75)
    out$L <- as.vector(Y - out$c0*a0*W - x68)
    out$dKdt <- rep(-D$Pb207U235,ns)
    out$dLdt <- rep(-D$Pb206U238,ns)
    out$dMdt <- rep(-l2*exp(l2*tt),ns)
    out$dKda0 <- rep(0,ns)
    out$dLda0 <- -out$c0*W
    out$dMda0 <- rep(0,ns)
    out$dKdb0 <- -out$c0*U*W
    out$dLdb0 <- rep(0,ns)
    out$dMdb0 <- rep(0,ns)
    out
}

get.Ew_2D <- function(w,Y,a0){
    if (w>0){
        ns <- length(Y)
        U <- iratio('U238U235')[1]
        sa0 <- (a0*w)^2
        J <- matrix(0,2*ns,1)
        J[1:ns,1] <- -U*Y
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
        J[1:ns,2] <- -U*Z   # dRdb0
        J[(ns+1):(2*ns),1] <- -Z # drda0
        out <- J %*% E %*% t(J)
    } else {
        out <- 0
    }
    out
}
get.Ew_Th <- function(w,W,a0,b0,c0){
    if (w>0){
        ns <- length(W)
        U <- settings('iratio','U238U235')[1]
        E <- diag(c(a0,b0)*w)^2
        J <- matrix(0,3*ns,2)
        J[1:ns,1] <- -c0*U*W         # dKda0
        J[(ns+1):(2*ns),2] <- -c0*W  # dLdb0
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
