#' Linear regression of U-Pb data with correlated errors, taking
#' into account decay constant uncertainties.
#'
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
#' @param ... optional arguments
#
# @param x a \eqn{3n}-element vector \eqn{[X Y Z]}, where \eqn{X},
#     \eqn{Y} and \eqn{Z} are three \eqn{n}-element vectors of
#     (isotopic ratio) values.
# @param covmat a \eqn{[3n x 3n]}-element covariance matrix of
#     \code{x}
#
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
    mswd <- mswd.lud(fit$par,x=x,anchor=anchor)
    out <- c(out,mswd)
    if (model==3){
        out$w <- c(fit$w,profile_LL_discordia_disp(fit,x=x,alpha=alpha))
        names(out$w) <- c('s','ll','ul')
    }
    if (x$format<4) parnames <- c('t[l]','76i')
    else parnames <- c('t','64i','74i')
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
    if (x$format<4 && anchored){
        out$df <- ns-1
    } else if (x$format<4) {
        out$df <- ns-2
    } else if (anchored && tanchored){
        out$df <- 2*ns-2
    } else if (anchored){
        out$df <- 2*ns-1
    } else {
        out$df <- 2*ns-3
    }
    out$mswd <- as.vector(SS/out$df)
    out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
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
    xy <- data2york(x,option=2)[,c('X','Y')]
    if (!anchor[[1]]) {
        xyfit <- stats::lm(xy[,'Y'] ~ xy[,'X'])
        intercept <- xyfit$coef[1]
        slope <- xyfit$coef[2]
        ta0b0 <- concordia.intersection.ab(intercept,slope,wetherill=FALSE,d=x$d)
    } else if (is.na(anchor[[2]])){
        if (x$format < 4){
            b0a0 <- settings('iratio','Pb207Pb206')[1]
        } else {
            a0 <- settings('iratio','Pb206Pb204')[1]
            b0 <- settings('iratio','Pb207Pb204')[1]
            b0a0 <- b0/a0
        }
        intercept <- b0a0
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
    if (x$format>3){
        U238Pb206 <- subset(get.U238Pb206.ratios(x),select='U238Pb206')
        Pb206U238 <- subset(get.Pb206U238.ratios(x),select='Pb206U238')
        Pb204U238 <- subset(get.Pb204U238.ratios(x),select='Pb204U238')
        Pb207Pb206 <- subset(get.Pb207Pb206.ratios(x),select='Pb207Pb206')
        Pb204Pb206 <- Pb204U238/Pb206U238
        if (!anchor[[1]]){
            lmfit <- stats::lm(Pb204Pb206 ~ U238Pb206)
            ta0b0[3] <- ta0b0[2]/lmfit$coef[1]    # 7/4
            ta0b0[2] <- 1/lmfit$coef[1]           # 6/4
        } else if (is.na(anchor[[2]])){
            ta0b0[2] <- a0
            ta0b0[3] <- b0
        } else if (is.numeric(anchor[[2]])){
            lmfit <- stats::lm(Pb204Pb206 ~
                               0 + I(U238Pb206-TW$x['U238Pb206']))
            slope <- lmfit$coef                   # 4/8
            Pb46 <- -TW$x['U238Pb206']*slope
            ta0b0[2] <- 1/Pb46                    # 6/4
            ta0b0[3] <- Pb46 * TW$x['Pb207Pb206'] # 7/4
        }
    }    
    ta0b0
}
get.ta0b0.model3 <- function(x,init,exterr=FALSE,anchor=list(FALSE,NA)){
    out <- list(w=0,par=init)
#    for (i in 1:5){ # loop for more accurate but slower and more unstable results
        out <- fit_ludwig_discordia(x,init=out$par,w=out$w,exterr=exterr,anchor=anchor)
        out$w <- stats::optimize(LL.lud.disp,interval=c(0,1),x=x,ta0b0=out$par,
                                 exterr=exterr,anchor=anchor,maximum=TRUE)$maximum
#    }
    out
}
fit_ludwig_discordia <- function(x,init,w=0,exterr=FALSE,anchor=list(FALSE,NA),...){
    optifix(parms=init,fn=LL.lud.UPb,method="BFGS",x=x,w=w,
            exterr=exterr,fixed=fixit(x,anchor))
}

LL.lud.disp <- function(w,x,ta0b0,exterr=FALSE,anchor=list(FALSE,NA)){
    # these two lines produce slightly more accurate but much slower results:
    # fit <- fit_ludwig_discordia(x,init=ta0b0,w=w,exterr=exterr,anchor=anchor)
    # LL.lud.UPb(ta0b0=fit$par,x=x,exterr=exterr,w=w,LL=TRUE)
    LL.lud.UPb(ta0b0=ta0b0,x=x,exterr=exterr,w=w,LL=TRUE)
}
LL.lud.UPb <- function(ta0b0,x,exterr=FALSE,w=0,LL=FALSE){
    if (x$format<4){
        return(LL.lud.2D(ta0b0,x=x,exterr=exterr,w=w,LL=LL))
    } else {
        return(LL.lud.3D(ta0b0,x=x,exterr=exterr,w=w,LL=LL))
    }
}
LL.lud.2D <- function(ta0,x,exterr=FALSE,w=0,LL=FALSE){
    l <- data2ludwig(x,tt=ta0[1],a0=ta0[2],exterr=exterr,w=w)
    R <- rbind(l$rx,l$ry)
    out <- t(R) %*% l$omega %*% R
    if (LL){
        k <- length(R)
        out <- -0.5*(out + k*log(2*pi) +
                     determinant(2*pi*l$omegainv,logarithm=TRUE)$modulus)
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
        k <- 2*ns
        out <- -0.5*(out + k*log(2*pi) +
                     determinant(l$omegainv,logarithm=TRUE)$modulus)
    }
    out
}

fisher.lud <- function(x,...){ UseMethod("fisher.lud",x) }
fisher.lud.default <- function(x,...){
    stop( "No default method available (yet)." )
}
fisher.lud.UPb <- function(x,fit,exterr=TRUE,anchor=list(FALSE,NA),...){
    ns <- length(x)
    if (x$format<4){
        fish <- fisher_lud_2D(x,fit)
        AA <- fish[1:ns,1:ns]
        BB <- fish[1:ns,(ns+1):(ns+2)]
        CC <- fish[(ns+1):(ns+2),1:ns]
        DD <- fish[(ns+1):(ns+2),(ns+1):(ns+2)]
    } else {
        fish <- fisher_lud_3D(x,fit)
        AA <- fish[1:ns,1:ns]
        BB <- fish[1:ns,(ns+1):(ns+3)]
        CC <- fish[(ns+1):(ns+3),1:ns]
        DD <- fish[(ns+1):(ns+3),(ns+1):(ns+3)]
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
    D <- wendt(tt,d=x$d)
    BB <- (exp(l5*tt)-1+D$d1)/U - a0*(exp(l8*tt)-1+D$d2)
    CC <- -t(l$X)%*%l$omega12 + t(l$x)%*%l$omega12 + t(l$x)%*%l$omega21 -
        t(l$Y)%*%l$omega22 + t(ones)%*%l$omega22*a0 + 2*t(l$x)%*%l$omega22*BB
    DD <- -t(l$X)%*%l$omega12%*%l$x + t(l$x)%*%l$omega12%*%l$x -
        t(l$Y)%*%l$omega22%*%l$x + t(ones)%*%l$omega22%*%l$x*a0 -
        t(l$x)%*%l$omega21%*%l$X + t(l$x)%*%l$omega21%*%l$x -
        t(l$x)%*%l$omega22%*%l$Y + t(l$x)%*%l$omega22%*%ones*a0 +
        2*t(l$x)%*%l$omega22%*%l$x*BB
    dBdt <- (l5*exp(l5*tt)+D$dd1dt)/U - a0*(l8*exp(l8*tt)+D$dd2dt)
    dBda0 <- -(exp(l8*tt)-1+D$d2)
    d2Bda0dt <- -(l8*exp(l8*tt)+D$dd2dt)
    d2Bdtda0 <- d2Bda0dt
    d2Bdt2 <- (l5*l5*exp(l5*tt)+D$d2d1dt2)/U - a0*(l8*l8*exp(l8*tt)+D$d2d2dt2)
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

data2ludwig <- function(x,...){ UseMethod("data2ludwig",x) }
data2ludwig.default <- function(x,...){ stop('default function undefined') }
data2ludwig.UPb <- function(x,tt,a0,b0=0,exterr=FALSE,w=0,...){
    if (x$format<4)
        out <- data2ludwig_2D(x,tt=tt,a0=a0,w=w,exterr=exterr)
    else
        out <- data2ludwig_3D(x,tt=tt,a0=a0,b0=b0,w=w,exterr=exterr)
    out
}
data2ludwig_2D <- function(x,tt,a0,w=0,exterr=FALSE){
    l5 <- settings('lambda','U235')
    l8 <- settings('lambda','U238')
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    D <- wendt(tt=tt,d=x$d)
    f1 <- exp(l5[1]*tt)-1 + D$d1
    f2 <- exp(l8[1]*tt)-1 + D$d2
    B <-  f1/U - a0*f2
    dBdl5 <- (tt*exp(l5[1]*tt)+D$dd1dl5)/U
    dBdl8 <- -a0*(tt*exp(l8[1]*tt)+D$dd2dl8)
    E <- matrix(0,2*ns+2,2*ns+2)
    J <- matrix(0,2*ns,2*ns+2)
    X <- matrix(0,ns,1)
    Y <- matrix(0,ns,1)
    for (i in 1:ns){
        XY <- tera.wasserburg(x,i)
        X[i] <- XY$x['U238Pb206']
        Y[i] <- XY$x['Pb207Pb206']
        E[(2*i-1):(2*i),(2*i-1):(2*i)] <- XY$cov
        E[2*i,2*i] <- E[2*i,2*i] + (a0*w)^2
        J[i,2*i-1] <- 1 # drx/dX
        J[ns+i,2*i] <- 1 # dry/dY
        if (exterr){
            J[ns+i,2*ns+1] <- -dBdl5*X[i]
            J[ns+i,2*ns+2] <- -dBdl8*X[i]
        }
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
    D <- wendt(tt=tt,d=x$d)
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        Z[i] <- wd$x['Pb204U238']
        R[i] <- wd$x['Pb207U235'] - exp(l5[1]*tt) + 1 - U*b0*Z[i] - D$d1
        r[i] <- wd$x['Pb206U238'] - exp(l8[1]*tt) + 1 - a0*Z[i] - D$d2
        Ew <- get.Ew(w=w,Z=Z[i],a0=a0,b0=b0,U=U)
        E[(3*i-2):(3*i),(3*i-2):(3*i)] <- wd$cov + Ew
        J[i,3*i-2] <- 1                                  # dRdX
        J[ns+i,3*i-1] <- 1                               # drdY
        J[2*ns+i,3*i] <- 1                               # dphidZ
        if (exterr){
            J[i,3*ns+1] <- -tt*exp(l5[1]*tt) - D$dd1dl5    # dRdl5
            J[ns+i,3*ns+2] <- -tt*exp(l8[1]*tt) - D$dd2dl8 # drdl8
        }
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
get.Ew <- function(w,Z,a0,b0,U){
    E <- diag(c(a0,b0)*w)^2
    J <- matrix(0,3,2)
    J[1,2] <- -U*Z # dRda0
    J[2,1] <- -Z   # dRdb0
    J %*% E %*% t(J)
}

fixit <- function(x,anchor=list(FALSE,NA)){
    if (x$format<4){
        if (!anchor[[1]]) out <- rep(FALSE,2)
        else if (is.numeric(anchor[[2]])) out <- c(TRUE,FALSE)
        else out <- c(FALSE,TRUE)
    } else {
        if (!anchor[[1]]) out <- rep(FALSE,3)
        else if (is.numeric(anchor[[2]])) out <- c(TRUE,FALSE,FALSE)
        else out <- c(FALSE,TRUE,TRUE)
    }
    out
}
