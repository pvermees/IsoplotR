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
#' @param alpha cutoff value for confidence intervals
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
ludwig.UPb <- function(x,exterr=FALSE,alpha=0.05,model=1,
                       anchor=list(FALSE,NA),...){
    fit <- get.ta0b0(x,exterr=exterr,model=model,anchor=anchor)
    out <- fit[c('par','w','model')]
    out$cov <- fisher.lud(x,fit=fit,anchor=anchor)
    mswd <- mswd.lud(fit$par,x=x,anchor=anchor)
    out <- c(out,mswd)
    if (model==3){
        out$w <- c(fit$w,
                   profile_LL_discordia_disp(fit,x=x,alpha=alpha))
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

get.ta0b0 <- function(x,exterr=FALSE,model=1,
                      anchor=list(FALSE,NA)){
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
    out <- fit_ludwig_discordia(x,init=init,w=0,exterr=exterr,
                                fixed=fixit(x,anchor))
    out$w <- 0
    out
}
get.ta0b0.model2 <- function(x,anchor=list(FALSE,NA)){
    xy <- data2york(x,wetherill=FALSE)[,c('X','Y')]
    if (!anchor[[1]]) {
        xyfit <- stats::lm(xy[,'Y'] ~ xy[,'X'])
        intercept <- xyfit$coef[1]
        slope <- xyfit$coef[2]
        ta0b0 <- concordia.intersection.ab(intercept,slope,wetherill=FALSE)$x
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
        ta0b0 <- concordia.intersection.ab(intercept,slope,wetherill=FALSE)$x
    } else if (is.numeric(anchor[[2]])){
        TW <- age_to_terawasserburg_ratios(anchor[[2]],st=0,exterr=FALSE)
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
    fit <- fit_ludwig_discordia(x,init=init,w=0,exterr=exterr,
                                fixed=fixit(x,anchor))
    ta0b0 <- fit$par
    w <- get_ludwig_disp(ta0b0,x,interval=get_lud_wrange(ta0b0,x))
    out <- fit_ludwig_discordia(x,init=ta0b0,w=w,exterr=exterr,
                                fixed=fixit(x,anchor))
    out$w <- w
    out
}
fit_ludwig_discordia <- function(x,init,w=0,exterr=FALSE,...){
    optifix(parms=init,fn=LL.lud.UPb,method="BFGS",x=x,
            w=w,exterr=exterr,...)
}
get_ludwig_disp <- function(ta0b0,x,interval){
    stats::optimize(LL.lud.UPb.disp,interval=interval,
                    x=x,ta0b0=ta0b0,maximum=TRUE)$maximum
}
get_lud_wrange <- function(ta0b0,x){
    c(0,1)
}    

LL.lud.UPb.disp <- function(w,x,ta0b0){
    LL.lud.UPb(ta0b0,x=x,exterr=FALSE,w=w,LL=TRUE)
}
LL.lud.UPb <- function(ta0b0,x,exterr=FALSE,w=0,LL=FALSE){
    if (x$format<4){
        return(LL.lud.2D(ta0b0,x=x,exterr=exterr,w=w,LL=LL))
    } else {
        return(LL.lud.3D(ta0b0,x=x,exterr=exterr,w=w,LL=LL))
    }
}
LL.lud.2D <- function(ta0,x,exterr=FALSE,w=0,LL=FALSE){
    d <- data2ludwig(x,tt=ta0[1],a0=ta0[2],exterr=exterr,w=w)
    R <- rbind(d$rx,d$ry)
    out <- t(R) %*% d$omega %*% R
    if (LL){
        k <- length(R)
        out <- -0.5*(out + k*log(2*pi) +
                     determinant(2*pi*d$omegainv,logarithm=TRUE)$modulus)
    }
    out
}
LL.lud.3D <- function(ta0b0,x,exterr=FALSE,w=0,LL=FALSE){
    d <- data2ludwig(x,tt=ta0b0[1],a0=ta0b0[2],
                     b0=ta0b0[3],exterr=exterr,w=w)
    phi <- d$phi
    R <- d$R
    r <- d$r
    out <- 0
    ns <- length(d$R)
    omega <- d$omega
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
                     determinant(d$omegainv,logarithm=TRUE)$modulus)
    }
    out
}

fisher.lud <- function(x,...){ UseMethod("fisher.lud",x) }
fisher.lud.default <- function(x,...){
    stop( "No default method available (yet)." )
}
fisher.lud.UPb <- function(x,fit,exterr=TRUE,
                           anchor=list(FALSE,NA),...){
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
    d <- data2ludwig_2D(x,tt=tt,a0=a0,w=0,exterr=fit$exterr)
    l5 <- settings('lambda','U235')[1]
    l8 <- settings('lambda','U238')[1]
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    ones <- matrix(1,ns,1)
    B <-  (exp(l5*tt)-1)/U - a0*(exp(l8*tt)-1)
    C <- -t(d$X)%*%d$omega12 + t(d$x)%*%d$omega12 + t(d$x)%*%d$omega21 -
        t(d$Y)%*%d$omega22 + t(ones)%*%d$omega22*a0 + 2*t(d$x)%*%d$omega22*B
    D <- -t(d$X)%*%d$omega12%*%d$x + t(d$x)%*%d$omega12%*%d$x -
        t(d$Y)%*%d$omega22%*%d$x + t(ones)%*%d$omega22%*%d$x*a0 -
        t(d$x)%*%d$omega21%*%d$X + t(d$x)%*%d$omega21%*%d$x -
        t(d$x)%*%d$omega22%*%d$Y + t(d$x)%*%d$omega22%*%ones*a0 +
        2*t(d$x)%*%d$omega22%*%d$x*B
    dBdt <- l5*exp(l5*tt)/U - a0*exp(l8*tt)*l8
    dBda0 <- -(exp(l8*tt)-1)
    d2Bdtda0 <- exp(l8*tt)*l8
    d2Bda0dt <- exp(l8*tt)*l8
    d2Bdt2 <- (l5^2)*exp(l5*tt)/U - a0*exp(l8*tt)*(l8^2)
    d2Bda02 <- 0
    dDda0 <- t(ones)%*%d$omega22%*%d$x +
        t(d$x)%*%d$omega22%*%ones + 2*t(d$x)%*%d$omega22%*%d$x*dBda0
    ns <- length(x)
    out <- matrix(0,ns+2,ns+2)
    out[1:ns,1:ns] <- # d2S/dx2
        2*(d$omega11 + d$omega21*B + d$omega12*B + d$omega22*B^2)

    out[ns+2,1:ns] <-   # d2S/dxda0
        2*(C*dBda0 + t(ones)%*%d$omega21 + t(ones)%*%d$omega22*B)
    out[ns+2,ns+2] <- # d2S/da02
        2*(t(ones)%*%d$omega22%*%ones +
           t(ones)%*%d$omega22%*%d$x*dBda0 +
           t(d$x)%*%d$omega22%*%ones*dBda0 +
           t(d$x)%*%d$omega22%*%d$x*(dBda0^2))
    out[ns+1,1:ns] <- 2*C*dBdt # d2S/dxdt
    out[ns+1,ns+1] <- # d2S/dt2
        D*d2Bdt2 + 2*t(d$x)%*%d$omega22%*%d$x*(dBdt^2)
    out[ns+1,ns+2] <- # d2S/dt0da0
        D*d2Bdtda0 + dDda0*dBdt
    out[1:ns,ns+1] <- t(out[ns+1,1:ns]) # d2S/dtdx
    out[1:ns,ns+2] <- t(out[ns+2,1:ns]) # d2S/da0dx
    out[ns+2,ns+1] <- out[ns+1,ns+2]
    out/2
}
fisher_lud_3D <- function(x,fit){
    tt <- fit$par[1]
    a0 <- fit$par[2]
    b0 <- fit$par[3]
    d <- data2ludwig_3D(x,tt=tt,a0=a0,b0=b0,w=0,exterr=fit$exterr)
    z <- d$z
    omega <- d$omega
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
    B <-  (exp(l5[1]*tt)-1)/U - a0*(exp(l8[1]*tt)-1)
    dBdl5 <- tt*exp(l5[1]*tt)/U
    dBdl8 <- -a0*tt*exp(l8[1]*tt)
    E <- matrix(0,2*ns+2,2*ns+2)
    J <- matrix(0,2*ns,2*ns+2)
    X <- matrix(0,ns,1)
    Y <- matrix(0,ns,1)
    for (i in 1:ns){
        XY <- tera.wasserburg(x,i,exterr=FALSE)
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
    ER <- J %*% E %*% t(J)
    out <- list()
    omega <- solve(ER)
    omega11 <- omega[1:ns,1:ns]
    omega12 <- omega[1:ns,(ns+1):(2*ns)]
    omega21 <- omega[(ns+1):(2*ns),1:ns]
    omega22 <- omega[(ns+1):(2*ns),(ns+1):(2*ns)]
    out$x <- solve(omega11 + B*(omega12+omega21) + omega22*B^2) %*%
        (omega11%*%X + B*omega12%*%X + omega21%*%(Y-a0) + B*omega22%*%(Y-a0))
    out$omegainv <- ER
    out$omega <- omega
    out$omega11 <- omega11
    out$omega12 <- omega12
    out$omega21 <- omega21
    out$omega22 <- omega22
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
    R <- rep(0,ns)
    r <- rep(0,ns)
    Z <- rep(0,ns)
    E <- matrix(0,3*ns,3*ns)
    if (exterr){
        P235 <- tt*exp(l5[1]*tt)
        P238 <- tt*exp(l8[1]*tt)
        E[1:ns,1:ns] <- (P235*l5[2])^2 # A
        E[(ns+1):(2*ns),(ns+1):(2*ns)] <- (P238*l8[2])^2 # B
    }
    for (i in 1:ns){
        d <- wetherill(x,i=i,exterr=FALSE)
        Z[i] <- d$x['Pb204U238']
        R[i] <- d$x['Pb207U235'] - exp(l5[1]*tt) + 1 - U*b0*Z[i]
        r[i] <- d$x['Pb206U238'] - exp(l8[1]*tt) + 1 - a0*Z[i]
        Ew <- get.Ew(w=w,Z=Z[i],a0=a0,b0=b0,U=U)
        E[i,i] <- E[i,i] + d$cov['Pb207U235','Pb207U235'] + Ew[1,1] # A
        E[ns+i,ns+i] <- E[ns+i,ns+i] + d$cov['Pb206U238','Pb206U238'] + Ew[2,2] # B
        E[2*ns+i,2*ns+i] <- d$cov['Pb204U238','Pb204U238'] # C
        E[i,ns+i] <- d$cov['Pb207U235','Pb206U238'] + Ew[1,2] # D
        E[ns+i,i] <- E[i,ns+i] + Ew[2,1]
        E[i,2*ns+i] <- d$cov['Pb207U235','Pb204U238'] # E
        E[2*ns+i,i] <- E[i,2*ns+i]
        E[ns+i,2*ns+i] <- d$cov['Pb206U238','Pb204U238'] # F
        E[2*ns+i,ns+i] <- E[ns+i,2*ns+i]
    }
    omega <- solve(E)
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
    J <- matrix(0,2,4)
    J[1,1] <- Z*U
    J[1,3] <- a0*U
    J[1,4] <- a0*Z
    J[2,2] <- Z
    J[2,3] <- b0
    E <- matrix(0,4,4)
    E[1,1] <- (a0*w)^2
    E[2,2] <- (b0*w)^2
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
