#' Linear regression of X,Y,Z-variables with correlated errors, taking
#' into account decay constant uncertainties.
#'
#' Implements the maximum likelihood algorithm of Ludwig (1998)
#'
#' @param x an object of class \code{UPb}
#' @param alpha cutoff value for confidence intervals
#' @param ... optional arguments
#'
# @param x a \eqn{3n}-element vector \eqn{[X Y Z]}, where \eqn{X},
#     \eqn{Y} and \eqn{Z} are three \eqn{n}-element vectors of
#     (isotopic ratio) values.
# @param covmat a \eqn{[3n x 3n]}-element covariance matrix of
#     \code{x}
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
#' }
#'
#' @examples
#' f <- system.file("UPb4.csv",package="IsoplotR")
#' d <- read.data(f,method="U-Pb",format=4)
#' fit <- ludwig(d)
#' @references
#' Ludwig, K.R., 1998. On the treatment of concordant uranium-lead
#' ages. Geochimica et Cosmochimica Acta, 62(4), pp.665-676.
#' @export
ludwig <- function(x,...){ UseMethod("ludwig",x) }
#' @rdname ludwig
#' @export
ludwig.default <- function(x,...){ stop( "No default method available (yet)." ) }
#' @param exterr propagate external sources of uncertainty (e.g., decay constant)?
#' @rdname ludwig
#' @export
ludwig.UPb <- function(x,exterr=FALSE,alpha=0.05,...){
    ta0 <- concordia.intersection.york(x,exterr=FALSE)$x
    if (x$format<4) init <- ta0
    else init <- c(ta0[1],10,10)
    fit <- stats::optim(init,fn=LL.lud.UPb,
                        method="BFGS",x=x,exterr=exterr)
    out <- list()
    out$par <- fit$par
    out$cov <- tryCatch({ # analytical
        fish <- fisher.lud(x,ta0b0=fit$par,exterr=exterr)
        ns <- length(x)
        # block inversion:
        AA <- fish[1:ns,1:ns]
        BB <- fish[1:ns,(ns+1):(ns+3)]
        CC <- fish[(ns+1):(ns+3),1:ns]
        DD <- fish[(ns+1):(ns+3),(ns+1):(ns+3)]
        solve(DD - CC %*% solve(AA) %*% BB)
    }, error = function(e){ # numerical
        fit <- stats::optim(init,fn=LL.lud.UPb,method="BFGS",
                            x=x,exterr=exterr,hessian=TRUE)
        solve(fit$hessian)
    })
    if (x$format<4) parnames <- c('t[l]','76i')
    else parnames <- c('t','64i','74i')
    names(out$par) <- parnames
    rownames(out$cov) <- parnames
    colnames(out$cov) <- parnames
    mswd <- mswd.lud(fit$par,x=x)
    out <- c(out,mswd)
    out
}

mswd.lud <- function(pars,x){
    ns <- length(x)
    # Mistake in Ludwig (1998)? Multiply the following by 2?
    SS <- LL.lud.UPb(pars,x=x,exterr=FALSE)
    out <- list()
    if (x$format<4) out$df <- ns-2
    else out$df <- 2*ns-2
    out$mswd <- as.vector(SS/out$df)
    out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    out
}

LL.lud.UPb <- function(pars,x,exterr=FALSE){
    if (x$format<4) return(LL.lud.2D(pars,x=x,exterr=exterr))
    else return(LL.lud.3D(pars,x=x,exterr=exterr))
}
LL.lud.2D <- function(ta0,x,exterr=FALSE){
    tt <- ta0[1]
    a0 <- ta0[2]
    l5 <- settings('lambda','U235')
    l8 <- settings('lambda','U238')
    U <- settings('iratio','U238U235')[1]
    XY <- data2york(x,wetherill=FALSE)
    b <- (exp(l5[1]*tt)-1)/U - a0*(exp(l8[1]*tt)-1) # slope
    xy <- get.york.xy(XY,a0,b) # get adjusted xi, yi
    ns <- length(x)
    
    v <- matrix(0,1,2*ns)
    v[1:ns] <- XY[,'X']-xy[,1]
    v[(ns+1):(2*ns)] <- XY[,'Y']-xy[,2]
    if (exterr){
        Ex <- matrix(0,2*ns+2,2*ns+2)
        Jv <- diag(1,2*ns,2*ns+2)
    } else {
        Ex <- matrix(0,2*ns,2*ns)
        Jv <- diag(1,2*ns,2*ns)
    }
    Ex[1:ns,1:ns] <- diag(XY[,'sX'])^2
    Ex[(ns+1):(2*ns),(ns+1):(2*ns)] <- diag(XY[,'sY'])^2
    Ex[1:ns,(ns+1):(2*ns)] <- diag(XY[,'rXY'])*diag(XY[,'sY'])*diag(XY[,'sY'])
    Ex[(ns+1):(2*ns),1:ns] <- Ex[1:ns,(ns+1):(2*ns)]
    if (exterr){
        Ex[2*ns+1,2*ns+1] <- l5[2]^2
        Ex[2*ns+2,2*ns+2] <- l8[2]^2
        Jv[,2*ns+1] <- -tt*exp(l5[1]*tt)*XY[,'X']/U
        Jv[,2*ns+2] <- a0*tt*exp(l8[1]*tt)*XY[,'X']
    }
    E <- Jv %*% Ex %*% t(Jv)
    S <- v %*% solve(E) %*% t(v)
    S/2
}
LL.lud.3D <- function(ta0b0,x,exterr=FALSE){
    tt <- ta0b0[1]
    a0 <- ta0b0[2]
    b0 <- ta0b0[3]
    d <- data2ludwig(x,a0=a0,b0=b0,tt=tt,exterr=exterr)
    phi <- d$phi
    R <- d$R
    r <- d$r
    SS <- 0
    ns <- length(d$R)
    if (exterr){
        omega <- d$omega
        for (i in 1:ns){
            i1 <- i
            i2 <- i + ns
            i3 <- i + 2*ns
            for (j in 1:ns){
                j1 <- j
                j2 <- j + ns
                j3 <- j + 2*ns
                SS <- SS + R[i]*R[j]*omega[i1,j1] + r[i]*r[j]*omega[i2,j2] +
                    phi[i]*phi[j]*omega[i3,j3] + 2*( R[i]*r[j]*omega[i1,j2] +
                    R[i]*phi[j]*omega[i1,j3] + r[i]*phi[j]*omega[i2,j3] )
            }
        }
    } else {
        for (i in 1:ns){
            omega <- d$omega[[i]]
            SS <- SS + omega[1,1]*R[i]^2 + omega[2,2]*r[i]^2 + omega[3,3]*phi[i]^2 +
                2*( R[i]*r[i]*omega[1,2] + R[i]*phi[i]*omega[1,3] +
                r[i]*phi[i]*omega[2,3] )
        }
    }
    SS/2
}

fisher.lud <- function(x,...){ UseMethod("fisher.lud",x) }
fisher.lud.default <- function(x,...){ stop( "No default method available (yet)." ) }
fisher.lud.UPb <- function(x,pars,exterr=FALSE,...){
    if (x$format<4) return(fisher.lud.2D(x,pars,exterr=FALSE,...))
    else return(fisher.lud.3D(x,pars,exterr=FALSE,...))
}
fisher.lud.2D <- function(x,ta0,exterr=FALSE,...){
    stop('not implemented yet')
}
fisher.lud.3D <- function(x,ta0b0,exterr=FALSE,...){
    tt <- ta0b0[1]
    a0 <- ta0b0[2]
    b0 <- ta0b0[3]
    d <- data2ludwig(x,a0=a0,b0=b0,tt=tt,z=z,exterr=exterr)
    z <- d$z
    omega <- d$omega
    if (exterr)
        out <- fisher_lud_with_decay_err(tt,a0=a0,b0=b0,z=z,omega=omega)
    else
        out <- fisher_lud_without_decay_err(tt,a0=a0,b0=b0,z=z,omega=omega)
    out
}
fisher_lud_without_decay_err <- function(tt,a0,b0,z,omega){
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
        O <- omega[[i]]
        d2L.dt2 <- d2L.dt2 +
            O[1,1]*Q235^2 + O[2,2]*Q238^2 + 2*Q235*Q238*O[1,2]
        d2L.da02 <- d2L.da02 + O[2,2]*z[i]^2
        d2L.db02 <- d2L.db02 + O[1,1]*(U*z[i])^2
        d2L.da0dt <- d2L.da0dt + z[i]*(Q238*O[2,2] + Q235*O[1,2])
        d2L.db0dt <- d2L.db0dt + U*z[i]*(Q235*O[1,1] + Q238*O[1,2])
        d2L.da0db0 <- d2L.da0db0 + U*O[1,2]*z[i]^2
        d2L.dz2 <- O[1,1]*(U*b0)^2 + O[2,2]*a0^2 + O[3,3] +
            2*(a0*U*b0*O[1,2] + U*b0*O[1,3] + a0*O[3,3])
        d2L.dzdt <- Q235*(U*b0*O[1,1] + O[1,3]) +
            Q238*(a0*O[2,2]+O[2,3]) +
            (Q238*U*b0 + Q235*a0)*O[1,2]
        d2L.dzda0 <- z[i]*(a0*O[2,2] + U*b0*O[1,2] + O[2,3])
        d2L.dzdb0 <- U*z[i]*(a0*O[1,2] + U*b0*O[1,1] + O[1,3])
        out[i,i] <- d2L.dz2
        out[i,ns+1] <- d2L.dzdt
        out[i,ns+2] <- d2L.dzda0
        out[i,ns+3] <- d2L.dzdb0
        out[ns+1,i] <- out[i,ns+1]
        out[ns+2,i] <- out[i,ns+2]
        out[ns+3,i] <- out[i,ns+3]
    }
    out[ns+1,ns+1] <- d2L.dt2
    out[ns+2,ns+2] <- d2L.da02
    out[ns+3,ns+3] <- d2L.db02
    out[ns+1,ns+2] <- d2L.da0dt
    out[ns+1,ns+3] <- d2L.db0dt
    out[ns+2,ns+3] <- d2L.da0db0
    out[ns+2,ns+1] <- out[ns+1,ns+2]
    out[ns+3,ns+1] <- out[ns+1,ns+3]
    out[ns+3,ns+2] <- out[ns+2,ns+3]
    out
}
fisher_lud_with_decay_err <- function(tt,a0,b0,z,omega){
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
            d2L.dzdz <-  (omega[i1,j1] + omega[i1,j3] + omega[i3,j1])*(U*b0)^2 +
                omega[i2,j2]*a0^2 + omega[i3,j3] + a0*(omega[i2,j3]+omega[i3,j2]) +
                                    a0*U*b0*(omega[i1,j2]+omega[i2,j1])
            out[i,j] <- d2L.dzdz # dij
            out[j,i] <- d2L.dzdz # dij
        }
        d2L.dz2 <- omega[i1,i1]*(U*b0)^2 + omega[i2,i2]*a0^2 + omega[i3,i3] +
                   2*(a0*U*b0*omega[i1,i2] + U*b0*omega[i1,i3] + a0*omega[i3,i3])
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
data2ludwig.UPb <- function(x,a0,b0,tt,exterr=FALSE,...){
    if (x$format < 4)
        stop('Ludwig regression is not possible for U-Pb data of format < 4.')
    if (exterr)
        out <- data2ludwig_with_decay_err(x,a0=a0,b0=b0,tt=tt)
    else
        out <- data2ludwig_without_decay_err(x,a0=a0,b0=b0,tt=tt)
    out
}
data2ludwig_without_decay_err <- function(x,a0,b0,tt){
    # initialise:
    l5 <- settings('lambda','U235')
    l8 <- settings('lambda','U238')
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    R <- rep(0,ns)
    r <- rep(0,ns)
    phi <- rep(0,ns)
    Z <- rep(0,ns)
    z <- rep(0,ns)
    omega <- list()
    E <- matrix(0,3,3)
    for (i in 1:ns){
        d <- wetherill(x,i=i,exterr=FALSE)
        Z[i] <- d$x['Pb204U238']
        R[i] <- d$x['Pb207U235'] - exp(l5[1]*tt) + 1 - U*b0*Z[i]
        r[i] <- d$x['Pb206U238'] - exp(l8[1]*tt) + 1 - a0*Z[i]
        E[1,1] <- d$cov['Pb207U235','Pb207U235']
        E[2,2] <- d$cov['Pb206U238','Pb206U238']
        E[3,3] <- d$cov['Pb204U238','Pb204U238']
        E[1,2] <- d$cov['Pb207U235','Pb206U238']
        E[1,3] <- d$cov['Pb207U235','Pb204U238']
        E[2,3] <- d$cov['Pb206U238','Pb204U238']
        E[2,1] <- E[1,2]
        E[3,1] <- E[1,3]
        E[3,2] <- E[2,3]
        O <- solve(E)
        omega[[i]] <- O
        # rearrange sum of squares:
        AA <- O[1,1]*(U*b0)^2 + O[2,2]*a0^2 + O[3,3] +
            2*(U*a0*b0*O[1,2] + U*b0*O[1,3] + a0*O[2,3])
        BB <- R[i]*U*b0*O[1,1] + r[i]*a0*O[2,2] +
             (R[i]*a0 + U*b0*r[i])*O[1,2] +
            R[i]*O[1,3] + r[i]*O[2,3]
        CC <- O[1,1]*R[i]^2 + O[2,2]*r[i]^2 + 2*R[i]*r[i]*O[1,2]
        z[i] <- Z[i] + BB/AA
        phi[i] <- Z[i] - z[i]
    }
    out <- list(R=R,r=r,phi=phi,z=z,omega=omega)
}
data2ludwig_with_decay_err <- function(x,a0,b0,tt){
    # initialise:
    l5 <- settings('lambda','U235')
    l8 <- settings('lambda','U238')
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    R <- rep(0,ns)
    r <- rep(0,ns)
    Z <- rep(0,ns)
    E <- matrix(0,3*ns,3*ns)
    if (TRUE){ # TRUE if propagating decay errors
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
        E[i,i] <- E[i,i] + d$cov['Pb207U235','Pb207U235'] # A
        E[ns+i,ns+i] <- E[ns+i,ns+i] + d$cov['Pb206U238','Pb206U238'] # B
        E[2*ns+i,2*ns+i] <- d$cov['Pb204U238','Pb204U238'] # C
        E[i,ns+i] <- d$cov['Pb207U235','Pb206U238'] # D
        E[ns+i,i] <- E[i,ns+i]
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
            W[i] <- W[i] - R[j] * (U*b0*omega[i1,j1] + a0*omega[i2,j1] + omega[i3,j1]) -
                           r[j] * (U*b0*omega[i1,j2] + a0*omega[i2,j2] + omega[i3,j2])
            # Ludwig (1998): V[i,j] <- U*b0*omega[i1,j3] + a0*omega[i2,j3] + omega[i3,j3]
            V[i,j] <- omega[i1,j1]*(U*b0)^2 + omega[i2,j2]*a0^2 + omega[i3,j3] +
                      2*(U*a0*b0*omega[i1,j2] + U*b0*omega[i1,j3] + a0*omega[i2,j3])
        }
    }
    phi <- solve(V,W)
    z <- Z - phi
    out <- list(R=R,r=r,phi=phi,z=z,omega=omega)
}
