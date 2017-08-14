#' Linear regression of X,Y,Z-variables with correlated errors, taking
#' into account decay constant uncertainties.
#'
#' Implements the maximum likelihood algorithm of Ludwig (1998)
#'
#' @param x an object of class \code{UPb} with \code{x$format > 3}.
#'
#' @param ... optional arguments
#'
# @param x a \eqn{3n}-element vector \eqn{[X Y Z]}, where \eqn{X},
#     \eqn{Y} and \eqn{Z} are three \eqn{n}-element vectors of
#     (isotopic ratio) values.
# @param covmat a \eqn{[3n x 3n]}-element covariance matrix of
#     \code{x}
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
ludwig.UPb <- function(x,exterr=FALSE,...){
    if (exterr){
        init <- ludwig(x,exterr=FALSE)$par
    } else {
        t0 <- concordia.intersection.york(x,wetherill=TRUE,exterr=FALSE)$x[1]
        init <- c(t0,10,10)
    }
    fit <- stats::optim(init,fn=LL.lud.UPb,method="BFGS",x=x,exterr=exterr)
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
    parnames <- c('t','64i','74i')
    names(out$par) <- parnames
    rownames(out$cov) <- parnames
    colnames(out$cov) <- parnames
    mswd <- mswd.lud(fit$par,x=x)
    out <- c(out,mswd)
    out
}

mswd.lud <- function(ta0b0,x){
    tt <- ta0b0[1]
    a0 <- ta0b0[2]
    b0 <- ta0b0[3]
    ns <- length(x)
    SS <- 2*LL.lud.UPb(ta0b0,x=x,exterr=FALSE)
    df <- 2*ns-2
    out <- list()
    out$mswd <- SS/df
    out$p.value <- as.numeric(1-stats::pchisq(SS,df))
    out
}

fisher.lud <- function(x,...){ UseMethod("fisher.lud",x) }
fisher.lud.default <- function(x,...){ stop( "No default method available (yet)." ) }
fisher.lud.UPb <- function(x,ta0b0,exterr=FALSE,...){
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

LL.lud.UPb <- function(ta0b0,x,exterr=FALSE){
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
