# Linear regression of X,Y,Z-variables with correlated errors, taking
# into account decay constant uncertainties
#
# Implements the maximum likelihood algorithm of Ludwig (1998)
#
# @param x a \eqn{3n}-element vector \eqn{[X Y Z]}, where \eqn{X},
#     \eqn{Y} and \eqn{Z} are three \eqn{n}-element vectors of
#     (isotopic ratio) values.
# @param covmat a \eqn{[3n x 3n]}-element covariance matrix of
#     \code{x}
# @references
# Ludwig, K.R., 1998. On the treatment of concordant uranium-lead
# ages. Geochimica et Cosmochimica Acta, 62(4), pp.665-676.
#
ludwig <- function(x,...){ UseMethod("ludwig",x) }
ludwig.default <- function(x,covmat,...){
}
ludwig.UPb <- function(x,exterr=FALSE,...){
    if (exterr){
        init <- ludwig(x,exterr=FALSE)$par
    } else {
        ns <- length(x)
        t0 <- discordia.age(x,wetherill=TRUE,exterr=FALSE)$x[1]
        init <- c(10,10,t0)
    }
    fit <- stats::optim(init,fn=S.lud,method="BFGS",x=x,exterr=exterr)
    fit
}

S.lud <- function(a0b0t,x,exterr=FALSE){
    ns <- length(x)
    a0 <- a0b0t[1]
    b0 <- a0b0t[2]
    tt <- a0b0t[3]
    d <- data2ludwig(x,a0,b0,tt,exterr=exterr)
    phi <- d$phi
    R <- d$R
    r <- d$r
    omega <- d$omega
    SS <- 0
    for (i in 1:ns){
        i1 <- i
        i2 <- i + ns
        i3 <- i + 2*ns
        for (j in 1:ns){
            j1 <- j
            j2 <- j + ns
            j3 <- j + 2*ns
            SS <- SS + R[i]*R[j]*omega[i1,j1] + r[i]*r[j]*omega[i2,j2] *
                  phi[i]*phi[j]*omega[i3,i3] + 2*( R[i]*r[j]*omega[i1,j2] +
                  R[i]*phi[j]*omega[i1,j3] + r[i]*phi[j]*omega[i2,j3] ) 
        }
    }
    SS
}

data2ludwig <- function(x,...){ UseMethod("data2ludwig",x) }
data2ludwig.default <- function(x,...){ stop('default function undefined') }
data2ludwig.UPb <- function(x,a0,b0,tt,exterr=FALSE,...){
    if (x$format < 4)
        stop('Ludwig regression is not possible for U-Pb data of format < 4.')
    # initialise:
    l5 <- settings('lambda','U235')
    l8 <- settings('lambda','U238')
    U <- settings('iratio','U238U235')[1]
    ns <- length(x)
    R <- rep(0,ns)
    r <- rep(0,ns)
    E <- matrix(0,3*ns,3*ns)
    # populate R, r, phi and omega
    if (exterr){
        P235 <- tt*exp(l5[1]*tt)
        P238 <- tt*exp(l8[1]*tt)
        E[1:ns,1:ns] <- (P235*l5[2])^2 # A 
        E[(ns+1):(2*ns),(ns+1):(2*ns)] <- (P238*l8[2])^2 # B
    }
    for (i in 1:ns){
        d <- wetherill(x,i=i,exterr=FALSE) # will add external errors based on Ludwig (1998)
        Zi <- d$x['Pb204U238']
        R[i] <- d$x['Pb207U235'] - exp(l5[1]*tt) + 1 - U*b0*Zi
        r[i] <- d$x['Pb206U238'] - exp(l8[1]*tt) + 1 - a0*Zi
        E[i,i] <- E[i,i] + d$cov['Pb207U235','Pb207U235'] # A
        E[ns+i,ns+i] <- E[ns+i,ns+i] + d$cov['Pb206U238','Pb206U238'] # B
        E[2*ns+i,2*ns+i] <- d$cov['Pb204U238','Pb204U238'] # C
        E[i,ns+i] <- d$cov['Pb207U235','Pb206U238'] # D
        E[ns+i,i] <- E[i,ns+i]
        E[i,2*ns+i] <- d$cov['Pb207U235','Pb204U238'] # E
        E[2*ns+i,i] <- E[ns+i,2*ns+i]
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
            W[i] <- W[i] - R[j] * (U*b0*omega[i1,j1] + a0*omega[i2,j1] + omega[i3,j1]) +
                r[j] * (U*b0*omega[i1,j2] + a0*omega[i2,j2] + omega[i3,j2])
            V[i,j] <- U*b0*omega[i1,j3] + a0*omega[i2,j3] + omega[i3,j3]
        }
    }
    phi <- solve(V,W)
    list(R=R,r=r,omega=omega,phi=phi)
}
