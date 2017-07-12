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
ludwig <- function(x){
    ns <- length(x)
    tt <- discordia.age(x,wetherill=TRUE,exterr=FALSE)
    fitXY <- york(x[,c(1,2,3,4,7)])
    a <- fitXY$a[1]
    b <- fitXY$b[1]
    fitXZ <- york(x[,c(1,2,5,6,8)])
    A <- fitXZ$a[1]
    B <- fitXZ$b[1]
    init <- c(a,b,A,B)
    dat <- matrix2covlist(x)
    fit <- stats::optim(init,fn=S.tit,gr=gr.tit,method="BFGS",dat)
    fish <- fisher.tit(fit$par,dat)
    covmat <- solve(fish)
}

get.omega.ludwig <- function(x,exterr=FALSE,tt=NA,...){
    if (x$format < 4)
        stop('Ludwig regression is not possible for U-Pb data of format < 4.')
    ns <- length(x)
    E <- matrix(0,3*ns,3*ns)
    if (exterr){
        l5 <- settings('lambda','U235')
        l8 <- settings('lambda','U238')
        P235 <- tt*exp(l5[1]*tt)
        P238 <- tt*exp(l8[1]*tt)
        E[1:ns,1:ns] <- (P235*l5[2])^2 # A 
        E[(ns+1):(2*ns),(ns+1):(2*ns)] <- (P238*l8[2])^2 # B
    }
    for (i in 1:ns){
        d <- wetherill(x,i=i,exterr=FALSE) # will add external errors based on Ludwig (1999)
        out$R[i] <- d$x['Pb207U235']-exp()
        out$r[i] <- d$x['Pb206U238']
        out$phi[i] <- d$x['Pb204U238']
        E[i,i] <- E[i,i] + d$cov['Pb207U235','Pb207U235'] # A
        E[ns+i,ns+i] <- E[ns+i,ns+i] + d$cov['Pb206U238','Pb206U238'] # B
        E[2*ns+i,2*ns+i] <- d$cov['Pb204U238','Pb204U238'] # C
        E[i,ns+i] <- d$cov['Pb207U235','Pb206U238'] # D
        E[ns+i,i] <- E[i,ns+i]
        E[ns+i,2*ns+i] <- d$cov['Pb207U235','Pb204U238'] # E
        E[2*ns+i,ns+i] <- E[ns+i,2*ns+i]
        E[ns+i,2*ns+i] <- d$cov['Pb206U238','Pb204U238'] # F
        E[2*ns+i,ns+i] <- E[ns+i,2*ns+i]
    }
    solve(E)
}

data2ludwig <- function(x,...){ UseMethod("data2ludwig",x) }
data2ludwig.default <- function(x,...){ stop('default function undefined') }
data2ludwig.UPb <- function(x,exterr=FALSE,tt=NA,...){

}
