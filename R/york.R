#' Linear regression of X,Y-variables with correlated errors
#'
#' Implements the unified regression algorithm of York et al. (2004)
#' which, although based on least squares, yields results that are
#' consistent with maximum likelihood estimates of Titterington and
#' Halliday (1979)
#'
#' @param x a 5-column matrix with the X-values, the analytical
#'     uncertainties of the X-values, the Y-values, the analytical
#'     uncertainties of the Y-values, and the correlation coefficients
#'     of the X- and Y-values.
#' @return a four-element list of vectors containing:
#'     \describe{
#'     \item{a}{the intercept of the straight line fit and its
#'     standard error}
#'     \item{b}{the slope of the fit and its standard
#'     error}
#'     \item{cov.ab}{the covariance of the slope and intercept}
#'     \item{mswd}{the mean square of the residuals (a.k.a
#'     `reduced Chi-square') statistic}
#'     }
#' @references
#' Titterington, D.M. and Halliday, A.N., 1979. On the fitting of
#' parallel isochrons and the method of maximum likelihood. Chemical
#' Geology, 26(3), pp.183-195.
#' 
#' York, Derek, et al. "Unified equations for the slope,
#' intercept, and standard errors of the best straight line."
#' American Journal of Physics 72.3 (2004): 367-375.
#'
#' @examples
#'    X <- c(1.550,12.395,20.445,20.435,20.610,24.900,
#'           28.530,50.540,51.595,86.51,106.40,157.35)
#'    Y <- c(.7268,.7849,.8200,.8156,.8160,.8322,
#'           .8642,.9584,.9617,1.135,1.230,1.490)
#'    n <- length(X)
#'    sX <- X*0.01
#'    sY <- Y*0.005
#'    rXY <- rep(0.8,n)
#'    dat <- cbind(X,sX,Y,sY,rXY)
#'    fit <- york(dat)
#'    covmat <- matrix(0,2,2)
#'    plot(range(X),fit$a[1]+fit$b[1]*range(X),type='l',ylim=range(Y))
#'    for (i in 1:n){
#'        covmat[1,1] <- sX[i]^2
#'        covmat[2,2] <- sY[i]^2
#'        covmat[1,2] <- rXY[i]*sX[i]*sY[i]
#'        covmat[2,1] <- covmat[1,2]
#'        ell <- ellipse(X[i],Y[i],covmat,alpha=0.05)
#'        polygon(ell)
#'    }
#' @export
york <- function(x){
    colnames(x) <- c('X','sX','Y','sY','rXY')
    ab <- stats::lm(x[,'Y'] ~ x[,'X'])$coefficients # initial guess
    a <- ab[1]
    b <- ab[2]
    wX <- 1/x[,'sX']^2
    wY <- 1/x[,'sY']^2
    for (i in 1:50){ # 50 = maximum number of iterations
        bold <- b
        alpha <- sqrt(wX*wY)
        W <- wX*wY/(wX+b*b*wY-2*b*x[,'rXY']*alpha)
        Xbar <- sum(W*x[,'X'],na.rm=TRUE)/sum(W,na.rm=TRUE)
        Ybar <- sum(W*x[,'Y'],na.rm=TRUE)/sum(W,na.rm=TRUE)
        U <- x[,'X']-Xbar
        V <- x[,'Y']-Ybar
        beta <- W*(U/wY+b*V/wX-(b*U+V)*x[,'rXY']/alpha)
        b <- sum(W*beta*V,na.rm=TRUE)/sum(W*beta*U,na.rm=TRUE)
        if ((bold/b-1)^2 < 1e-15) break # convergence reached
    }
    a <- Ybar-b*Xbar
    X <- Xbar + beta
    xbar <- sum(W*X,na.rm=TRUE)/sum(W,na.rm=TRUE)
    u <- X-xbar
    sb <- sqrt(1/sum(W*u^2,na.rm=TRUE))
    sa <- sqrt(1/sum(W,na.rm=TRUE)+(xbar*sb)^2)
    out <- list()
    out$a <- c(a,sa)
    out$b <- c(b,sb)
    out$cov.ab <- -Xbar*sb^2
    mswd <- get.york.mswd(x,a,b)
    out <- c(out,mswd)
    out
}

get.york.mswd <- function(x,a,b){
    xy <- get.york.xy(x,a,b)
    X2 <- 0
    nn <- length(x[,'X'])
    for (i in 1:nn){
        E <- cor2cov2(x[i,'sX'],x[i,'sY'],x[i,'rXY'])
        X <- matrix(c(x[i,'X']-xy[i,1],x[i,'Y']-xy[i,2]),1,2)
        if (!any(is.na(X)))
            X2 <- X2 + 0.5*X %*% solve(E) %*% t(X)
    }
    df <- (2*nn-2)
    out <- list()
    out$mswd <- as.numeric(X2/df)
    out$p.value <- as.numeric(1-stats::pchisq(X2,df))
    out
}

# get fitted X and X given a dataset x=cbind(X,sX,Y,sY,rXY),
# an intercept a and slope b. This function is useful
# for evaluating log-likelihoods of derived quantities
get.york.xy <- function(x,a,b){
    wX <- 1/x[,'sX']^2
    wY <- 1/x[,'sY']^2
    alpha <- sqrt(wX*wY)
    W <- wX*wY/(wX+b*b*wY-2*b*x[,'rXY']*alpha)
    Xbar <- sum(W*x[,'X'],na.rm=TRUE)/sum(W,na.rm=TRUE)
    Ybar <- a + b*Xbar
    U <- x[,'X']-Xbar
    V <- x[,'Y']-Ybar
    beta <- W*(U/wY+b*V/wX-(b*U+V)*x[,'rXY']/alpha)
    out <- cbind(Xbar+beta,Ybar+b*beta)
    out
}

data2york <- function(x,...){ UseMethod("data2york",x) }
data2york.default <- function(x,...){ stop('default function undefined') }
data2york.UPb <- function(x,wetherill=TRUE,...){
    ns <- length(x)
    out <- matrix(0,ns,5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    if (wetherill){
        if (x$format %in% c(1,3,4)){
            out <- x$x[,c('Pb207U235','errPb207U235',
                          'Pb206U238','errPb206U238','rhoXY')]
        } else if (x$format %in% c(2,5,6)){
            for (i in 1:ns){
                samp <- wetherill(x,i=i,exterr=FALSE)
                out[i,1] <- samp$x['Pb207U235']
                out[i,2] <- sqrt(samp$cov['Pb207U235','Pb207U235'])
                out[i,3] <- samp$x['Pb206U238']
                out[i,4] <- sqrt(samp$cov['Pb206U238','Pb206U238'])
                out[i,5] <- stats::cov2cor(samp$cov)[1,2]
            }
        }
    } else {
        if (x$format %in% c(2,5)){
            out <- x$x[,c('U238Pb206','errU238Pb206',
                          'Pb207Pb206','errPb207Pb206','rhoXY')]
        } else if (x$format %in% c(1,3,4,6)){
            for (i in 1:ns){
                samp <- tera.wasserburg(x,i=i,exterr=FALSE)
                out[i,1] <- samp$x['U238Pb206']
                out[i,2] <- sqrt(samp$cov['U238Pb206','U238Pb206'])
                out[i,3] <- samp$x['Pb207Pb206']
                out[i,4] <- sqrt(samp$cov['Pb207Pb206','Pb207Pb206'])
                out[i,5] <- stats::cov2cor(samp$cov)[1,2]
            }
        }
    }
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out
}
data2york.ArAr <- function(x,inverse=TRUE,...){
    if (inverse)
        out <- ArAr.inverse.ratios(x)
    else
        out <- ArAr.normal.ratios(x)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out
}
data2york.PbPb <- function(x,inverse=TRUE,...){
    if (inverse)
        out <- PbPb.inverse.ratios(x)
    else
        out <- PbPb.normal.ratios(x)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out
}
data2york.PD <- function(x,exterr=FALSE,common=FALSE,...){
    if (x$format==1){
        out <- x$x
    } else if (x$format==2){
        out <- ppm2ratios(x,exterr=exterr,common=common)
    }
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out
}
