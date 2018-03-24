#' Linear regression of X,Y-variables with correlated errors
#'
#' Implements the unified regression algorithm of York et al. (2004)
#' which, although based on least squares, yields results that are
#' consistent with maximum likelihood estimates of Titterington and
#' Halliday (1979)
#'
#' @details
#' Given n pairs of (approximately) collinear measurements \eqn{X_i}
#' and \eqn{Y_i} (for \eqn{1 \leq i \leq n}), their uncertainties
#' \eqn{s[X_i]} and \eqn{s[Y_i]}, and their covariances
#' cov[\eqn{X_i,Y_i}], the \code{york} function finds the best fitting
#' straight line using the least-squares algorithm of York et
#' al. (2004). This algorithm is modified from an earlier method
#' developed by York (1968) to be consistent with the maximum
#' likelihood approach of Titterington and Halliday (1979). It
#' computes the MSWD as a measure of under/overdispersion.
#' Overdispersed datasets (MSWD>1) can be dealt with in the same three
#' ways that are described in the documentation of the
#' \code{\link{isochron}} function.
#'
#' @param x a 5-column matrix with the X-values, the analytical
#'     uncertainties of the X-values, the Y-values, the analytical
#'     uncertainties of the Y-values, and the correlation coefficients
#'     of the X- and Y-values.
#' @param alpha cutoff value for confidence intervals
#' @return A four-element list of vectors containing:
#'
#'     \describe{
#'
#'     \item{a}{the intercept of the straight line fit and its
#'     standard error}
#'
#'     \item{b}{the slope of the fit and its standard error}
#'
#'     \item{cov.ab}{the covariance of the slope and intercept}
#'
#'     \item{mswd}{the mean square of the residuals (a.k.a
#'     `reduced Chi-square') statistic}
#'
#'     \item{df}{degrees of freedom of the linear fit \eqn{(2n-2)}}
#'
#'     \item{p.value}{p-value of a Chi-square value with \code{df}
#'     degrees of freedom}
#'
#'     }
#' @seealso \code{\link{york}}, \code{\link{isochron}}, \code{\link{ludwig}}
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
york <- function(x,alpha=0.05){
    colnames(x) <- c('X','sX','Y','sY','rXY')
    ab <- stats::lm(x[,'Y'] ~ x[,'X'])$coefficients # initial guess
    a <- ab[1]
    b <- ab[2]
    if (any(is.na(ab)))
        stop('Cannot fit a straight line through these data')
    wX <- 1/x[,'sX']^2
    wY <- 1/x[,'sY']^2
    for (i in 1:50){ # 50 = maximum number of iterations
        bold <- b
        A <- sqrt(wX*wY)
        W <- wX*wY/(wX+b*b*wY-2*b*x[,'rXY']*A)
        Xbar <- sum(W*x[,'X'],na.rm=TRUE)/sum(W,na.rm=TRUE)
        Ybar <- sum(W*x[,'Y'],na.rm=TRUE)/sum(W,na.rm=TRUE)
        U <- x[,'X']-Xbar
        V <- x[,'Y']-Ybar
        B <- W*(U/wY+b*V/wX-(b*U+V)*x[,'rXY']/A)
        b <- sum(W*B*V,na.rm=TRUE)/sum(W*B*U,na.rm=TRUE)
        if ((bold/b-1)^2 < 1e-15) break # convergence reached
    }
    a <- Ybar-b*Xbar
    X <- Xbar + B
    xbar <- sum(W*X,na.rm=TRUE)/sum(W,na.rm=TRUE)
    u <- X-xbar
    sb <- sqrt(1/sum(W*u^2,na.rm=TRUE))
    sa <- sqrt(1/sum(W,na.rm=TRUE)+(xbar*sb)^2)
    out <- get.york.mswd(x,a,b)
    out$a <- c(a,sa)
    out$b <- c(b,sb)
    out$cov.ab <- -Xbar*sb^2
    names(out$a) <- c('a','s[a]')
    names(out$b) <- c('b','s[b]')
    out$type <- 'york'
    out
}

get.york.mswd <- function(x,a,b){
    xy <- get.york.xy(x,a,b)
    X2 <- 0
    ns <- length(x[,'X'])
    for (i in 1:ns){
        E <- cor2cov2(x[i,'sX'],x[i,'sY'],x[i,'rXY'])
        X <- matrix(c(x[i,'X']-xy[i,1],x[i,'Y']-xy[i,2]),1,2)
        if (!any(is.na(X)))
            X2 <- X2 + 0.5*X %*% solve(E) %*% t(X)
    }
    out <- list()
    out$df <- ns-2
    out$mswd <- as.numeric(X2/out$df)
    out$p.value <- as.numeric(1-stats::pchisq(X2,out$df))
    out
}

york.1966.zero.intercept <- function(x,alpha=0.05){
    colnames(x)[1:4] <- c('X','sX','Y','sY')
    b <- stats::lm(x[,'Y'] ~ 0 + x[,'X'])$coefficients # initial guess
    wX <- 1/x[,'sX']^2
    wY <- 1/x[,'sY']^2
    for (i in 1:10){
        W <- wX*wY/(wX+wY*b^2)
        Xbar <- sum(W*x[,'X'])/sum(W)
        Ybar <- sum(W*x[,'Y'])/sum(W)
        b <- Ybar/Xbar
    }
    U <- x[,'X']-Xbar
    V <- x[,'Y']-Ybar
    sb <- (sum(W*(b*U-V)^2)/sum(W*U^2))/(nrow(x)-1)
    c(b,sb)
}

# get fitted X and X given a dataset x=cbind(X,sX,Y,sY,rXY),
# an intercept a and slope b. This function is useful
# for evaluating log-likelihoods of derived quantities
get.york.xy <- function(x,a,b){
    wX <- 1/x[,'sX']^2
    wY <- 1/x[,'sY']^2
    A <- sqrt(wX*wY)
    W <- wX*wY/(wX+b*b*wY-2*b*x[,'rXY']*A)
    Xbar <- sum(W*x[,'X'],na.rm=TRUE)/sum(W,na.rm=TRUE)
    Ybar <- a + b*Xbar
    U <- x[,'X']-Xbar
    V <- x[,'Y']-Ybar
    B <- W*(U/wY+b*V/wX-(b*U+V)*x[,'rXY']/A)
    out <- cbind(Xbar+B,Ybar+b*B)
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
            out <- subset(x$x,select=c('Pb207U235','errPb207U235',
                                       'Pb206U238','errPb206U238','rhoXY'))
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
            out <- subset(x$x,select=c('U238Pb206','errU238Pb206',
                                       'Pb207Pb206','errPb207Pb206','rhoXY'))
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
data2york.UThHe <- function(x,...){
    ns <- length(x)
    out <- matrix(0,ns,5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    R <- iratio('U238U235')
    L8 <- lambda('U238')
    L5 <- lambda('U235')
    L2 <- lambda('Th232')
    L7 <- lambda('Sm147')
    f147 <- f147Sm()
    P <- rep(0,ns)
    sP <- rep(0,ns)
    J <- matrix(0,1,9)
    E <- matrix(0,9,9)
    for (i in 1:ns){
        P[i] <- 8*L8[1]*x[i,'U']*R[1]/(1+R[1]) +
            7*L5[1]*x[i,'U']/(1+R[1]) +
            6*L2[1]*x[i,'Th']
        J[1,1] <- 8*L8[1]*R[1]/(1+R[1]) + 7*L5[1]/(1+R[1])  # dP.dU
        J[1,2] <- 6*L2[1]                                   # dP.dTh
        J[1,4] <- 8*x[i,'U']*R[1]/(1+R[1])                  # dP.dL8
        J[1,5] <- 7*x[i,'U']/(1+R[1])                       # dP.dL5
        J[1,6] <- 6*x[i,'Th']                               # dP.dL2
        J[1,8] <- (8*L8[1]-7*L5[1])*x[i,'U']/(1+R[1])^2     # dP.dR
        E[1,1] <- x[i,'errU']^2
        E[2,2] <- x[i,'errTh']^2
        E[4,4] <- L8[2]^2
        E[5,5] <- L5[2]^2
        E[6,6] <- L2[2]^2
        E[7,7] <- L7[2]^2
        E[8,8] <- R[2]^2
        E[9,9] <- f147[2]^2
        if (doSm(x)) {
            P <- P + f147[1]*L7[1]*x[i,'Sm']
            J[1,3] <- f147[1]*L7[1]       # dP.dSm
            J[1,7] <- f147[1]*x[i,'Sm']   # dP.dL7
            J[1,9] <- L7[1]*x[i,'Sm']     # dP.df147
            E[3,3] <- x[i,'errSm']^2
        }
        sP[i] <- sqrt(J %*% E %*% t(J))
    }
    out[,'X'] <- P
    out[,'sX'] <- sP
    out[,'Y'] <- x[,'He']
    out[,'sY'] <- x[,'errHe']
    out
}
data2york.ThU <- function(x,type=2,...){
    if (x$format %in% c(1,3) & type==1){
        out <- subset(x$x,select=c('U238Th232','errU238Th232',
                                   'Th230Th232','errTh230Th232','rho'))
    } else if (x$format %in% c(2,4) & type==2){
        out <- subset(x$x,select=c('Th232U238','errTh232U238',
                                   'Th230U238','errTh230U238','rho'))
    } else if (x$format %in% c(2,4) & type==1){
        out <- ThConversionHelper(x)
        colnames(out) <- c('U238Th232','errU238Th232',
                           'Th230Th232','errTh230Th232','rho')
    } else if (x$format %in% c(1,3) & type==2){
        out <- ThConversionHelper(x)
        colnames(out) <- c('Th232U238','errTh232U238',
                           'Th230U238','errTh230U238','rho')
    } else {
        stop('Incorrect data format and/or plot type')
    }
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out
}

ThConversionHelper <- function(x){
    ns <- length(x)
    J <- matrix(0,2,2)
    E <- matrix(0,2,2)
    out <- matrix(0,ns,5)
    for (i in 1:ns){
        out[i,1] <- 1/x$x[i,1]
        out[i,3] <- x$x[i,3]/x$x[i,1]
        J[1,1] <- -out[i,1]/x$x[i,1]
        J[2,1] <- -out[i,3]/x$x[i,1]
        J[2,2] <- 1/x$x[i,1]
        E[1,1] <- x$x[i,2]^2
        E[2,2] <- x$x[i,4]^2
        E[1,2] <- x$x[i,2]*x$x[i,4]*x$x[i,5]
        E[2,1] <- E[1,2]
        covmat <- J %*% E %*% t(J)
        out[i,2] <- sqrt(covmat[1,1])
        out[i,4] <- sqrt(covmat[2,2])
        out[i,5] <- covmat[1,2]/(out[i,2]*out[i,4])
    }
    out
}
