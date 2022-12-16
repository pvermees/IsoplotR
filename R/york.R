#' @title
#' Linear regression of X,Y-variables with correlated errors
#'
#' @description
#' Implements the unified regression algorithm of York et al. (2004)
#' which, although based on least squares, yields results that are
#' consistent with maximum likelihood estimates of Titterington and
#' Halliday (1979).
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
#' @param x a 4 or 5-column matrix with the X-values, the analytical
#'     uncertainties of the X-values, the Y-values, the analytical
#'     uncertainties of the Y-values, and (optionally) the correlation
#'     coefficients of the X- and Y-values.
#' @return A seven-element list of vectors containing:
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
#'     \item{df}{degrees of freedom of the linear fit \eqn{(n-2)}}
#'
#'     \item{p.value}{p-value of a Chi-square value with \code{df}
#'     degrees of freedom}
#'
#'     }
#' 
#' @seealso \code{\link{data2york}}, \code{\link{titterington}},
#'     \code{\link{isochron}}, \code{\link{ludwig}}
#' 
#' @references
#' Titterington, D.M. and Halliday, A.N., 1979. On the fitting of
#' parallel isochrons and the method of maximum likelihood. Chemical
#' Geology, 26(3), pp.183-195.
#'
#' York, Derek, et al., 2004. Unified equations for the slope,
#' intercept, and standard errors of the best straight line.  American
#' Journal of Physics 72.3, pp.367-375.
#'
#' @examples
#' X <- c(1.550,12.395,20.445,20.435,20.610,24.900,
#'        28.530,50.540,51.595,86.51,106.40,157.35)
#' Y <- c(.7268,.7849,.8200,.8156,.8160,.8322,
#'        .8642,.9584,.9617,1.135,1.230,1.490)
#' n <- length(X)
#' sX <- X*0.01
#' sY <- Y*0.005
#' rXY <- rep(0.8,n)
#' dat <- cbind(X,sX,Y,sY,rXY)
#' fit <- york(dat)
#' scatterplot(dat,fit=fit)
#' @export
york <- function(x){
    if (ncol(x)==4) x <- cbind(x,0)
    colnames(x) <- c('X','sX','Y','sY','rXY')
    X <- x[,'X']
    Y <- x[,'Y']
    ab <- stats::lm(Y ~ X)$coefficients # initial guess
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
        U <- X-Xbar
        V <- Y-Ybar
        B <- W*(U/wY+b*V/wX-(b*U+V)*x[,'rXY']/A)
        b <- sum(W*B*V,na.rm=TRUE)/sum(W*B*U,na.rm=TRUE)
        if ((bold/b-1)^2 < 1e-15) break # convergence reached
    }
    a <- Ybar-b*Xbar
    xadj <- Xbar+B
    xbar <- sum(W*xadj,na.rm=TRUE)/sum(W,na.rm=TRUE)
    u <- xadj-xbar
    sb <- sqrt(1/sum(W*u^2,na.rm=TRUE))
    sa <- sqrt(1/sum(W,na.rm=TRUE)+(xbar*sb)^2)
    out <- list()
    out$a <- c(a,sa)
    out$b <- c(b,sb)
    out$cov.ab <- -xbar*sb^2
    names(out$a) <- c('a','s[a]')
    names(out$b) <- c('b','s[b]')
    out$type <- 'york'
    # compute MSWD:
    X2 <- sum(W*(Y-b*X-a)^2)
    out$df <- nrow(x)-2
    if (out$df>0){
        out$mswd <- as.numeric(X2/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(X2,out$df))
    } else {
        out$mswd <- 1
        out$p.value <- 1
    }
    out
}

# get fitted X and Y given a dataset x=cbind(X,sX,Y,sY,rXY),
# an intercept a and slope b. This function is useful
# for evaluating log-likelihoods of derived quantities
get.york.xy <- function(XY,a,b){
    X <- XY[,'X']
    sX <- XY[,'sX']
    Y <- XY[,'Y']
    sY <- XY[,'sY']
    sXY <- XY[,'rXY']*sX*sY
    O <- invertcovmat(sx=sX,sy=sY,sxy=sXY)
    N <- O[,'xx']*X + O[,'xy']*b*X + O[,'xy']*(Y-a) + b*(Y-a)*O[,'yy']
    D <- O[,'xx'] + 2*b*O[,'xy'] + O[,'yy']*b^2
    x <- N/D
    y <- a + b*x
    dNda <- - O[,'xy'] - b*O[,'yy']
    dNdb <- O[,'xy']*X + (Y-a)*O[,'yy']
    dDda <- 0
    dDdb <- 2*O[,'xy'] + 2*O[,'yy']*b
    dxda <- (dNda*D-N*dDda)/D^2
    dxdb <- (dNdb*D-N*dDdb)/D^2
    dyda <- 1 + b*dxda
    dydb <- x + b*dxdb
    cbind('x'=x,'y'=y,'dxda'=dxda,'dxdb'=dxdb,'dyda'=dyda,'dydb'=dydb)
}

#' @title Prepare geochronological data for York regression
#'
#' @description
#' Takes geochronology data as input and produces a five-column table
#' with the variables, their uncertainties and error correlations as
#' output. These can subsequently be used for York regression.
#'
#' @param x a five or six column matrix OR an object of class
#'     \code{UPb}, \code{PbPb}, \code{ThPb}, \code{ArAr}, \code{ThU},
#'     \code{UThHe}, or \code{PD} (which includes objects of class
#'     \code{RbSr}, \code{SmNd}, \code{LuHf} and \code{ReOs}),
#'     generated by the \code{read.data(...)} function
#' @param format one of
#'
#' \code{1} or \code{2}: \code{X}, \code{s[X]}, \code{Y}, \code{s[Y]},
#' \code{rho}; where \code{rho} is the error correlation between
#' \code{X} and \code{Y}; or
#'
#' \code{3}: \code{X/Z}, \code{s[X/Z]}, \code{Y/Z}, \code{s[Y/Z]},
#' \code{X/Y}, \code{s[X/Y]}; for which the error correlations are
#' automatically computed from the redundancy of the three ratios.
#'
#' @param ... optional arguments
#' 
#' @return a five-column table that can be used as input for
#'     \code{\link{york}}-regression.
#' 
#' @examples
#' f <- system.file("RbSr1.csv",package="IsoplotR")
#' dat <- read.csv(f)
#' yorkdat <- data2york(dat)
#' fit <- york(yorkdat)
#' 
#' @seealso \code{\link{york}}
#' @rdname data2york
#' @export
data2york <- function(x,...){ UseMethod("data2york",x) }
#' @rdname data2york
#' @export
data2york.default <- function(x,format=1,...){
    cnames <- c('X','sX','Y','sY','rXY')
    if (format==3){
        X <- cbind(x[,1:4],get.cor.div(x[,1],x[,2],x[,3],
                                       x[,4],x[,5],x[,6]))
        opt <- NULL
    } else {
        X <- x
        opt <- 5
    }
    insert.data(x=X,cnames=cnames,opt=opt)
}

#' @param option one of
#'
#' \code{1}: Wetherill concordia ratios: \code{X=07/35},
#' \code{sX=s[07/35]}, \code{Y=06/38}, \code{sY=s[06/38]}, \code{rXY}.
#'
#' \code{2}: Tera-Wasserburg ratios: \code{X=38/06},
#' \code{sX=s[38/06]}, \code{Y=07/06}, \code{sY=s[07/06]},
#' \code{rho=rXY}.
#'
#' \code{3}: \code{X=38/06}, \code{sX=s[38/06]}, \code{Y=04/06},
#' \code{sY=s[04/06]}, \code{rho=rXY} (only valid if \code{x$format=4,5}
#' or \code{6}).
#'
#' \code{4}: \code{X=35/07}, \code{sX=s[35/07]}, \code{Y=04/07},
#' \code{sY=s[04/07]}, \code{rho=rXY} (only valid if \code{x$format=4,5}
#' or \code{6}).
#'
#' \code{5}: U-Th-Pb concordia ratios: \code{X=06/38},
#' \code{sX=s[06/38]}, \code{Y=08/32}, \code{sY=s[08/32]},
#' \code{rho=rXY} (only valid if \code{x$format=7,8}).
#'
#' \code{6}: \code{X=38/06}, \code{sX=s[38/06]}, \code{Y=08/06},
#' \code{sY=s[08/06]}, \code{rho=rXY} (only valid if \code{x$format=7,8}).
#' 
#' \code{7}: \code{X=35/07}, \code{sX=s[35/07]}, \code{Y=08/07},
#' \code{sY=s[08/07]}, \code{rho=rXY} (only valid if \code{x$format=7,8}).
#' 
#' \code{8}: \code{X=32/08}, \code{sX=s[32/08]}, \code{Y=06/08},
#' \code{sY=s[06/08]}, \code{rho=rXY} (only valid if \code{x$format=7,8}).
#'
#' \code{9}: \code{X=32/08}, \code{sX=s[32/08]}, \code{Y=07/08},
#' \code{sY=s[07/08]}, \code{rho=rXY} (only valid if \code{x$format=7,8}).
#'
#' @param tt the age of the sample. This is only used if
#'     \code{x$format=7} or \code{8}, in order to calculate the
#'     inherited \eqn{{}^{208}}Pb/\eqn{{}^{232}}Th ratio.
#'
#' @rdname data2york
#' @export
data2york.UPb <- function(x,option=1,tt=0,...){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (option==1){ # 06/38 vs. 07/35
        for (i in 1:ns){
            wd <- wetherill(x,i=i)
            out[i,] <- data2york_UPb_helper(wd,i1='Pb207U235',i2='Pb206U238')
        }
    } else if (option==2){ # 07/06 vs. 38/06
        for (i in 1:ns){
            td <- tera.wasserburg(x,i=i)
            out[i,] <- data2york_UPb_helper(td,i1='U238Pb206',i2='Pb207Pb206')
        }
    } else if (option==3 && x$format%in%c(4,5,6)){ # 04/06 vs 38/06
        for (i in 1:ns){
            ir <- get.UPb.isochron.ratios.204(x,i)
            out[i,] <- data2york_UPb_helper(ir,i1='U238Pb206',i2='Pb204Pb206')
        }
    } else if (option==4 && x$format%in%c(4,5,6)){ # 04/07 vs. 35/07
        for (i in 1:ns){
            ir <- get.UPb.isochron.ratios.204(x,i)
            out[i,] <- data2york_UPb_helper(ir,i1='U235Pb207',i2='Pb204Pb207')
        }
    } else if (option==5 && x$format%in%c(7,8)){ # 08/32 vs. 06/38
        for (i in 1:ns){
            wd <- wetherill(x,i=i)
            out[i,] <- data2york_UPb_helper(wd,i1='Pb206U238',i2='Pb208Th232')
        }        
    } else if (option==6 && x$format%in%c(7,8)){ # 08/06 vs. 38/06
        for (i in 1:ns){
            ir <- get.UPb.isochron.ratios.208(x,i,tt=tt)
            out[i,] <- data2york_UPb_helper(ir,i1='U238Pb206',i2='Pb208cPb206')
        }
    } else if (option==7 && x$format%in%c(7,8)){ # 08/07 vs. 35/07
        for (i in 1:ns){
            ir <- get.UPb.isochron.ratios.208(x,i,tt=tt)
            out[i,] <- data2york_UPb_helper(ir,i1='U235Pb207',i2='Pb208cPb207')
        }
    } else if (option==8 && x$format%in%c(7,8)){ # 06/08 vs. 32/08
        for (i in 1:ns){
            ir <- get.UPb.isochron.ratios.208(x,i,tt=tt)
            out[i,] <- data2york_UPb_helper(ir,i1='Th232Pb208',i2='Pb206cPb208')
        }
    } else if (option==9 && x$format%in%c(7,8)){ # 07/08 vs. 32/08
        for (i in 1:ns){
            ir <- get.UPb.isochron.ratios.208(x,i,tt=tt)
            out[i,] <- data2york_UPb_helper(ir,i1='Th232Pb208',i2='Pb207cPb208')
        }
    } else {
        stop('Incompatible input format and concordia type.')
    }
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out
}
data2york_UPb_helper <- function(x,i1=1,i2=2){
    X <- x$x[i1]
    sX <- sqrt(x$cov[i1,i1])
    Y <- x$x[i2]
    sY <- sqrt(x$cov[i2,i2])
    rXY <- x$cov[i1,i2]/(sX*sY)
    c(X,sX,Y,sY,rXY)
}
#' @param inverse toggles between normal and inverse isochron
#'     ratios. \code{data2york} returns five columns \code{X},
#'     \code{s[X]}, \code{Y}, \code{s[Y]} and \code{r[X,Y]}.
#'
#' If \code{inverse=TRUE}, then \code{X} =
#' \eqn{{}^{204}}Pb/\eqn{{}^{206}}Pb and \code{Y} =
#' \eqn{{}^{207}}Pb/\eqn{{}^{206}}Pb (if \code{x} has class
#' \code{PbPb}), or \code{X} = \eqn{{}^{232}}Th/\eqn{{}^{208}}Pb and
#' \code{Y} = \eqn{{}^{204}}Pb/\eqn{{}^{208}}Pb (if \code{x} has class
#' \code{ThPb}), or \code{X} = \eqn{{}^{39}}Ar/\eqn{{}^{40}}Ar and
#' \code{Y} = \eqn{{}^{36}}Ar/\eqn{{}^{40}}Ar (if \code{x} has class
#' \code{ArAr}), or \code{X} = \eqn{{}^{40}}K/\eqn{{}^{40}}Ca and
#' \code{Y} = \eqn{{}^{44}}Ca/\eqn{{}^{40}}Ca (if \code{x} has class
#' \code{KCa}), or \code{X} = \eqn{{}^{87}}Rb/\eqn{{}^{87}}Sr and
#' \code{Y} = \eqn{{}^{86}}Sr/\eqn{{}^{87}}Sr (if \code{x} has class
#' \code{RbSr}), or \code{X} = \eqn{{}^{147}}Sm/\eqn{{}^{143}}Nd and
#' \code{Y} = \eqn{{}^{144}}Nd/\eqn{{}^{143}}Nd (if \code{x} has class
#' \code{SmNd}), or \code{X} = \eqn{{}^{187}}Re/\eqn{{}^{187}}Os and
#' \code{Y} = \eqn{{}^{188}}Os/\eqn{{}^{187}}Os (if \code{x} has class
#' \code{ReOs}), or \code{X} = \eqn{{}^{176}}Lu/\eqn{{}^{176}}Hf and
#' \code{Y} = \eqn{{}^{177}}Hf/\eqn{{}^{176}}Hf (if \code{x} has class
#' \code{LuHf}).
#' 
#' If \code{inverse=FALSE}, then \code{X} =
#' \eqn{{}^{206}}Pb/\eqn{{}^{204}}Pb and \code{Y} =
#' \eqn{{}^{207}}Pb/\eqn{{}^{204}}Pb (if \code{x} has class
#' \code{PbPb}), or \code{X} = \eqn{{}^{232}}Th/\eqn{{}^{204}}Pb and
#' \code{Y} = \eqn{{}^{208}}Pb/\eqn{{}^{204}}Pb (if \code{x} has class
#' \code{ThPb}), or \code{X} = \eqn{{}^{39}}Ar/\eqn{{}^{36}}Ar and
#' \code{Y} = \eqn{{}^{40}}Ar/\eqn{{}^{36}}Ar (if \code{x} has class
#' \code{ArAr}), or \code{X} = \eqn{{}^{40}}K/\eqn{{}^{44}}Ca and
#' \code{Y} = \eqn{{}^{40}}Ca/\eqn{{}^{44}}Ca (if \code{x} has class
#' \code{KCa}), or \code{X} = \eqn{{}^{87}}Rb/\eqn{{}^{86}}Sr and
#' \code{Y} = \eqn{{}^{87}}Sr/\eqn{{}^{86}}Sr (if \code{x} has class
#' \code{RbSr}), or \code{X} = \eqn{{}^{147}}Sm/\eqn{{}^{144}}Nd and
#' \code{Y} = \eqn{{}^{143}}Nd/\eqn{{}^{144}}Nd (if \code{x} has class
#' \code{SmNd}), or \code{X} = \eqn{{}^{187}}Re/\eqn{{}^{188}}Os and
#' \code{Y} = \eqn{{}^{187}}Os/\eqn{{}^{188}}Os (if \code{x} has class
#' \code{ReOs}), or \code{X} = \eqn{{}^{176}}Lu/\eqn{{}^{177}}Hf and
#' \code{Y} = \eqn{{}^{176}}Hf/\eqn{{}^{177}}Hf (if \code{x} has class
#' \code{LuHf}).
#'
#' @rdname data2york
#' @export
data2york.ArAr <- function(x,inverse=TRUE,...){
    out <- data2york(x$x,format=x$format,...)
    invert <- (inverse & x$format==1) | (!inverse & x$format%in%c(2,3))
    if (invert) out <- normal2inverse(out)
    out
}
#' @rdname data2york
#' @export
data2york.ThPb <- function(x,inverse=FALSE,...){
    out <- data2york(x$x,format=x$format,...)
    invert <- (inverse & x$format==1) | (!inverse & x$format%in%c(2,3))
    if (invert) out <- normal2inverse(out)
    out
}
#' @rdname data2york
#' @export
data2york.KCa <- function(x,inverse=FALSE,...){
    out <- data2york(x$x,format=x$format,...)
    invert <- (inverse & x$format%in%c(1,3)) | (!inverse & x$format==2)
    if (invert) out <- normal2inverse(out)
    out
}
#' @rdname data2york
#' @export
data2york.PbPb <- function(x,inverse=TRUE,...){
    out <- data2york(x$x,format=x$format,...)
    invert <- (inverse & x$format%in%c(1,3)) | (!inverse & x$format==2)
    if (invert){ # swap columns for normal2inverse function
        out[,c('X','sX','Y','sY','rXY')] <- out[,c('Y','sY','X','sX','rXY')]
        out <- normal2inverse(out)
        out[,c('X','sX','Y','sY','rXY')] <- out[,c('Y','sY','X','sX','rXY')]
    }
    out
}

#' @param exterr If \code{TRUE}, propagates the external uncertainties
#'     (e.g. decay constants) into the output errors.
#' @rdname data2york
#' @export
data2york.PD <- function(x,exterr=FALSE,inverse=FALSE,...){
    if (x$format<3){
        X <- x$x
        format <- x$format
    } else {
        X <- ppm2ratios(x,exterr=exterr)
        format <- 1
    }
    out <- data2york(X,format=format,...)
    invert <- (inverse & format==1) | (!inverse & format==2)
    if (invert) out <- normal2inverse(out)
    out
}
#' @rdname data2york
#' @export
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
#' @param type Return `Rosholt' or `Osmond' ratios?
#'
#' Rosholt (\code{type=1}) returns \code{X=8/2}, \code{sX=s[8/2]},
#' \code{Y=0/2}, \code{sY=s[0/2]}, \code{rXY}.
#'
#' Osmond (\code{type=2}) returns \code{X=2/8}, \code{sX=s[2/8]},
#' \code{Y=0/8}, \code{sY=s[0/8]}, \code{rXY}.
#'
#' @param generic If \code{TRUE}, uses the following column headers:
#'     \code{X}, \code{sX}, \code{Y}, \code{sY}, \code{rXY}.
#'
#' If \code{FALSE} and \code{type=1}, uses \code{U238Th232},
#'     \code{errU238Th232}, \code{Th230Th232}, \code{errTh230Th232},
#'     \code{rho}
#'
#' or if \code{FALSE} and \code{type=2}, uses \code{Th232U238},
#'     \code{errTh232U238}, \code{Th230U238}, \code{errTh230U238},
#'     \code{rho}.
#' @rdname data2york
#' @export
data2york.ThU <- function(x,type=2,generic=TRUE,...){
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
    if (generic) colnames(out) <- c('X','sX','Y','sY','rXY')
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

normal2inverse <- function(x){
    out <- x
    iX <- 1
    isX <- 2
    iY <- 3
    isY <- 4
    irXY <- 5
    out[,iX] <- x[,iX]/x[,iY]
    out[,iY] <- 1/x[,iY]
    E11 <- x[,isX]^2
    E22 <- x[,isY]^2
    E12 <- x[,irXY]*x[,isX]*x[,isY]
    J11 <- 1/x[,iY]
    J12 <- -out[,iX]/x[,iY]
    J21 <- rep(0,nrow(x))
    J22 <- -out[,iY]/x[,iY]
    err <- errorprop(J11,J12,J21,J22,E11,E22,E12)
    out[,isX] <- sqrt(err[,'varX'])
    out[,isY] <- sqrt(err[,'varY'])
    out[,irXY] <- err[,'cov']/(out[,isX]*out[,isY])
    out
}
