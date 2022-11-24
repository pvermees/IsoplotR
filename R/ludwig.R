#' @title
#' Linear regression of U-Pb data with correlated errors, taking
#' into account decay constant uncertainties.
#'
#' @description
#' Implements the maximum likelihood algorithm for Total-Pb/U isochron
#' regression of Ludwig (1998) and extends the underlying methodology
#' to accommodate U-Th-Pb data and initial U-series disequilibrium.
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
#' ages.
#'
#' @param x an object of class \code{UPb}
#' 
#' @param exterr propagate external sources of
#' uncertainty (i.e. decay constants)?
#' 
#' @param model one of three regression models:
#'
#' \code{1}: fit a discordia line through the data using the maximum
#' likelihood algorithm of Ludwig (1998), which assumes that the
#' scatter of the data is solely due to the analytical
#' uncertainties. In this case, \code{IsoplotR} will either calculate
#' an upper and lower intercept age (for Wetherill concordia), or a
#' lower intercept age and common
#' (\eqn{^{207}}Pb/\eqn{^{206}}Pb)\eqn{_\circ}-ratio intercept (for
#' Tera-Wasserburg). If the p-value for the chi-square test is less
#' than \code{alpha()}, then the analytical uncertainties are augmented
#' by a factor \eqn{\sqrt{MSWD}}.
#'
#' \code{2}: fit a discordia line ignoring the analytical uncertainties
#'
#' \code{3}: fit a discordia line using a modified maximum likelihood
#' algorithm that includes accounts for any overdispersion by adding a
#' geological (co)variance term.
#' 
#' @param anchor
#' control parameters to fix the intercept age or common Pb
#' composition of the isochron fit. This can be a scalar or a vector.
#'
#' If \code{anchor[1]=0}: do not anchor the isochron.
#'
#' If \code{anchor[1]=1}: fix the common Pb composition at the values
#' stored in \code{settings('iratio',...)}.
#'
#' If \code{anchor[1]=2}: force the isochron line to intersect the
#' concordia line at an age equal to \code{anchor[2]}.
#'
#' @param ... optional arguments
#' 
#' @return
#' \describe{
#'
#' \item{par}{a vector with the lower concordia intercept, the common
#' Pb ratios, the dispersion parameter (if \code{model=3}), and the
#' initial \eqn{{}^{234}}U/\eqn{{}^{238}}U and
#' \eqn{{}^{230}}Th/\eqn{{}^{238}}U activity ratio (in the presence of
#' initial disequilibrium).}
#'
#' \item{cov}{the covariance matrix of \code{par}}
#'
#' \item{df}{the degrees of freedom of the model fit (\eqn{n-2} if
#' \code{x$format<4} or \eqn{2n-3} if \code{x$format>3}, where \eqn{n}
#' is the number of aliquots).}
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
#' 
#' @references
#' Ludwig, K.R., 1998. On the treatment of concordant uranium-lead
#' ages. Geochimica et Cosmochimica Acta, 62(4), pp.665-676.
#'
#' Ludwig, K.R. and Titterington, D.M., 1994. Calculation of
#' \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#' Cosmochimica Acta, 58(22), pp.5031-5042.
#'
#' @seealso \code{\link{concordia}}, \code{\link{titterington}},
#'     \code{\link{isochron}}
#' @rdname ludwig
#' @export
ludwig <- function(x,...){ UseMethod("ludwig",x) }
#' @rdname ludwig
#' @export
ludwig <- function(x,model=1,anchor=0,exterr=FALSE,joint=TRUE,type=1,...){
    init <- init.ludwig(x,model=model,anchor=anchor)
    fit1 <- stats::optim(par=init$pars,fn=LL.ludwig,method="L-BFGS-B",
                         lower=init$lower,upper=init$upper,hessian=FALSE,
                         x=x,anchor=anchor,model=model,exterr=exterr)
    fit2 <- fit1
    am1 <- anchormerge(fit1$par,x,anchor=anchor)
    fit2$par <- am1$x
    hess <- stats::optimHess(par=fit2$par,fn=LL.ludwig,
                             x=x,anchor=anchor,model=model,
                             exterr=exterr,hessian=TRUE)
    if (det(hess)<0){
        warning('Ill-conditioned Hessian, replaced by ',
                'nearest positive definite matrix')
        hess <- nearPD(hess)
    }
    fit2$cov <- MASS::ginv(hess)
    dimnames(fit2$cov) <- dimnames(hess)
    am2 <- anchormerge(fit1$par,x,anchor=anchor,dontchecksx=TRUE)
    anchored <- names(am2$x)[names(am2$x) %ni% names(fit1$par)]
    out <- fit2
    out$par <- am2$x
    out$cov <- matrix(0,length(am2$x),length(am2$x))
    rownames(out$cov) <- colnames(out$cov) <- names(am2$x)
    out$cov[names(fit2$par),names(fit2$par)] <- fit2$cov
    out$cov[anchored,anchored] <- 0
    diag(out$cov)[anchored] <- am2$sx[anchored]^2
    mswd <- mswd.lud(out$par,x=x,anchor=anchor,exterr=exterr)
    out$model <- model
    c(out,mswd)
}

anchormerge <- function(pars,x,anchor=0,dontchecksx=FALSE){
    X <- pars
    sX <- 0*pars
    if (anchor[1]==2){
        nonzerosx <- (anchor[1]==2 && length(anchor)>2 && anchor[3]>0)
        if (dontchecksx || nonzerosx){
            X['t'] <- anchor[2]
            sX['t'] <- ifelse(length(anchor)>2,anchor[3],0)
        }
    }
    if (anchor[1]==1){
        if (x$format<4){
            nonzerosx <- (iratio('Pb207Pb206')[2]>0)
            if (dontchecksx || nonzerosx){
                X['a0'] <- iratio('Pb207Pb206')[1]
                sX['a0'] <- iratio('Pb207Pb206')[2]
            }
        } else if (x$format<7){
            nonzerosx <- (iratio('Pb206Pb204')[2]>0 && iratio('Pb207Pb204')[2]>0)
            if (dontchecksx || nonzerosx){
                X['a0'] <- iratio('Pb206Pb204')[1]
                X['b0'] <- iratio('Pb207Pb204')[1]
                sX['a0'] <- iratio('Pb206Pb204')[2]
                sX['b0'] <- iratio('Pb207Pb204')[2]
            }
        } else {
            nonzerosx <- (iratio('Pb208Pb206')[2]>0 && iratio('Pb208Pb207')[2]>0)
            if (dontchecksx || nonzerosx){
                X['a0'] <- 1/iratio('Pb208Pb206')[1]
                X['b0'] <- 1/iratio('Pb208Pb207')[1]
                sX['a0'] <- iratio('Pb208Pb206')[2]*X['a0']^2
                sX['b0'] <- iratio('Pb208Pb207')[2]*X['b0']^2
            }
        }
    }
    checkit <- function(x,dontchecksx=FALSE){
        x$option==1 && (dontchecksx || (!is.null(x$sx) && x$sx>0))
    }
    if (checkit(x$d$U48,dontchecksx)){
        X['U48i'] <- x$d$U48$x
        sX['U48i'] <- x$d$U48$sx
    }
    if (checkit(x$d$ThU,dontchecksx)){
        X['ThUi'] <- x$d$ThU$x
        sX['ThUi'] <- x$d$ThU$sx
    }
    if (checkit(x$d$RaU,dontchecksx)){
        X['RaUi'] <- x$d$RaU$x
        sX['RaUi'] <- x$d$RaU$sx
    }
    if (checkit(x$d$PaU,dontchecksx)){
        X['PaUi'] <- x$d$PaU$x
        sX['PaUi'] <- x$d$PaU$sx
    }
    sortednames <- c('t','a0','b0','w','U48i','ThUi','RaUi','PaUi')
    pnames <- names(X)
    snames <- sortednames[which(sortednames %in% pnames)]
    out <- list()
    out$x <- X[snames]
    out$sx <- sX[snames]
    out
}

init.ludwig <- function(x,model=1,anchor=0,joint=TRUE,type=1){
    pars <- lower <- upper <- vector()
    if (x$format<4){
        yd <- data2york(x,option=2)
        yfit <- york(yd)
        PbU0 <- abs(yfit$b[1]/yfit$a[1])
        if (anchor[1]==1){
            pars['t'] <- get.Pb206U238.age(x=PbU0)[1]
            lower['t'] <- pars['t']/10
            upper['t'] <- pars['t']*10
        } else if (anchor[1]==2){
            pars['a0'] <- yfit$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
        } else {
            pars['t'] <- get.Pb206U238.age(x=PbU0)[1]
            lower['t'] <- pars['t']/10
            upper['t'] <- pars['t']*10            
            pars['a0'] <- yfit$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10            
        }
    } else if (x$format<7){
        yda <- data2york(x,option=3)
        yfita <- york(yda)
        Pb6U8 <- abs(yfita$b[1]/yfita$a[1])
        ydb <- data2york(x,option=4)
        yfitb <- york(ydb)
        if (anchor[1]==1){
            pars['t'] <- get.Pb206U238.age(x=Pb6U8)[1]
            lower['t'] <- pars['t']/10
            upper['t'] <- pars['t']*10
        } else if (anchor[1]==2){
            pars['a0'] <- 1/yfita$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
            pars['b0'] <- 1/yfitb$a[1]
            lower['b0'] <- pars['b0']/10
            upper['b0'] <- pars['b0']*10
        } else {
            pars['t'] <- get.Pb206U238.age(x=Pb6U8)[1]
            lower['t'] <- pars['t']/10
            upper['t'] <- pars['t']*10
            pars['a0'] <- 1/yfita$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
            pars['b0'] <- 1/yfitb$a[1]
            lower['b0'] <- pars['b0']/10
            upper['b0'] <- pars['b0']*10
        }
    } else {
        yd <- data2york(x,option=2)
        yfit <- york(yd)
        Pb6U8 <- abs(yfit$b[1]/yfit$a[1])
        tt <- get.Pb206U238.age(x=Pb6U8)[1]
        yda <- data2york(x,option=6,tt=tt)
        yfita <- york(yda)
        ydb <- data2york(x,option=7,tt=tt)
        yfitb <- york(ydb)
        if (anchor[1]==1){
            pars['t'] <- tt
            lower['t'] <- tt/10
            upper['t'] <- tt*10
        } else if (anchor[1]==2){
            pars['a0'] <- 1/yfita$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
            pars['b0'] <- 1/yfitb$a[1]
            lower['b0'] <- pars['b0']/10
            upper['b0'] <- pars['b0']*10
        } else {
            pars['t'] <- tt
            lower['t'] <- tt/10
            upper['t'] <- tt*10
            pars['a0'] <- 1/yfita$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
            pars['b0'] <- 1/yfitb$a[1]
            lower['b0'] <- pars['b0']/10
            upper['b0'] <- pars['b0']*10
        }
    }
    if (model==3){
        pars['w'] <- 0.01
        lower['w'] <- 0
        upper['w'] <- 10
    }
    if (x$d$U48$option==2){
        if (x$d$U48$sx>0){
            pars['U48i'] <- 1
            lower['U48i'] <- 0
            upper['U48i'] <- 20
        } else {
            stop('Zero uncertainty of measured 234/238 activity ratio')
        }
    }
    if (x$d$ThU$option==2){
        if (x$d$ThU$sx>0){
            pars['ThUi'] <- 1
            lower['ThUi'] <- 0
            upper['ThUi'] <- 20
        } else {
            stop('Zero uncertainty of measured 230/238 activity ratio')
        }
    }
    list(pars=pars,lower=lower,upper=upper)
}

LL.ludwig <- function(pars,x,model=1,exterr=FALSE,anchor=0,hessian=FALSE){
    X <- x
    if ('U48i' %in% names(pars)){
        X$d$U48$x <- pars['U48i']
        X$d$U48$option <- 1
    }
    if ('ThUi' %in% names(pars)){
        X$d$ThU$x <- pars['ThUi']
        X$d$ThU$option <- 1
    }
    if ('RaUi' %in% names(pars)){
        X$d$RaU$x <- pars['RaUi']
        X$d$RaU$option <- 1
    }
    if ('PaUi' %in% names(pars)){
        X$d$PaU$x <- pars['PaUi']
        X$d$PaU$option <- 1
    }
    LL <- 0
    if (x$format<4){
        if (anchor[1]==1){
            tt <- pars['t']
            if (hessian && iratio('Pb207Pb206')[2]>0){
                a0 <- pars['a0']
                LL <- LL - stats::dnorm(x=a0,
                                        mean=iratio('Pb207Pb206')[1],
                                        sd=iratio('Pb207Pb206')[2],
                                        log=TRUE)
            } else {
                a0 <- iratio('Pb207Pb206')[1]
            }
        } else if (anchor[1]==2){
            a0 <- pars['a0']
            if (hessian && length(anchor)>2 && anchor[3]>0){
                tt <- pars['t']
                st <- anchor[3]
                LL <- LL - stats::dnorm(x=tt,mean=anchor[2],sd=st,log=TRUE)
            } else {
                tt <- anchor[2]
            }
        } else {
            tt <- pars['t']
            a0 <- pars['a0']
        }
        ta0b0w <- c('t'=unname(tt),'a0'=unname(a0))
    } else if (x$format<7){
        if (anchor[1]==1){
            tt <- pars['t']
            if (hessian && iratio('Pb206Pb204')[2]>0){
                a0 <- pars['a0']
                LL <- LL - stats::dnorm(x=a0,
                                        mean=iratio('Pb206Pb204')[1],
                                        sd=iratio('Pb206Pb204')[2],
                                        log=TRUE)
            } else {
                a0 <- iratio('Pb206Pb204')[1]
            }
            if (hessian && iratio('Pb207Pb204')[2]>0){
                b0 <- pars['b0']
                LL <- LL - stats::dnorm(x=b0,
                                        mean=iratio('Pb207Pb204')[1],
                                        sd=iratio('Pb207Pb204')[2],
                                        log=TRUE)
            } else {
                b0 <- iratio('Pb207Pb204')[1]
            }
        } else if (anchor[1]==2){
            a0 <- pars['a0']
            b0 <- pars['b0']
            if (hessian && length(anchor)>2 && anchor[3]>0){
                tt <- pars['t']
                st <- anchor[3]
                LL <- LL - stats::dnorm(x=tt,mean=anchor[2],sd=st,log=TRUE)
            } else {
                tt <- anchor[2]
            }
        } else {
            tt <- pars['t']
            a0 <- pars['a0']
            b0 <- pars['b0']
        }
        ta0b0w <- c('t'=unname(tt),'a0'=unname(a0),'b0'=unname(b0))
    } else {
        if (anchor[1]==1){
            tt <- pars['t']
            if (hessian && iratio('Pb208Pb206')[2]>0){
                a0 <- pars['a0']
                LL <- LL - stats::dnorm(x=1/a0,
                                        mean=iratio('Pb208Pb206')[1],
                                        sd=iratio('Pb208Pb206')[2],
                                        log=TRUE)
            } else {
                a0 <- 1/iratio('Pb208Pb206')[1]
            }
            if (hessian && iratio('Pb208Pb207')[2]>0){
                b0 <- pars['b0']
                LL <- LL - stats::dnorm(x=1/b0,
                                        mean=iratio('Pb208Pb207')[1],
                                        sd=iratio('Pb208Pb207')[2],
                                        log=TRUE)
            } else {
                b0 <- 1/iratio('Pb208Pb207')[1]
            }
        } else if (anchor[1]==2){
            a0 <- pars['a0']
            b0 <- pars['b0']
            if (hessian && length(anchor)>2 && anchor[3]>0){
                tt <- pars['t']
                st <- anchor[3]
                LL <- LL - stats::dnorm(x=tt,mean=anchor[2],sd=st,log=TRUE)
            } else {
                tt <- anchor[2]
            }
        } else {
            tt <- pars['t']
            a0 <- pars['a0']
            b0 <- pars['b0']
        }
        ta0b0w <- c('t'=unname(tt),'a0'=unname(a0),'b0'=unname(b0))
    }
    if (model==3){
        ta0b0w['w'] <- pars['w']
    }
    if (model==2){
        LL <- LL + LL.ludwig.model2(ta0b0w,X,anchor=anchor,exterr=exterr)
    } else {
        LL <- LL + data2ludwig(X,ta0b0w,anchor=anchor,exterr=exterr)$LL
    }
    if (x$d$U48$option==2){
        pred <- mclean(tt=tt,d=X$d)
        LL <- LL - stats::dnorm(x=pred$U48,mean=x$d$U48$x,sd=x$d$U48$sx,log=TRUE)
    } else if (hessian && x$d$U48$option==1 && x$d$U48$sx>0){
        U48i <- pars['U48i']
        LL <- LL - stats::dnorm(x=U48i,mean=x$d$U48$x,sd=x$d$U48$sx,log=TRUE)        
    }
    if (x$d$ThU$option==2){
        pred <- mclean(tt=tt,d=X$d)
        LL <- LL - stats::dnorm(x=pred$ThU,mean=x$d$ThU$x,sd=x$d$ThU$sx,log=TRUE)        
    } else if (hessian && x$d$ThU$option==1 && x$d$ThU$sx>0){
        ThUi <- pars['ThUi']
        LL <- LL - stats::dnorm(x=ThUi,mean=x$d$ThU$x,sd=x$d$ThU$sx,log=TRUE)
    }
    if (hessian && x$d$RaU$option==1 && x$d$RaU$sx>0){
        RaUi <- pars['RaUi']
        LL <- LL - stats::dnorm(x=RaUi,mean=x$d$RaU$x,sd=x$d$RaU$sx,log=TRUE)
    }
    if (hessian && x$d$PaU$option==1 && x$d$PaU$sx>0){
        PaUi <- pars['PaUi']
        LL <- LL - stats::dnorm(x=PaUi,mean=x$d$PaU$x,sd=x$d$PaU$sx,log=TRUE)
    }
    LL
}

data2ludwig <- function(x,ta0b0w,anchor=0,exterr=FALSE){
    out <- list()
    U <- iratio('U238U235')[1]
    tt <- ifelse(anchor[1]==2,anchor[2],ta0b0w['t'])
    disp <- 'w' %in% names(ta0b0w)
    ns <- length(x)
    zeros <- rep(0,ns)
    X <- zeros
    Y <- zeros
    K0 <- zeros
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    if (x$format%in%c(1,2,3)){
        a0 <- ifelse(anchor[1]==1,iratio('Pb207Pb206')[1],ta0b0w['a0'])
        NP <- 2 # number of fit parameters (tt, a0)
        NR <- 2 # number of isotopic ratios (X, Y)
    } else if (x$format%in%c(4,5,6)){
        a0 <- ifelse(anchor[1]==1,iratio('Pb206Pb204')[1],ta0b0w['a0'])
        b0 <- ifelse(anchor[1]==1,iratio('Pb207Pb204')[1],ta0b0w['b0'])
        Z <- zeros
        L0 <- zeros
        NP <- 3 # tt, a0, b0
        NR <- 3 # X, Y, Z
    } else if (x$format%in%c(7,8)){
        a0 <- ifelse(anchor[1]==1,1/iratio('Pb208Pb206')[1],ta0b0w['a0'])
        b0 <- ifelse(anchor[1]==1,1/iratio('Pb208Pb207')[1],ta0b0w['b0'])
        Z <- zeros
        W <- zeros
        L0 <- zeros
        NP <- 3 # tt, a0, b0
        NR <- 4 # X, Y, Z, W
    } else {
        stop('Incorrect input format.')
    }
    if (disp) w <- ta0b0w['w']
    E <- matrix(0,NR*ns+7,NR*ns+7)
    J <- matrix(0,NP*ns,NR*ns+7)
    J[1:(NP*ns),1:(NP*ns)] <- diag(NP*ns)
    nc <- length(D$ThUi) # nc>1 if each aliquot has its own diseq correction
    j <- 1
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        if (x$format%in%c(4,5,6)){
            Z[i] <- wd$x['Pb204U238']
        } else if (x$format>6){
            Z[i] <- wd$x['Pb208Th232']
            W[i] <- wd$x['Th232U238']
        }
        E[(0:(NR-1))*ns+i,(0:(NR-1))*ns+i] <- wd$cov
        if (nc>1) j <- i
        J[i,NR*ns+2] <- -D$dPb207U235dl35[j]     #dKdl35
        J[i,NR*ns+5] <- -D$dPb207U235dl31[j]     #dKdl31
        J[ns+i,NR*ns+1] <- -D$dPb206U238dl38[j]  #dLdl38
        J[ns+i,NR*ns+3] <- -D$dPb206U238dl34[j]  #dLdl34
        J[ns+i,NR*ns+6] <- -D$dPb206U238dl30[j]  #dLdl30
        J[ns+i,NR*ns+7] <- -D$dPb206U238dl26[j]  #dLdl26
        if (x$format>6) J[2*ns+i,NR*ns+4] <- -D$dPb208Th232dl32[j] # dMdl32
    }
    E[NR*ns+1:7,NR*ns+1:7] <- getEl()
    ED <- J%*%E%*%t(J)
    if (disp){ # fit overdispersion
        Ew <- get.Ewd(w=w,format=x$format,ns=ns,D=D)
        ED <- ED + Ew
    }
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    if (x$format<4){
        O <- blockinverse(AA=ED[i1,i1],BB=ED[i1,i2],
                          CC=ED[i2,i1],DD=ED[i2,i2],doall=TRUE)
    } else {
        i3 <- (2*ns+1):(3*ns)
        O <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                             DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                             GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    }
    if (x$format%in%c(1,2,3)){
        K0 <- X - D$Pb207U235 + a0*U*(D$Pb206U238 - Y)
        A <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*a0*U +
               K0%*%(O[i1,i2]+t(O[i2,i1])))
        B <- -(a0*U*(O[i1,i1]+t(O[i1,i1]))*a0*U +
               (O[i2,i2]+t(O[i2,i2])) +
               a0*U*(O[i1,i2]+t(O[i1,i2])) +
               (O[i2,i1]+t(O[i2,i1]))*a0*U)
        L <- as.vector(solve(B,A))
        c0 <- Y - D$Pb206U238 - L
        K <- X - D$Pb207U235 - a0*U*c0
        KLM <- c(K,L)
    } else if (x$format%in%c(4,5,6)){
        K0 <- X - D$Pb207U235 - U*b0*Z
        L0 <- Y - D$Pb206U238 - a0*Z
        V <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*U*b0 +
               L0%*%(O[i1,i2]+t(O[i2,i1]))*U*b0 +
               K0%*%(O[i1,i2]+t(O[i2,i1]))*a0 +
               L0%*%(O[i2,i2]+t(O[i2,i2]))*a0 +
               K0%*%(O[i1,i3]+t(O[i3,i1])) +
               L0%*%(O[i2,i3]+t(O[i3,i2])))
        W <- -(U*b0*(O[i1,i1]+t(O[i1,i1]))*U*b0 +
               U*b0*(O[i1,i2]+t(O[i1,i2]))*a0 +
               U*b0*(O[i1,i3]+t(O[i1,i3])) +
               a0*(O[i2,i1]+t(O[i2,i1]))*U*b0 +
               a0*(O[i2,i2]+t(O[i2,i2]))*a0 +
               a0*(O[i2,i3]+t(O[i2,i3])) +
               (O[i3,i1]+t(O[i3,i1]))*U*b0 +
               (O[i3,i2]+t(O[i3,i2]))*a0 +
               (O[i3,i3]+t(O[i3,i3])))
        M <- as.vector(solve(W,V))
        c0 <- as.vector(Z - M)
        K <- as.vector(X - D$Pb207U235 - U*b0*c0)
        L <- as.vector(Y - D$Pb206U238 - a0*c0)
        KLM <- c(K,L,M)
    } else if (x$format%in%c(7,8)){
        Wd <- diag(W)
        K0 <- X - D$Pb207U235 - (Z-D$Pb208Th232)*U*W*b0
        L0 <- Y - D$Pb206U238 - (Z-D$Pb208Th232)*W*a0
        AA <- (Wd%*%O[i1,i1]%*%Wd)*(U*b0)^2 +
            (Wd%*%O[i2,i2]%*%Wd)*a0^2 + O[i3,i3] +
            U*a0*b0*Wd%*%(O[i1,i2]+O[i2,i1])%*%Wd +
            U*b0*(Wd%*%O[i1,i3]+O[i3,i1]%*%Wd) +
            a0*(Wd%*%O[i2,i3]+O[i3,i2]%*%Wd)
        BT <- t(U*b0*K0%*%O[i1,i1]%*%Wd +
                a0*L0%*%O[i2,i2]%*%Wd +
                a0*K0%*%O[i1,i2]%*%Wd +
                U*b0*L0%*%O[i2,i1]%*%Wd +
                K0%*%O[i1,i3] + L0%*%O[i2,i3])
        CC <- U*b0*Wd%*%O[i1,i1]%*%K0 +
            a0*Wd%*%O[i2,i2]%*%L0 +
            a0*Wd%*%O[i2,i1]%*%K0 +
            U*b0*Wd%*%O[i1,i2]%*%L0 +
            O[i3,i1]%*%K0 + O[i3,i2]%*%L0
        M <- as.vector(solve(-(AA+t(AA)),(BT+CC)))
        c0 <- as.vector(Z - D$Pb208Th232 - M)
        K <- as.vector(X - D$Pb207U235 - c0*b0*U*W)
        L <- as.vector(Y - D$Pb206U238 - c0*a0*W)
        KLM <- c(K,L,M)
    }
    out$c0 <- c0
    out$SS <- KLM%*%O%*%KLM
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- (NP*ns*log(2*pi) + detED + out$SS)/2
    out
}

get.Ewd <- function(w=0,format=1,ns=1,D=mclean()){
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    if (format<4) {
        J <- matrix(0,2*ns,ns)
    } else {
        J <- matrix(0,3*ns,ns)
        i3 <- (2*ns+1):(3*ns)
    }
    diag(J[i1,i1]) <- -D$dPb207U235dt      # dKdt
    diag(J[i2,i1]) <- -D$dPb206U238dt      # dLdt
    if (format>6) {
        diag(J[i3,i1]) <- -D$dPb208Th232dt # dMdt
    }
    dEdx <- w^2
    dEdx*J%*%t(J)
}

LL.ludwig.model2 <- function(ta0b0,x,anchor=0,exterr=FALSE){
    tt <- ifelse(anchor[1]==2,anchor[2],ta0b0['t'])
    nn <- length(x)
    if (x$format<4){
        a0 <- ifelse(anchor[1]==1,iratio('Pb207Pb206')[1],ta0b0['a0'])
        xy <- data2york(x,option=2)[,c('X','Y'),drop=FALSE]
        xr <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        yr <- age_to_Pb207Pb206_ratio(tt,st=0,d=x$d)[1]
        xx <- xy[,'X',drop=FALSE]
        yy <- xy[,'Y',drop=FALSE]
        dem <- deming(a=a0,b=(yr-a0)/xr,x=xx,y=yy)
        SS <- sum(dem$d^2) # = sum(d^2)
        if (exterr){
            D <- mclean(tt,d=x$d,exterr=exterr)
            dbdxr <- (a0-yr)/xr^2
            dbdyr <- 1/xr
            dxrdPbU <- -xr^2
            dyrdPbPb <- 1
            dddPbU <- dem$dddb*dbdxr*dxrdPbU
            dddPbPb <- dem$dddb*dbdyr*dyrdPbPb
            dPbUdl <- dPbPbdl <- rep(0,7)
            dPbUdl[1] <- D$dPb206U238dl38
            dPbUdl[3] <- D$dPb206U238dl34
            dPbUdl[6] <- D$dPb206U238dl30
            dPbUdl[7] <- D$dPb206U238dl26
            dPbPbdl[1] <- D$dPb207Pb206dl38
            dPbPbdl[2] <- D$dPb207Pb206dl35
            dPbPbdl[3] <- D$dPb207Pb206dl34
            dPbPbdl[5] <- D$dPb207Pb206dl31
            dPbPbdl[6] <- D$dPb207Pb206dl30
            dPbPbdl[7] <- D$dPb207Pb206dl26
            J <- dddPbU%*%dPbUdl + dddPbPb%*%dPbPbdl
            covmat <- diag(SS/(nn-2),nn,nn) + J%*%getEl()%*%t(J)
            LL <- LL.norm(as.vector(dem$d),covmat)
        } else {
            LL <- SS2LL(SS=SS,nn=nn)
        }
    } else {
        if (x$format<7){
            a0 <- ifelse(anchor[1]==1,iratio('Pb206Pb204')[1],ta0b0['a0'])
            b0 <- ifelse(anchor[1]==1,iratio('Pb207Pb204')[1],ta0b0['b0'])
            xy <- get.UPb.isochron.ratios.204(x)
        } else {
            a0 <- ifelse(anchor[1]==1,1/iratio('Pb208Pb206')[1],ta0b0['a0'])
            b0 <- ifelse(anchor[1]==1,1/iratio('Pb208Pb207')[1],ta0b0['b0'])
            xy <- get.UPb.isochron.ratios.208(x,tt=tt)[,1:4]
        }
        x6 <- xy[,1,drop=FALSE] # U238Pb206
        y6 <- xy[,2,drop=FALSE] # Pb204Pb206 or Pb208cPb206
        x7 <- xy[,3,drop=FALSE] # U235Pb207
        y7 <- xy[,4,drop=FALSE] # Pb204Pb207 or Pb208cPb207
        r86 <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        r57 <- age_to_U235Pb207_ratio(tt,st=0,d=x$d)[1]
        dem6 <- deming(a=1/a0,b=-1/(a0*r86),x=x6,y=y6)
        dem7 <- deming(a=1/b0,b=-1/(b0*r57),x=x7,y=y7)
        E <- stats::cov(cbind(dem6$d,dem7$d))
        if (exterr){
            D <- mclean(tt=tt,d=x$d,exterr=exterr)
            dbd86 <- 1/(a0*r86^2)
            dbd57 <- 1/(b0*r57^2)
            dr68d86 <- -r86^2
            dr57d75 <- -r57^2
            dd6d68 <- dem6$dddb*dbd86*dr68d86
            dd7d75 <- dem7$dddb*dbd57*dr57d75
            d68dl <- matrix(0,1,7)
            d68dl[1] <- D$dPb206U238dl38
            d68dl[3] <- D$dPb206U238dl34
            d68dl[6] <- D$dPb206U238dl30
            d68dl[7] <- D$dPb206U238dl26
            d75dl <- matrix(0,1,7)
            d75dl[2] <- D$dPb207U235dl35
            d75dl[5] <- D$dPb207U235dl31
            J <- matrix(0,2*nn,7)
            J[1:nn,] <- dd6d68%*%d68dl
            J[nn+(1:nn),] <- dd7d75%*%d75dl
            EE <- matrix(0,2*nn,2*nn)
            diag(EE)[1:nn] <- E[1,1]
            diag(EE)[nn+(1:nn)] <- E[2,2]
            diag(EE[1:nn,nn+(1:nn)]) <- E[1,2]
            diag(EE[nn+(1:nn),1:nn]) <- E[2,1]
            covmat <- EE + J %*% getEl() %*% t(J)
            dd <- c(dem6$d,dem7$d)
            LL <- LL.norm(dd,covmat)
        } else {
            dd <- cbind(dem6$d,dem7$d)
            LL <- sum(apply(dd,1,LL.norm,covmat=E))
        }
    }
    LL
}

mswd.lud <- function(ta0b0,x,anchor=0,exterr=FALSE){
    out <- list()
    ns <- length(x)
    anchored <- (anchor[1]>0)
    tanchored <- (anchor[1]==2 && length(anchor)>1 && is.numeric(anchor[2]))
    if (x$format<4){
        if (anchored) out$df <- ns-1
        else out$df <- ns-2
    } else {
        if (anchored){
            if (tanchored) out$df <- 2*ns-2
            else out$df <- 2*ns-1
        } else {
            out$df <- 2*ns-3
        }
    }
    SS <- data2ludwig(x,ta0b0w=ta0b0,anchor=anchor)$SS
    if (out$df>0){
        out$mswd <- as.vector(SS/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    } else {
        out$mswd <- 1
        out$p.value <- 1
    }
    out$n <- ns
    out
}
