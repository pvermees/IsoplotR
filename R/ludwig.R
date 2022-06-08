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
#' \item{LL}{the log likelihood of the discordia fit}
#'
#' \item{par}{a vector with the lower concordia intercept, the common
#' Pb ratios and (if \code{model=3}) the dispersion parameter}
#'
#' \item{cov}{the covariance matrix of \code{par}}
#'
#' \item{logpar}{the logarithms of \code{par}}
#'
#' \item{logcov}{the logarithms of \code{cov}}
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
ludwig.default <- function(x,exterr=FALSE,model=1,anchor=0,...){
    fit <- fit.lta0b0w(x,exterr=exterr,model=model,anchor=anchor,...)
    fit$model <- model
    fit$exterr <- exterr
    if (model==2){
        fit$logpar <- fit$par
        np <- length(fit$par)
        fit$logcov <- matrix(0,np,np)
        fit$logcov[!fit$fixed,!fit$fixed] <- solve(fit$hessian)
    } else {
        fit$LL <- fit$value
        fit$logpar <- fit$par
        fit$logcov <- fisher.lud(fit)
    }
    parnames <- c('log(t)','log(a0)','log(b0)','log(w)')
    if (model!=3) parnames <- parnames[-4]
    if (x$format < 4) parnames <- parnames[-3]
    names(fit$logpar) <-  rownames(fit$logcov) <- colnames(fit$logcov) <- parnames
    out <- exponentiate_ludwig(fit,format=x$format)
    out$n <- length(x)
    mswd <- mswd.lud(out$logpar,x=x,anchor=anchor)
    c(out,mswd)
}

fit.lta0b0w <- function(x,exterr=FALSE,model=1,anchor=0,...){
    if (measured.disequilibrium(x$d) & anchor[1]<1){
        out <- fit.lta0b0w2step(x,exterr=exterr,model=model,hessian=TRUE,...)
    } else {
        fixed <- fixit(x,anchor=anchor,model=model)
        if (model==3){
            model1fit <- fit.lta0b0w(x,model=1,anchor=anchor)
            init <- model1fit$par
            LL0 <- LL.lud.w(-5,lta0b0=init,x=x)
            LLw <- LL.lud.w(-4,lta0b0=init,x=x)
            if (LLw<LL0){ # zero dispersion
                out <- model1fit
                out$model <- model
                out$par <- c(out$par,-Inf)
                out$fixed[length(fixed)] <- TRUE
                return(out)
            } else {
                w <- stats::optimise(LL.lud.w,interval=init[1]+c(-10,0),
                                     lta0b0=init,x=x,maximum=TRUE)$maximum
            }
            init <- c(init,w)
        } else if (anchor[1] %in% c(1,2)){
            init <- anchored.lta0b0.init(x=x,anchor=anchor)
        } else {
            init <- get.lta0b0.init(x)
        }
        lower <- (init-1)[!fixed]
        upper <- (init+1)[!fixed]
        if (model==2){
            out <- optifix(parms=init,fn=LL.lud.model2,method="L-BFGS-B",
                           x=x,fixed=fixed,lower=lower,upper=upper,
                           exterr=exterr,hessian=TRUE)
        } else {
            out <- optifix(parms=init,fn=LL.lud,#gr=LL.lud.gr,
                           method="L-BFGS-B",x=x,exterr=exterr,fixed=fixed,
                           lower=lower,upper=upper,control=list(fnscale=-1),...)
        }
        out$x <- x
        out$model <- model
        out$exterr <- exterr
        out$fixed <- fixed
    }
    out    
}

# special case for measured disequilibrium
fit.lta0b0w2step <- function(x,exterr=FALSE,model=1,...){
    # 1. fit without disequilibrium
    X <- x
    X$d <- diseq()
    fit1 <- fit.lta0b0w(X,exterr=exterr,model=model)
    # 2. check that the measured disequilibrium is physically possible
    X$d <- replace.impossible.diseq(tt=exp(fit1$par[1]),d=x$d)
    # 3. model-2 fit
    init <- fit1$par
    np2 <- ifelse(x$format<4,2,3)
    fit2 <- stats::optim(init[1:np2],fn=LL.lud.model2,method="L-BFGS-B",
                         x=X,lower=init-1,upper=init+1,exterr=exterr,...)
    if (model==2){
        out <- fit2
        out$fixed <- rep(FALSE,np2)
    } else {
        init[1:np2] <- fit2$par
        out <- optifix(parms=init,fn=LL.lud,#gr=LL.lud.gr,
                       method="L-BFGS-B",x=X,exterr=exterr,fixed=fit1$fixed,
                       lower=init-1,upper=init+1,control=list(fnscale=-1),...)
    }
    out$x <- X
    out$exterr <- exterr
    out$model <- model
    out
}

get.lta0b0.init <- function(x){
    out <- list()
    xy <- data2york(x,option=2)
    # First, get the approximate intercept age using a model-2 regression
    fit <- stats::lm(xy[,'Y'] ~ xy[,'X'])
    a <- fit$coef[1]
    b <- fit$coef[2]
    covmat <- stats::vcov(fit)
    # no disequilibrium correction:
    tint <- concordia.intersection.ab(a,b,covmat=covmat)
    tt <- tint['t[l]']
    # Then, estimate the common Pb intercept(s)
    minval <- 0.01
    if (x$format<4){
        a0 <- a
        expinit <- c(tt,a0)
        labels <- c('log(t)','log(a0)')
    } else {
        if (x$format<7) option <- 3
        else option <- 6
        xy <- data2york(x,option=option,tt=tt) # 04-08c/06 vs 38/06
        fit <- stats::lm(xy[,'Y'] ~ xy[,'X'])
        a0 <- 1/fit$coef[1]
        xy <- data2york(x,option=option+1,tt=tt) # 04-08c/07 vs 38/07
        fit <- stats::lm(xy[,'Y'] ~ xy[,'X'])
        b0 <- 1/fit$coef[1]
        expinit <- c(tt,a0,b0)
        labels <- c('log(t)','log(a0)','log(b0)')
    }
    expinit[expinit<=0] <- 1e-5
    out <- log(expinit)
    names(out) <- labels
    out
}

anchored.lta0b0.init <- function(x,anchor=1){
    if (x$format<4) np <- 2
    else np <- 3
    init <- rep(0,np)
    names(init) <- c('log(t)','log(a0)','log(b0)')[1:np]
    if (anchor[1]==1){ # fix common Pb composition
        xy <- data2york(x,option=2)
        if (x$format%in%c(1,2,3)){
            i76 <- iratio('Pb207Pb206')[1]
            init['log(a0)'] <- log(i76)
        } else if (x$format%in%c(4,5,6)){
            i64 <- iratio('Pb206Pb204')[1]
            i74 <- iratio('Pb207Pb204')[1]
            init['log(a0)'] <- log(i64)
            init['log(b0)'] <- log(i74)
            i76 <- i74/i64
        } else if (x$format%in%c(7,8)){
            i86 <- iratio('Pb208Pb206')[1]
            i87 <- iratio('Pb208Pb207')[1]
            init['log(a0)'] <- -log(i86)
            init['log(b0)'] <- -log(i87)
            i76 <- i86/i87
        } else {
            stop('incorrect input format')
        }
        fit <- stats::lm(I(xy[,'Y']-i76) ~ 0 + xy[,'X'])
        b <- fit$coef
        covmat <- matrix(0,2,2)
        covmat[2,2] <- stats::vcov(fit)
        tint <- concordia.intersection.ab(i76,b,covmat=covmat,
                                          d=mediand(x$d))[1]
        init['log(t)'] <- log(tint)
    } else if (anchor[1]==2){ # fix age
        init['log(t)'] <- log(anchor[2])
        if (x$format<4){
            xy <- data2york(x,option=2)
            TW <- age_to_terawasserburg_ratios(anchor[2],st=0,
                                               exterr=FALSE,d=mediand(x$d))
            b <- stats::lm(I(xy[,'Y']-TW$x['Pb207Pb206']) ~
                               0 + I(xy[,'X']-TW$x['U238Pb206']))$coef
            init['log(a0)'] <- log(TW$x['Pb207Pb206'] - b*TW$x['U238Pb206'])
        } else {
            r86 <- age_to_U238Pb206_ratio(anchor[2],st=0,d=mediand(x$d))[1]
            if (x$format<7){
                xy6 <- data2york(x,option=3)
                xy7 <- data2york(x,option=4)
            } else {
                xy6 <- data2york(x,option=6,tt=anchor[2])
                xy7 <- data2york(x,option=7,tt=anchor[2])
            }
            b <- stats::lm(xy6[,'Y'] ~ 0 + I(xy6[,'X']-r86))$coef
            init['log(a0)'] <- -(log(-b)+log(r86))
            r57 <- age_to_U235Pb207_ratio(anchor[2],st=0,d=mediand(x$d))[1]
            b <- stats::lm(xy7[,'Y'] ~ 0 + I(xy7[,'X']-r57))$coef
            init['log(b0)'] <- -(log(-b)+log(r57))
        }
    } else { 
        stop("Invalid discordia regression anchor.")
    }
    init
}

exponentiate_ludwig <- function(fit,format){
    out <- fit
    np <- length(fit$logpar)
    J <- matrix(0,np,np)
    out$par <- exp(fit$logpar)
    parnames <- c('t','a0','b0','w')
    if (fit$model!=3) parnames <- parnames[-4]
    if (format < 4) parnames <- parnames[-3]
    names(out$par) <- parnames
    rownames(J) <- parnames
    diag(J) <- exp(fit$logpar[1:np])
    out$cov <- J %*% fit$logcov %*% t(J)
    out
}

mswd.lud <- function(lta0b0,x,anchor=0){
    if (measured.disequilibrium(x$d))
        x$d <- replace.impossible.diseq(tt=exp(lta0b0[1]),d=x$d)
    ns <- length(x)
    out <- list()
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
    SS <- data2ludwig(x,lta0b0w=lta0b0)$SS
    if (out$df>0){
        out$mswd <- as.vector(SS/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    } else {
        out$mswd <- 1
        out$p.value <- 1
    }
    out
}

replace.impossible.diseq <- function(tt,d){
    out <- d
    D <- mclean(tt=tt,d=d)
    if (D$U48i<0){
        out$U48$x <- 0
        out$U48$option <- 1
    }
    if (D$ThUi<0){
        out$ThU$x <- 0
        out$ThU$option <- 1
    }
    out
}

LL.lud.model2 <- function(lta0b0,x,exterr=FALSE){
    tt <- exp(lta0b0[1])
    a0 <- exp(lta0b0[2])
    nn <- length(x)
    if (x$format<4){
        xy <- data2york(x,option=2)[,c('X','Y'),drop=FALSE]
        xr <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        yr <- age_to_Pb207Pb206_ratio(tt,st=0,d=x$d)[1]
        yp <- a0+(yr-a0)*xy[,'X']/xr
        dy <- yp-xy[,'Y']
        SS <- sum(dy^2)
        if (exterr){
            D <- mclean(tt,d=x$d,exterr=exterr)
            dypdxr <- (a0-yr)*xy[,'X',drop=FALSE]/xr^2
            dypdyr <- xy[,'X',drop=FALSE]/xr
            dxrdPbU <- -xr^2
            dyrdPbPb <- 1
            dPbUdl <- rep(0,7)
            dPbPbdl <- rep(0,7)
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
            J <- (dypdxr*dxrdPbU)%*%dPbUdl + (dypdyr*dyrdPbPb)%*%dPbPbdl
            covmat <- diag(SS/(nn-2),nn,nn) + J%*%getEl()%*%t(J)
            LL <- LL.norm(dy,covmat)
        } else {
            LL <- SS2LL(SS=SS,nn=nn)
        }
    } else {
        b0 <- exp(lta0b0[3])
        if (x$format<7) xy <- get.UPb.isochron.ratios.204(x)
        else xy <- get.UPb.isochron.ratios.208(x,tt=tt)[,1:4]
        x6 <- xy[,1,drop=FALSE] # U238Pb206
        y6 <- xy[,2,drop=FALSE] # Pb204Pb206 or Pb208cPb206
        x7 <- xy[,3,drop=FALSE] # U235Pb207
        y7 <- xy[,4,drop=FALSE] # Pb204Pb207 or Pb208cPb207
        r86 <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        r57 <- age_to_U235Pb207_ratio(tt,st=0,d=x$d)[1]
        y6p <- (r86-x6)/(a0*r86)
        y7p <- (r57-x7)/(b0*r57)
        dy <- cbind(y6-y6p,y7-y7p)
        E <- stats::cov(dy)
        if (exterr){
            D <- mclean(tt=tt,d=x$d,exterr=exterr)
            dy6pd68 <- -x6/a0
            dy7pd75 <- -x7/b0
            d68dl <- matrix(0,1,7)
            d68dl[1] <- D$dPb206U238dl38
            d68dl[3] <- D$dPb206U238dl34
            d68dl[6] <- D$dPb206U238dl30
            d68dl[7] <- D$dPb206U238dl26
            d75dl <- matrix(0,1,7)
            d75dl[2] <- D$dPb207U235dl35
            d75dl[5] <- D$dPb207U235dl31
            J <- matrix(0,2*nn,7)
            J[1:nn,] <- dy6pd68%*%d68dl
            J[nn+(1:nn),] <- dy7pd75%*%d75dl
            EE <- matrix(0,2*nn,2*nn)
            diag(EE)[1:nn] <- E[1,1]
            diag(EE)[nn+(1:nn)] <- E[2,2]
            diag(EE[1:nn,nn+(1:nn)]) <- E[1,2]
            diag(EE[nn+(1:nn),1:nn]) <- E[2,1]
            covmat <- EE + J %*% getEl() %*% t(J)
            LL <- LL.norm(c(y6-y6p,y7-y7p),covmat)
        } else {
            LL <- sum(apply(dy,1,LL.norm,covmat=E))
        }
    }
    LL
}

LL.lud <- function(lta0b0w,x,exterr=FALSE,LL=TRUE){
    l <- data2ludwig(x,lta0b0w=lta0b0w,exterr=exterr)
    if (LL) out <- l$LL
    else out <- l$SS
    out
}
LL.lud.gr <- function(lta0b0w,x,exterr=FALSE){
    data2ludwig(x,lta0b0w=lta0b0w,exterr=exterr,jacobian=TRUE)$jacobian
}
# LL to initialise w:
LL.lud.w <- function(w,lta0b0,x){
    lta0b0w <- c(lta0b0,w)
    LL.lud(lta0b0w=lta0b0w,x=x,LL=TRUE)
}

fisher.lud <- function(fit){
    fish <- data2ludwig(fit$x,lta0b0w=fit$par,
                        exterr=fit$exterr,hessian=TRUE)$hessian
    ns <- length(fit$x)
    np <- length(fit$par)
    i1 <- 1:ns
    i2 <- ns+(1:np)
    out <- matrix(0,ns+np,ns+np)
    anchorfish(AA=fish[i1,i1],BB=fish[i1,i2],
               CC=fish[i2,i1],DD=fish[i2,i2],fixed=fit$fixed)
}
anchorfish <- function(AA,BB,CC,DD,fixed=rep(FALSE,nrow(DD))){
    ns <- nrow(AA)
    np <- nrow(DD)
    i <- (1:np)[!fixed]
    bb <- BB[,i,drop=FALSE]
    cc <- CC[i,,drop=FALSE]
    dd <- DD[i,i,drop=FALSE]
    out <- matrix(0,np,np)
    out[i,i] <- blockinverse(AA,bb,cc,dd)
    out
}

data2ludwig <- function(x,lta0b0w,exterr=FALSE,jacobian=FALSE,hessian=FALSE){
    out <- list()
    U <- iratio('U238U235')[1]
    lt <- lta0b0w[1]
    a0 <- lta0b0w[2]
    ns <- length(x)
    zeros <- rep(0,ns)
    X <- zeros
    Y <- zeros
    K0 <- zeros
    D <- mclean(tt=exp(lt),d=x$d,exterr=exterr)
    if (x$format%in%c(1,2,3)){
        NP <- 2 # number of fit parameters (lt, a0)
        NR <- 2 # number of isotopic ratios (X, Y)
    } else if (x$format%in%c(4,5,6)){
        b0 <- lta0b0w[3]
        Z <- zeros
        L0 <- zeros
        NP <- 3 # lt, a0, b0
        NR <- 3 # X, Y, Z
    } else if (x$format%in%c(7,8)){
        b0 <- lta0b0w[3]
        Z <- zeros
        W <- zeros
        L0 <- zeros
        NP <- 3 # lt, a0, b0
        NR <- 4 # X, Y, Z, W
    } else {
        stop('Incorrect input format.')
    }
    np <- min(NP+1,length(lta0b0w))  # np = NP (no w) or NP+1 (with w)
    if (np==(NP+1)) w <- lta0b0w[np] # model 3
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
    if (np==(NP+1)){ # fit overdispersion
        Ew <- get.Ew(w=w,format=x$format,ns=ns,D=D)
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
        K0 <- X - D$Pb207U235 + exp(a0)*U*(D$Pb206U238 - Y)
        A <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*exp(a0)*U +
               K0%*%(O[i1,i2]+t(O[i2,i1])))
        B <- -(exp(a0)*U*(O[i1,i1]+t(O[i1,i1]))*exp(a0)*U +
               (O[i2,i2]+t(O[i2,i2])) +
               exp(a0)*U*(O[i1,i2]+t(O[i1,i2])) +
               (O[i2,i1]+t(O[i2,i1]))*exp(a0)*U)
        L <- as.vector(solve(B,A))
        c0 <- Y - D$Pb206U238 - L
        K <- X - D$Pb207U235 - exp(a0)*U*c0
        KLM <- c(K,L)
    } else if (x$format%in%c(4,5,6)){
        K0 <- X - D$Pb207U235 - U*exp(b0)*Z
        L0 <- Y - D$Pb206U238 - exp(a0)*Z
        V <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*U*exp(b0) +
               L0%*%(O[i1,i2]+t(O[i2,i1]))*U*exp(b0) +
               K0%*%(O[i1,i2]+t(O[i2,i1]))*exp(a0) +
               L0%*%(O[i2,i2]+t(O[i2,i2]))*exp(a0) +
               K0%*%(O[i1,i3]+t(O[i3,i1])) +
               L0%*%(O[i2,i3]+t(O[i3,i2])))
        W <- -(U*exp(b0)*(O[i1,i1]+t(O[i1,i1]))*U*exp(b0) +
               U*exp(b0)*(O[i1,i2]+t(O[i1,i2]))*exp(a0) +
               U*exp(b0)*(O[i1,i3]+t(O[i1,i3])) +
               exp(a0)*(O[i2,i1]+t(O[i2,i1]))*U*exp(b0) +
               exp(a0)*(O[i2,i2]+t(O[i2,i2]))*exp(a0) +
               exp(a0)*(O[i2,i3]+t(O[i2,i3])) +
               (O[i3,i1]+t(O[i3,i1]))*U*exp(b0) +
               (O[i3,i2]+t(O[i3,i2]))*exp(a0) +
               (O[i3,i3]+t(O[i3,i3])))
        M <- as.vector(solve(W,V))
        c0 <- as.vector(Z - M)
        K <- as.vector(X - D$Pb207U235 - U*exp(b0)*c0)
        L <- as.vector(Y - D$Pb206U238 - exp(a0)*c0)
        KLM <- c(K,L,M)
    } else if (x$format%in%c(7,8)){
        Wd <- diag(W)
        K0 <- X - D$Pb207U235 - (Z-D$Pb208Th232)*U*W*exp(b0)
        L0 <- Y - D$Pb206U238 - (Z-D$Pb208Th232)*W*exp(a0)
        AA <- (Wd%*%O[i1,i1]%*%Wd)*(U*exp(b0))^2 +
            (Wd%*%O[i2,i2]%*%Wd)*exp(a0)^2 + O[i3,i3] +
            U*exp(a0)*exp(b0)*Wd%*%(O[i1,i2]+O[i2,i1])%*%Wd +
            U*exp(b0)*(Wd%*%O[i1,i3]+O[i3,i1]%*%Wd) +
            exp(a0)*(Wd%*%O[i2,i3]+O[i3,i2]%*%Wd)
        BT <- t(U*exp(b0)*K0%*%O[i1,i1]%*%Wd +
                exp(a0)*L0%*%O[i2,i2]%*%Wd +
                exp(a0)*K0%*%O[i1,i2]%*%Wd +
                U*exp(b0)*L0%*%O[i2,i1]%*%Wd +
                K0%*%O[i1,i3] + L0%*%O[i2,i3])
        CC <- U*exp(b0)*Wd%*%O[i1,i1]%*%K0 +
            exp(a0)*Wd%*%O[i2,i2]%*%L0 +
            exp(a0)*Wd%*%O[i2,i1]%*%K0 +
            U*exp(b0)*Wd%*%O[i1,i2]%*%L0 +
            O[i3,i1]%*%K0 + O[i3,i2]%*%L0
        M <- as.vector(solve(-(AA+t(AA)),(BT+CC)))
        c0 <- as.vector(Z - D$Pb208Th232 - M)
        K <- as.vector(X - D$Pb207U235 - c0*exp(b0)*U*W)
        L <- as.vector(Y - D$Pb206U238 - c0*exp(a0)*W)
        KLM <- c(K,L,M)
    }
    out$c0 <- c0
    out$SS <- KLM%*%O%*%KLM
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- -(NP*ns*log(2*pi) + detED + out$SS)/2
    if (jacobian | hessian){
        JKLM <- matrix(0,NP*ns,ns+NP)
        out$jacobian <- rep(0,np) # derivatives of LL w.r.t lt, a0, (b0) (and w)
        if (x$format<4){
            colnames(JKLM) <- c(paste0('c0[',i1,']'),'lt','a0')
            rownames(JKLM) <- c(paste0('K[',i1,']'),paste0('L[',i1,']'))
            names(out$jacobian) <- c('lt','a0','w')[1:np]
        } else {
            colnames(JKLM) <- c(paste0('c0[',i1,']'),'lt','a0','b0')
            rownames(JKLM) <- c(paste0('K[',i1,']'),
                                paste0('L[',i1,']'),
                                paste0('M[',i1,']'))
            names(out$jacobian) <- c('lt','a0','b0','w')[1:np]
        }
        dtdlt <- exp(lt)
        dKdt <- -D$dPb207U235dt
        dKda0 <- 0; dKdb0 <- 0; dKdc0 <- 0
        dLdt <- -D$dPb206U238dt
        dLda0 <- 0; dLdb0 <- 0; dLdc0 <- 0
        dMdt <- 0; dMda0 <- 0; dMdb0 <- 0; dMdc0 <- 0
        if (x$format%in%c(1,2,3)){
            dKda0 <- -U*c0*exp(a0)
            dKdc0 <- -U*exp(a0)
            dLdc0 <- -1
        } else if (x$format%in%c(4,5,6)){
            dKdb0 <- -U*c0*exp(b0)
            dKdc0 <- -U*exp(b0)
            dLda0 <- -c0*exp(a0)
            dLdc0 <- -exp(a0)
            dMdc0 <- -1
        } else if (x$format%in%c(7,8)){
            dKdb0 <- -c0*U*W*exp(b0)
            dKdc0 <- -U*W*exp(b0)
            dLda0 <- -c0*W*exp(a0)
            dLdc0 <- -W*exp(a0)
            dMdt <- -D$dPb208Th232dt
            dMdc0 <- -1
        }
        diag(JKLM[i1,i1]) <- dKdc0
        JKLM[i1,'lt']  <- dKdt*dtdlt
        JKLM[i1,'a0'] <- dKda0
        diag(JKLM[i2,i1]) <- dLdc0
        JKLM[i2,'lt']  <- dLdt*dtdlt
        JKLM[i2,'a0'] <- dLda0
        if (x$format>3){
            JKLM[i1,'b0'] <- dKdb0
            JKLM[i2,'b0'] <- dLdb0
            diag(JKLM[i3,i1]) <- dMdc0
            JKLM[i3,'lt']  <- dMdt*dtdlt
            JKLM[i3,'a0']  <- dMda0
            JKLM[i3,'b0']  <- dMdb0
            out$jacobian['b0'] <- -KLM%*%O%*%JKLM[,'b0']
        }
        out$jacobian['lt'] <- -KLM%*%O%*%JKLM[,'lt']
        out$jacobian['a0'] <- -KLM%*%O%*%JKLM[,'a0']
        if (np==(NP+1)){
            dEDdw <- get.Ew(w=w,format=x$format,ns=ns,D=D,deriv=1)
            dlnDetEDdw <- trace(O%*%dEDdw)
            dOdw <- -O%*%dEDdw%*%O
            out$jacobian['w'] <- -(dlnDetEDdw + KLM%*%dOdw%*%KLM)/2
        }
    }
    if (hessian){
        d2tdlt2 <- exp(lt)
        d2Kdt2 <- -D$d2Pb207U235dt2
        d2Kdlt2 <- d2Kdt2*dtdlt^2 + dKdt*d2tdlt2
        d2Kda02 <- 0; d2Kdb02 <- 0; d2Kdc02 <- 0
        d2Kdltda0 <- 0; d2Kdltdb0 <- 0; d2Kdltdc0 <- 0
        d2Kda0db0 <- 0; d2Kda0dc0 <- 0; d2Kdb0dc0 <- 0
        d2Ldt2 <- -D$d2Pb206U238dt2
        d2Ldlt2 <- d2Ldt2*dtdlt^2 + dLdt*d2tdlt2
        d2Lda02 <- 0; d2Ldb02 <- 0; d2Ldc02 <- 0
        d2Ldltda0 <- 0; d2Ldltdb0 <- 0; d2Ldltdc0 <- 0
        d2Lda0db0 <- 0; d2Lda0dc0 <- 0; d2Ldb0dc0 <- 0
        if (x$format%in%c(1,2,3)){
            d2Kda02 <- -U*c0*exp(a0)
            d2Kda0dc0 <- -U*exp(a0)
        } else if (x$format%in%c(4,5,6)){
            d2Kdb02 <- -U*c0*exp(b0)
            d2Kdb0dc0 <- -U*exp(b0)
            d2Lda02 <- -c0*exp(a0)
            d2Lda0dc0 <- -exp(a0)
            d2Mdlt2 <- 0
        } else if (x$format%in%c(7,8)){
            d2Kdb02 <- -c0*U*W*exp(b0)
            d2Kdb0dc0 <- -U*W*exp(b0)
            d2Lda02 <- -c0*W*exp(a0)
            d2Lda0dc0 <- -W*exp(a0)
            d2Mdt2 <- -D$d2Pb208Th232dt2
            d2Mdlt2 <- d2Mdt2*dtdlt^2 + dMdt*d2tdlt2
        }
        ZEROS <- rep(0,NP*ns)
        d2KLMdlt2 <- ZEROS; d2KLMda02 <- ZEROS; d2KLMdb02 <- ZEROS
        d2KLMda0dc0 <- matrix(0,NP*ns,ns); d2KLMdb0dc0 <- matrix(0,NP*ns,ns)
        d2KLMdlt2[i1] <- d2Kdlt2
        d2KLMdlt2[i2] <- d2Ldlt2
        d2KLMda02[i1] <- d2Kda02
        diag(d2KLMda0dc0[i1,i1]) <- d2Kda0dc0
        diag(d2KLMda0dc0[i2,i1]) <- d2Lda0dc0
        out$hessian <- matrix(0,ns+np,ns+np)
        if (x$format<4){
            hesnames <- c(paste0('c0[',1:ns,']'),'lt','a0','w')[1:(ns+np)]
        } else {
            hesnames <- c(paste0('c0[',1:ns,']'),'lt','a0','b0','w')[1:(ns+np)]
            d2KLMdb02[i1] <- d2Kdb02
            diag(d2KLMdb0dc0[i1,i1]) <- d2Kdb0dc0
            d2KLMdlt2[i3] <- d2Mdlt2
        }
        rownames(out$hessian) <- hesnames
        colnames(out$hessian) <- hesnames
        out$hessian[i1,i1] <- t(JKLM[,i1])%*%O%*%JKLM[,i1]                # d2dc02
        out$hessian['lt','lt'] <-
            t(JKLM[,'lt'])%*%O%*%JKLM[,'lt'] + KLM%*%O%*%d2KLMdlt2        # d2dt2
        out$hessian['a0','a0'] <-
            t(JKLM[,'a0'])%*%O%*%JKLM[,'a0'] + KLM%*%O%*%d2KLMda02        # d2da02
        out$hessian['lt','a0'] <- t(JKLM[,'lt'])%*%O%*%JKLM[,'a0']        # d2dltda0
        out$hessian['lt',i1] <- t(JKLM[,'lt'])%*%O%*%JKLM[,i1]            # d2dc0dt
        out$hessian['a0',i1] <-
            t(JKLM[,'a0'])%*%O%*%JKLM[,i1] + KLM%*%O%*%d2KLMda0dc0        # d2dc0da0
        out$hessian['a0','lt'] <- out$hessian['lt','a0']                  # d2da0dlt
        out$hessian[i1,'lt'] <- out$hessian['lt',i1]                      # d2dltdc0
        out$hessian[i1,'a0'] <- out$hessian['a0',i1]                      # d2da0dc0
        if (x$format>3){ # with b0
            out$hessian['lt','b0'] <- t(JKLM[,'lt'])%*%O%*%JKLM[,'b0']    # d2dltdb0
            out$hessian['a0','b0'] <- t(JKLM[,'a0'])%*%O%*%JKLM[,'b0']    # d2da0db0
            out$hessian['b0','b0'] <-
                t(JKLM[,'b0'])%*%O%*%JKLM[,'b0'] + KLM%*%O%*%d2KLMdb02    # d2db02
            out$hessian['b0',i1] <-
                t(JKLM[,'b0'])%*%O%*%JKLM[,i1] + KLM%*%O%*%d2KLMdb0dc0    # d2dc0db0
            out$hessian['b0','lt'] <- out$hessian['lt','b0']              # d2db0dlt
            out$hessian['b0','a0'] <- out$hessian['a0','b0']              # d2db0da0
            out$hessian[i1,'b0'] <- out$hessian['b0',i1]                  # d2db0dc0
        }
        if (np==(NP+1)){ # with w
            d2EDdw2 <- get.Ew(w=w,format=x$format,ns=ns,D=D,deriv=2)
            d2Odw2 <- -(dOdw%*%dEDdw%*%O + O%*%d2EDdw2%*%O + O%*%dEDdw%*%dOdw)
            d2lnDetEDdw2 <- trace(dOdw%*%dEDdw + O%*%d2EDdw2)
            out$hessian['w','w'] <- (d2lnDetEDdw2 + KLM%*%d2Odw2%*%KLM)/2 # d2dw2
            out$hessian['w',i1] <- KLM%*%dOdw%*%JKLM[,i1]                 # d2dwdc0
            out$hessian['w','lt'] <- KLM%*%dOdw%*%JKLM[,'lt']             # d2dwdt
            out$hessian['w','a0'] <- KLM%*%dOdw%*%JKLM[,'a0']             # d2dwda0
            out$hessian[i1,'w'] <- out$hessian['w',i1]                    # d2dc0dw
            out$hessian['lt','w'] <- out$hessian['w','lt']                # d2dtdw
            out$hessian['a0','w'] <- out$hessian['w','a0']                # d2da0dw
            if (x$format>3){
                out$hessian['w','b0'] <- KLM%*%dOdw%*%JKLM[,'b0']         # d2dwdb0
                out$hessian['b0','w'] <- out$hessian['w','b0']            # d2db0dw
            }
        }
        if (any(eigen(out$hessian)$values<0)){
            warning('Non positive definite Hessian matrix in Ludwig regression.')
            out$hessian <- nearPD(out$hessian)
        }
    }
    out
}

get.Ew <- function(w=0,format=1,ns=1,D=mclean(),deriv=0){
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
    if (deriv==1) {
        dEdx <- 2*exp(w)^2
    } else if (deriv==2) {
        dEdx <- 4*exp(w)^2
    } else {
        dEdx <- exp(w)^2
    }
    dEdx*J%*%t(J)
}

fixit <- function(x,anchor=0,model=1,joint=TRUE){
    if (x$format>3 & joint) NP <- 3
    else NP <- 2
    if (model==3) np <- NP+1
    else np <- NP
    out <- rep(FALSE,np)
    if (anchor[1]>0){ # anchor t or a0(,b0)
        if (anchor[1]==2) out[1] <- TRUE # fix t
        else out[2:NP] <- TRUE # fix a0(,b0)
    }
    out
}
