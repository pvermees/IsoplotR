#' @title
#' Linear regression of U-Pb data with correlated errors, taking
#' into account decay constant uncertainties.
#'
#' @description
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
#' ages.
#'
#' @param x an object of class \code{UPb}
#' 
#' @param alpha cutoff value for confidence intervals
#' 
#' @param ... optional
#' 
#' @return
#' \describe{
#'
#' \item{par}{a vector with the lower concordia intercept,
#' the common Pb ratios and the dispersion parameter.}
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
#' uncertainty (i.e. decay constants)?
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
#' The first element is a boolean flag indicating whether the
#' discordia line should be anchored. If this is \code{FALSE}, then
#' the second item is ignored and both the common Pb composition and
#' age are estimated.
#'
#' If the first element is \code{TRUE} and the second element is
#' \code{NA}, then the common Pb composition is fixed at the values
#' stored in \code{settings('iratio',...)}.
#'
#' If the first element is \code{TRUE} and the second element is a
#' number, then the discordia line is forced to intersect the
#' concordia line at an age equal to that number.
#'
#' @seealso \code{\link{concordia}}, \code{\link{titterington}},
#'     \code{\link{isochron}}
#' @rdname ludwig
#' @export
ludwig.UPb <- function(x,exterr=FALSE,alpha=0.05,model=1,
                       anchor=list(FALSE,NA)){
    if (x$format %in% c(1,2,3)) parnames <- c('t','76i')
    else if (x$format %in% c(4,5,6)) parnames <- c('t','64i','74i')
    else if (x$format %in% c(7,8)) parnames <- c('t','68i','78i')
    else stop("Illegal input format.")
    if (model==3) parnames <- c(parnames,'w')
    out <- get.lta0b0w(x,exterr=exterr,model=model,anchor=anchor)
    names(out$par) <- parnames
    rownames(out$cov) <- parnames
    colnames(out$cov) <- parnames
    out$n <- length(x)
    mswd <- mswd.lud(out$par,x=x,anchor=anchor)
    c(out,mswd)
}

mswd.lud <- function(lta0b0,x,anchor=list(FALSE,NA)){
    ns <- length(x)
    out <- list()
    anchored <- anchor[[1]]
    tanchored <- is.numeric(anchor[[2]])
    if (x$format<4){
        if (anchored) out$df <- ns-1
        else out$df <- ns-2
    } else {
        if (anchored && tanchored) out$df <- 2*ns-1
        else if (tanchored) out$df <- 2*ns-2
        else out$df <- 2*ns-3
    }
    SS <- data2ludwig(x,lta0b0=lta0b0)$SS
    if (out$df>0){
        out$mswd <- as.vector(SS/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    } else {
        out$mswd <- 1
        out$p.value <- 1
    }
    out
}

get.lta0b0w <- function(x,exterr=FALSE,model=1,
                        anchor=list(FALSE,NA),w=NA,...){
    out <- list()
    fit2 <- get.lta0b0.model2(x,anchor=anchor)
    lta0b0 <- fit2$lta0b0
    if (model==2){
        out$par <- lta0b0
        out$cov <- fit2$cov
    } else {
        init <- lta0b0
        lower <- init-5
        upper <- init+5
        if (model==3){
            if (is.na(w)) ww <- init[1]-5
            init <- c(init,ww)
            lower <- c(lower,ww-10)
            upper <- c(upper,upper[1])
        }
        fit <- optifix(parms=init,fn=LL.lud.UPb,gr=LL.lud.UPb.gr,
                       method="L-BFGS-B",x=x,exterr=exterr,
                       fixed=fixit(x,anchor=anchor,model=model,w=w),
                       lower=lower,upper=upper,
                       control=list(fnscale=-1),...)
        out <- list()
        out$LL <- fit$value
        out$par <- fit$par
        out$cov <- fisher.lud(x,fit=fit,exterr=exterr,anchor=anchor)
    }
    out$model <- model
    out$exterr <- exterr
    out
}

get.lta0b0.model2 <- function(x,anchor=list(FALSE,NA),...){
    xy <- data2york(x,option=2)[,c('X','Y'),drop=FALSE]
    fit <- lm(xy[,'Y'] ~ xy[,'X'])
    x0 <- -fit$coef[1]/fit$coef[2]
    tt <- get.Pb206U238.age(1/x0)[1]
    if (x$format%in%c(1,2,3)){
        np <- 2
        a0 <- fit$coef[1]
        init <- log(c(tt,a0))
    } else if (x$format%in%c(4,5,6)){
        np <- 3
        xyz <- get.UPb.isochron.ratios.204(x)
        fit <- lm(xyz[,'Pb204Pb206'] ~ xyz[,'U238Pb206'])
        a0 <- 1/fit$coef[1]
        fit <- lm(xyz[,'Pb204Pb207'] ~ xyz[,'U235Pb207'])
        b0 <- 1/fit$coef[1]
        init <- log(c(tt,a0,b0))
    } else if (x$format%in%c(7,8)){
        np <- 3
        xyz <- get.UPb.isochron.ratios.208(x,tt=tt)
        fit6 <- lm(xyz[,'Pb208cPb206'] ~ xyz[,'U238Pb206'])
        a0 <- 1/fit6$coef[1]
        fit7 <- lm(xyz[,'Pb208cPb207'] ~ xyz[,'U235Pb207'])
        b0 <- 1/fit7$coef[1]
        init <- log(c(tt,a0,b0))
    } else {
        stop('Incorrect input format.')
    }
    fit <- optifix(parms=init,fn=SS.model2,method="L-BFGS-B",
                   x=x,fixed=fixit(x,anchor=anchor,model=2),
                   lower=init-5,upper=init+5,hessian=TRUE,...)
    out <- list()
    out$lta0b0 <- fit$par
    out$cov <- np*fit$value/(length(x)-2)*solve(fit$hessian) # from R-intro
    out
}

SS.model2 <- function(lta0b0,x){
    tt <- exp(lta0b0[1])
    a0 <- exp(lta0b0[2])
    out <- list()
    if (x$format<4){
        xy <- data2york(x,option=2)[,c('X','Y'),drop=FALSE]
        xr <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        yr <- age_to_Pb207Pb206_ratio(tt,st=0,d=x$d)[1]
        yp <- yr+(exp(a0)-yr)*(xr-xy[,'X'])/xr
        out <- sum((yp-xy[,'Y'])^2)
    } else {
        b0 <- exp(lta0b0[3])
        ns <- length(x)
        if (x$format<7){
            xy <- get.UPb.isochron.ratios.204(x)
        } else {
            xy <- get.UPb.isochron.ratios.208(x,tt=tt)[,1:4]
        }
        x6 <- xy[,1] # U238Pb206
        y6 <- xy[,2] # Pb204Pb206 or Pb208cPb206
        x7 <- xy[,3] # U235Pb207
        y7 <- xy[,4] # Pb204Pb207 or Pb208cPb207
        r86 <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        r57 <- age_to_U235Pb207_ratio(tt,st=0,d=x$d)[1]
        y6p <- (r86-x6)/(r86*exp(a0))
        SS6 <- sum((y6p-y6)^2)
        y7p <- (r57-x7)/(r57*exp(b0))
        SS7 <- sum((y7p-y7)^2)
        out <- SS6 + SS7
    }
    out
}

LL.lud.UPb <- function(lta0b0w,x,exterr=FALSE,LL=TRUE){
    l <- data2ludwig(x,lta0b0w=lta0b0w,exterr=exterr)
    if (LL) out <- l$LL
    else out <- l$SS
    out
}

LL.lud.UPb.gr <- function(lta0b0w,x,exterr=FALSE){
    data2ludwig(x,lta0b0w=lta0b0w,exterr=exterr,jacobian=TRUE)$jacobian
}

fisher.lud <- function(x,fit,exterr=TRUE,anchor=list(FALSE,NA),...){
    fish <- data2ludwig(x,lta0w=fit$par,exterr=exterr,hessian=TRUE)$hessian
    ns <- length(x)
    np <- length(fit$par)
    i1 <- 1:ns
    i2 <- (ns+1):(ns+np)
    anchorfish(AA=fish[i1,i1],BB=fish[i1,i2],
               CC=fish[i2,i1],DD=fish[i2,i2],anchor=anchor)
}
anchorfish <- function(AA,BB,CC,DD,anchor=list(FALSE,NA)){
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

data2ludwig <- function(x,lta0b0w,exterr=FALSE,jacobian=FALSE,hessian=FALSE){
    if (x$format %in% c(1,2,3))
        out <- data2ludwig_2D(x,lta0w=lta0b0w,exterr=exterr,
                              jacobian=jacobian,hessian=hessian)
    else if (x$format %in% c(4,5,6))
        out <- data2ludwig_3D(x,lta0b0w=lta0b0w,exterr=exterr,
                              jacobian=jacobian,hessian=hessian)
    else if (x$format %in% c(7,8))
        out <- data2ludwig_Th(x,lta0b0w=lta0b0w,exterr=exterr,
                              jacobian=jacobian,hessian=hessian)
    else stop('Incorrect input format.')
    out
}
data2ludwig_2D <- function(x,lta0w,exterr=FALSE,jacobian=FALSE,hessian=FALSE){
    out <- list()
    U <- iratio('U238U235')[1]
    np <- min(3,length(lta0w)) # number of parameters
    lt <- lta0w[1]
    a0 <- lta0w[2]
    if (np==3) w <- ta0w[3] # model 3 regression
    ns <- length(x)
    zeros <- rep(0,ns)
    E <- matrix(0,2*ns+6,2*ns+6)
    J <- matrix(0,2*ns,2*ns+6)
    J[1:(2*ns),1:(2*ns)] <- diag(2*ns)
    X <- zeros
    Y <- zeros
    D <- mclean(tt=exp(lt),d=x$d,exterr=exterr)
    for (i in 1:ns){
        wd <- wetherill(x,i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        E[(0:1)*ns+i,(0:1)*ns+i] <- wd$cov
        J[i,2*ns+2] <- -D$dPb207U235dl35     # dKdl35
        J[i,2*ns+4] <- -D$dPb207U235dl31     # dKdl31
        J[ns+i,2*ns+1] <- -D$dPb206U238dl38  # dLdl31
        J[ns+i,2*ns+3] <- -D$dPb206U238dl34  # dLdl34
        J[ns+i,2*ns+5] <- -D$dPb206U238dl30  # dLdl30
        J[ns+i,2*ns+6] <- -D$dPb206U238dl26  # dLdl26
    }
    E[2*ns+1,2*ns+1] <- lambda('U238')[2]^2
    E[2*ns+2,2*ns+2] <- lambda('U235')[2]^2
    E[2*ns+3,2*ns+3] <- (lambda('U234')[2]*1000)^2
    E[2*ns+4,2*ns+4] <- (lambda('Pa231')[2]*1000)^2
    E[2*ns+5,2*ns+5] <- (lambda('Th230')[2]*1000)^2
    E[2*ns+6,2*ns+6] <- (lambda('Ra226')[2]*1000)^2
    ED <- J%*%E%*%t(J)
    if (np==3){
        Ew <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D)
        ED <- ED + Ew
    }
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    O <- blockinverse(AA=ED[i1,i1],BB=ED[i1,i2],
                      CC=ED[i2,i1],DD=ED[i2,i2],doall=TRUE)
    K0 <- X - D$Pb207U235 + a0*U*(D$Pb206U238 - Y)
    A <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*exp(a0)*U + K0%*%(O[i1,i2]+t(O[i2,i1])))
    B <- -(exp(a0)*U*(O[i1,i1]+t(O[i1,i1]))*exp(a0)*U + (O[i2,i2]+t(O[i2,i2])) +
           exp(a0)*U*(O[i1,i2]+t(O[i1,i2])) + (O[i2,i1]+t(O[i2,i1]))*exp(a0)*U)
    L <- as.vector(solve(B,A))
    c0 <- Y - D$Pb206U238 - L
    K <- X - D$Pb207U235 - exp(a0)*U*c0
    KL <- c(K,L)
    out$SS <- KL%*%O%*%KL
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- -(2*ns*log(2*pi) + detED + out$SS)/2
    dtdlt <- lt
    if (jacobian | hessian){
        JKL <- matrix(0,2*ns,ns+2) # derivatives of KL w.r.t. c0, t, and a0
        colnames(JKL) <- c(paste0('c0[',i1,']'),'lt','a0')
        rownames(JKL) <- c(paste0('K[',i1,']'),paste0('L[',i1,']'))
        diag(JKL[i1,i1]) <- -U*exp(a0)        # dKdc0
        diag(JKL[i2,i1]) <- -1                # dLdc0
        JKL[i1,'lt']  <- -D$dPb207U235dt*dtdlt # dKdlt
        JKL[i2,'lt']  <- -D$dPb206U238dt*dtdlt # dLdlt
        JKL[i1,'a0'] <- -U*c0                 # dKda0
        out$jacobian <- rep(0,np) # derivatives of LL w.r.t t, a0 (and w)
        names(out$jacobian) <- c('lt','a0','w')[1:np]
        out$jacobian['lt'] <- -KL%*%O%*%JKL[,'lt']
        out$jacobian['a0'] <- -KL%*%O%*%JKL[,'a0']
        if (np==3){
            dEDdw <- get.Ew(w=w,format=x$format,ns=ns,lt=lt,D=D,deriv=1)
            dlnDetEDdw <- trace(O%*%dEDdw)
            dOdw <- -O%*%dEDdw%*%O
            out$jacobian['w'] <- -(dlnDetEDdw + KL%*%dOdw%*%KL)/2
        }
    }
    if (hessian){ # second derivatives of the NEGATIVE log-likelihood
        d2tdlt2 <- lt
        out$hessian <- matrix(0,ns+np,ns+np)
        hesnames <- c(paste0('c0[',1:ns,']'),'lt','a0','w')[1:(ns+np)]
        rownames(out$hessian) <- hesnames
        colnames(out$hessian) <- hesnames
        d2KLdt2 <- rep(0,2*ns)
        d2KLdt2[i1] <- -D$d2Pb207U235dt2         # d2Kdt2
        d2KLdt2[i2] <- -D$d2Pb206U238dt2         # d2Ldt2
        d2KLdlt2 <- d2Kdt2*dtdlt^2 + dKdt*d2tdlt2
        d2KLdc0da0 <- matrix(0,2*ns,ns)
        diag(d2KLdc0da0[i1,i1]) <- -U*exp(a0)    # d2Kdc0da0
        out$hessian[i1,i1] <- t(JKL[,i1])%*%O%*%JKL[,i1]                  # d2dc02
        out$hessian['lt',i1] <- t(JKL[,'lt'])%*%O%*%JKL[,i1]              # d2dc0dlt
        out$hessian['a0',i1] <-
            t(JKL[,'a0'])%*%O%*%JKL[,i1] + KL%*%O%*%d2KLdc0da0            # d2dc0da0
        out$hessian['lt','lt'] <-
            t(JKL[,'lt'])%*%O%*%JKL[,'lt'] + KL%*%O%*%d2KLdt2             # d2dt2
        out$hessian['a0','a0'] <- t(JKL[,'a0'])%*%O%*%JKL[,'a0']          # d2da02
        out$hessian['lt','a0'] <- t(JKL[,'lt'])%*%O%*%JKL[,'a0']          # d2dtda0
        out$hessian['a0','lt'] <- out$hessian['lt','a0']                  # d2da0dt
        out$hessian[i1,'lt'] <- out$hessian['lt',i1]                      # d2dc0dt
        out$hessian[i1,'a0'] <- out$hessian['a0',i1]                      # d2dc0da0
        if (np==3){
            d2EDdw2 <- get.Ew(w=w,format=x$format,ns=ns,lt=lt,D=D,deriv=2)
            d2lnDetEDdw2 <- trace(O%*%d2EDdw2) - trace(O%*%dEDdw%*%O%*%dEDdw) 
            d2Odw2 <- -(dOdw%*%dEDdw%*%O + O%*%d2EDdw2%*%O + O%*%dEDdw%*%dOdw)
            out$hessian['w','w'] <- (d2lnDetEDdw2 + KL%*%d2Odw2%*%KL)/2   # d2dw2
            out$hessian['w',i1] <- KL%*%dOdw%*%JKL[,i1]                   # d2dwdc0
            out$hessian['w','lt'] <- KL%*%dOdw%*%JKL[,'lt']               # d2dwdt
            out$hessian['w','a0'] <- KL%*%dOdw%*%JKL[,'a0']               # d2dwda0
            out$hessian[i1,'w'] <- out$hessian['w',i1]                    # d2dc0dw
            out$hessian['lt','w'] <- out$hessian['w','lt']                # d2dtdw
            out$hessian['a0','w'] <- out$hessian['w','a0']                # d2da0dw
        }
    }
    out
}
# rederived from Ludwig (1998):
data2ludwig_3D <- function(x,ta0b0w,exterr=FALSE,jacobian=FALSE,hessian=FALSE){
    out <- list()
    U <- iratio('U238U235')[1]
    np <- min(4,length(ta0b0w)) # number of parameters
    tt <- ta0b0w[1]
    a0 <- ta0b0w[2]
    b0 <- ta0b0w[3]
    if (np==4) w <- ta0b0w[4] # model 3 regression
    ns <- length(x)
    zeros <- rep(0,ns)
    E <- matrix(0,3*ns+6,3*ns+6)
    J <- matrix(0,3*ns,3*ns+6)
    J[1:(3*ns),1:(3*ns)] <- diag(3*ns)
    X <- zeros
    Y <- zeros
    Z <- zeros
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        Z[i] <- wd$x['Pb204U238']
        E[(0:2)*ns+i,(0:2)*ns+i] <- wd$cov
        J[i,3*ns+2] <- -D$dPb207U235dl35     #dKdl35
        J[i,3*ns+4] <- -D$dPb207U235dl31     #dKdl31
        J[ns+i,3*ns+1] <- -D$dPb206U238dl38  #dLdl38
        J[ns+i,3*ns+3] <- -D$dPb206U238dl34  #dLdl34
        J[ns+i,3*ns+5] <- -D$dPb206U238dl30  #dLdl30
        J[ns+i,3*ns+6] <- -D$dPb206U238dl26  #dLdl26
    }
    E[3*ns+1,3*ns+1] <- lambda('U238')[2]^2
    E[3*ns+2,3*ns+2] <- lambda('U235')[2]^2
    E[3*ns+3,3*ns+3] <- (lambda('U234')[2]*1000)^2
    E[3*ns+4,3*ns+4] <- (lambda('Pa231')[2]*1000)^2
    E[3*ns+5,3*ns+5] <- (lambda('Th230')[2]*1000)^2
    E[3*ns+6,3*ns+6] <- (lambda('Ra226')[2]*1000)^2
    ED <- J%*%E%*%t(J)
    if (np==4){
        Ew <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D)
        ED <- ED + Ew
    }
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    i3 <- (2*ns+1):(3*ns)
    O <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                         DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                         GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    o11 <- O[i1,i1]; o12 <- O[i1,i2]; o13 <- O[i1,i3]
    o21 <- O[i2,i1]; o22 <- O[i2,i2]; o23 <- O[i2,i3]
    o31 <- O[i3,i1]; o32 <- O[i3,i2]; o33 <- O[i3,i3]
    K0 <- X - D$Pb207U235 - U*b0*Z
    L0 <- Y - D$Pb206U238 - a0*Z
    V <- t(K0%*%(o11+t(o11))*U*b0 + L0%*%(o12+t(o21))*U*b0 + K0%*%(o12+t(o21))*a0 +
           L0%*%(o22+t(o22))*a0 + K0%*%(o13+t(o31)) + L0%*%(o23+t(o32)))
    W <- -(U*b0*(o11+t(o11))*U*b0 + U*b0*(o12+t(o12))*a0 + U*b0*(o13+t(o13)) +
           a0*(o21+t(o21))*U*b0 + a0*(o22+t(o22))*a0 + a0*(o23+t(o23)) +
           (o31+t(o31))*U*b0 + (o32+t(o32))*a0 + (o33+t(o33)))
    M <- as.vector(solve(W,V))
    c0 <- as.vector(Z - M)
    K <- as.vector(X - D$Pb207U235 - U*b0*c0)
    L <- as.vector(Y - D$Pb206U238 - a0*c0)
    KLM <- c(K,L,M)
    out$SS <- KLM%*%O%*%KLM
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- -(3*ns*log(2*pi) + detED + out$SS)/2
    if (jacobian | hessian){
        JKLM <- matrix(0,3*ns,ns+3) # derivatives of KLM w.r.t. c0, t, a0 and b0
        colnames(JKLM) <- c(paste0('c0[',i1,']'),'t','a0','b0')
        rownames(JKLM) <- c(paste0('K[',i1,']'),paste0('L[',i1,']'),paste0('M[',i1,']'))
        diag(JKLM[i1,i1]) <- -U*b0       # dKdc0
        diag(JKLM[i2,i1]) <- -a0         # dLdc0
        diag(JKLM[i3,i1]) <- -1          # dMdc0
        JKLM[i1,'t']  <- -D$dPb207U235dt # dKdt
        JKLM[i2,'t']  <- -D$dPb206U238dt # dLdt
        JKLM[i2,'a0'] <- -c0             # dLda0
        JKLM[i1,'b0'] <- -U*c0           # dKdb0
        out$jacobian <- rep(0,np) # derivatives of LL w.r.t t, a0, b0 (and w)
        names(out$jacobian) <- c('t','a0','b0','w')[1:np]
        out$jacobian['t'] <- -KLM%*%O%*%JKLM[,'t']
        out$jacobian['a0'] <- -KLM%*%O%*%JKLM[,'a0']
        out$jacobian['b0'] <- -KLM%*%O%*%JKLM[,'b0']
        if (np==4){
            dEDdw <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D,deriv=1)
            dlnDetEDdw <- trace(O%*%dEDdw)
            dOdw <- -O%*%dEDdw%*%O
            out$jacobian['w'] <- -(dlnDetEDdw + KLM%*%dOdw%*%KLM)/2
        }
    }
    if (hessian){
        out$hessian <- matrix(0,ns+np,ns+np)
        hesnames <- c(paste0('c0[',1:ns,']'),'t','a0','b0','w')[1:(ns+np)]
        rownames(out$hessian) <- hesnames
        colnames(out$hessian) <- hesnames
        d2KLMdt2 <- rep(0,3*ns)
        d2KLMdt2[i1] <- -D$d2Pb207U235dt2 # d2Kdt2
        d2KLMdt2[i2] <- -D$d2Pb206U238dt2 # d2Ldt2
        d2KLMdc0da0 <- matrix(0,3*ns,ns)
        diag(d2KLMdc0da0[i2,i1]) <- -1    # d2Ldc0da0
        d2KLMdc0db0 <- matrix(0,3*ns,ns)
        diag(d2KLMdc0db0[i1,i1]) <- -U    # d2Kdc0db0
        out$hessian[i1,i1] <- t(JKLM[,i1])%*%O%*%JKLM[,i1]                # d2dc02
        out$hessian['t','t'] <-
            t(JKLM[,'t'])%*%O%*%JKLM[,'t'] + KLM%*%O%*%d2KLMdt2           # d2dt2
        out$hessian['t','a0'] <- t(JKLM[,'t'])%*%O%*%JKLM[,'a0']          # d2dtda0
        out$hessian['t','b0'] <- t(JKLM[,'t'])%*%O%*%JKLM[,'b0']          # d2dtdb0
        out$hessian['a0','a0'] <- t(JKLM[,'a0'])%*%O%*%JKLM[,'a0']        # d2da02
        out$hessian['b0','b0'] <- t(JKLM[,'b0'])%*%O%*%JKLM[,'b0']        # d2db02
        out$hessian['a0','b0'] <- t(JKLM[,'a0'])%*%O%*%JKLM[,'b0']        # d2da0db0
        out$hessian['t',i1] <- t(JKLM[,'t'])%*%O%*%JKLM[,i1]              # d2dc0dt
        out$hessian['a0',i1] <-
            t(JKLM[,'a0'])%*%O%*%JKLM[,i1] + KLM%*%O%*%d2KLMdc0da0        # d2dc0da0
        out$hessian['b0',i1] <-
            t(JKLM[,'b0'])%*%O%*%JKLM[,i1] + KLM%*%O%*%d2KLMdc0db0        # d2dc0db0
        out$hessian['a0','t'] <- out$hessian['t','a0']                    # d2da0dt
        out$hessian['b0','t'] <- out$hessian['t','b0']                    # d2db0dt
        out$hessian['b0','a0'] <- out$hessian['a0','b0']                  # d2db0da0
        out$hessian[i1,'t'] <- out$hessian['t',i1]                        # d2dtdc0
        out$hessian[i1,'a0'] <- out$hessian['a0',i1]                      # d2da0dc0
        out$hessian[i1,'b0'] <- out$hessian['b0',i1]                      # d2db0dc0
        if (np==4){
            d2EDdw2 <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D,deriv=2)
            d2lnDetEDdw2 <- trace(O%*%d2EDdw2) - trace(O%*%dEDdw%*%O%*%dEDdw) 
            d2Odw2 <- -(dOdw%*%dEDdw%*%O + O%*%d2EDdw2%*%O + O%*%dEDdw%*%dOdw)
            out$hessian['w','w'] <- (d2lnDetEDdw2 + KLM%*%d2Odw2%*%KLM)/2 # d2dw2
            out$hessian['w',i1] <- KLM%*%dOdw%*%JKLM[,i1]                 # d2dwdc0
            out$hessian['w','t'] <- KLM%*%dOdw%*%JKLM[,'t']               # d2dwdt
            out$hessian['w','a0'] <- KLM%*%dOdw%*%JKLM[,'a0']             # d2dwda0
            out$hessian['w','b0'] <- KLM%*%dOdw%*%JKLM[,'b0']             # d2dwdb0
            out$hessian[i1,'w'] <- out$hessian['w',i1]                    # d2dc0dw
            out$hessian['t','w'] <- out$hessian['w','t']                  # d2dtdw
            out$hessian['a0','w'] <- out$hessian['w','a0']                # d2da0dw
            out$hessian['b0','w'] <- out$hessian['w','b0']                # d2db0dw
        }
    }
    out
}
data2ludwig_Th <- function(x,ta0b0w,exterr=FALSE,jacobian=FALSE,hessian=FALSE){
    out <- list()
    np <- min(4,length(ta0b0w)) # number of parameters
    tt <- ta0b0w[1]
    a0 <- ta0b0w[2]
    b0 <- ta0b0w[3]
    if (np==4) w <- ta0b0w[4] # model 3 regression
    U <- iratio('U238U235')[1]
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    ns <- length(x)
    zeros <- rep(0,ns)
    X <- zeros
    Y <- zeros
    Z <- zeros
    W <- zeros
    K0 <- zeros
    L0 <- zeros
    E <- matrix(0,4*ns+7,4*ns+7)
    J <- matrix(0,3*ns,4*ns+7)
    J[1:(3*ns),1:(3*ns)] <- diag(3*ns)
    for (i in 1:ns){
        wd <- wetherill(x,i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        Z[i] <- wd$x['Pb208Th232']
        W[i] <- wd$x['Th232U238']
        K0[i] <- X[i] - D$Pb207U235 - (Z[i]-D$Pb208Th232)*exp(b0)*U*W[i]
        L0[i] <- Y[i] - D$Pb206U238 - (Z[i]-D$Pb208Th232)*exp(a0)*W[i]
        E[(0:3)*ns+i,(0:3)*ns+i] <- wd$cov
        J[i,4*ns+2] <- -D$dPb207U235dl35       # dKdl35
        J[i,4*ns+5] <- -D$dPb207U235dl31       # dKdl31
        J[ns+i,4*ns+1] <- -D$dPb206U238dl38    # dLdl38
        J[ns+i,4*ns+3] <- -D$dPb206U238dl34    # dLdl34
        J[ns+i,4*ns+6] <- -D$dPb206U238dl30    # dLdl30
        J[ns+i,4*ns+7] <- -D$dPb206U238dl26    # dLdl26
        J[2*ns+i,4*ns+4] <- -D$dPb208Th232dl32 # dMdl32
    }
    E[4*ns+1,4*ns+1] <- lambda('U238')[2]^2
    E[4*ns+2,4*ns+2] <- lambda('U235')[2]^2
    E[4*ns+3,4*ns+3] <- (lambda('U234')[2]*1000)^2
    E[4*ns+4,4*ns+4] <- lambda('Th232')[2]^2
    E[4*ns+5,4*ns+5] <- (lambda('Pa231')[2]*1000)^2
    E[4*ns+6,4*ns+6] <- (lambda('Th230')[2]*1000)^2
    E[4*ns+7,4*ns+7] <- (lambda('Ra226')[2]*1000)^2
    ED <- J%*%E%*%t(J)
    if (np==4){
        Ew <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D)
        ED <- ED + Ew
    }
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    i3 <- (2*ns+1):(3*ns)
    O <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                         DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                         GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    o11 <- O[i1,i1]; o12 <- O[i1,i2]; o13 <- O[i1,i3]
    o21 <- O[i2,i1]; o22 <- O[i2,i2]; o23 <- O[i2,i3]
    o31 <- O[i3,i1]; o32 <- O[i3,i2]; o33 <- O[i3,i3]
    Wd <- diag(W)
    AA <- (Wd%*%o11%*%Wd)*(U*exp(b0))^2 + (Wd%*%o22%*%Wd)*exp(a0)^2 + o33 +
        U*exp(a0)*exp(b0)*Wd%*%(o12+o21)%*%Wd +
        U*exp(b0)*(Wd%*%o13+o31%*%Wd) + exp(a0)*(Wd%*%o23+o32%*%Wd)
    BT <- t(U*exp(b0)*K0%*%o11%*%Wd + exp(a0)*L0%*%o22%*%Wd +
            exp(a0)*K0%*%o12%*%Wd + U*exp(b0)*L0%*%o21%*%Wd +
            K0%*%o13 + L0%*%o23)
    CC <- U*exp(b0)*Wd%*%o11%*%K0 + exp(a0)*Wd%*%o22%*%L0 +
        exp(a0)*Wd%*%o12%*%K0 + U*exp(b0)*Wd%*%o21%*%L0 +
        o13%*%K0 + o23%*%L0
    M <- as.vector(solve(-(AA+t(AA)),(BT+CC)))
    c0 <- as.vector(Z - D$Pb208Th232 - M)
    K <- as.vector(X - D$Pb207U235 - c0*exp(b0)*U*W)
    L <- as.vector(Y - D$Pb206U238 - c0*exp(a0)*W)
    KLM <- c(K,L,M)
    out$SS <- KLM%*%O%*%KLM
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- -(3*log(2*pi) + detED + out$SS)/2
    if (jacobian | hessian){
        JKLM <- matrix(0,3*ns,ns+3)
        colnames(JKLM) <- c(paste0('c0[',i1,']'),'t','a0','b0')
        rownames(JKLM) <- c(paste0('K[',i1,']'),paste0('L[',i1,']'),paste0('M[',i1,']'))
        diag(JKLM[i1,i1]) <- -exp(b0)*U*W  # dKdc0
        diag(JKLM[i2,i1]) <- -exp(a0)*W    # dLdc0
        diag(JKLM[i3,i1]) <- -1            # dMdc0
        JKLM[i1,'t']  <- -D$dPb207U235dt   # dKdt
        JKLM[i2,'t']  <- -D$dPb206U238dt   # dLdt
        JKLM[i3,'t']  <- -D$dPb208Th232dt  # dMdt
        JKLM[i2,'a0'] <- -c0*W*exp(a0)     # dLda0
        JKLM[i1,'b0'] <- -c0*U*W*exp(b0)   # dKdb0
        out$jacobian <- rep(0,np)
        names(out$jacobian) <- c('t','a0','b0','w')[1:np]
        out$jacobian['t'] <- -KLM%*%O%*%JKLM[,'t']
        out$jacobian['a0'] <- -KLM%*%O%*%JKLM[,'a0']
        out$jacobian['b0'] <- -KLM%*%O%*%JKLM[,'b0']
        if (np==4){
            dEDdw <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D,deriv=1)
            dlnDetEDdw <- trace(O%*%dEDdw)
            dOdw <- -O%*%dEDdw%*%O
            out$jacobian['w'] <- -(dlnDetEDdw + KLM%*%dOdw%*%KLM)/2
        }
    }
    if (hessian){
        out$hessian <- matrix(0,ns+np,ns+np)
        hesnames <- c(paste0('c0[',1:ns,']'),'t','a0','b0','w')[1:(ns+np)]
        rownames(out$hessian) <- hesnames
        colnames(out$hessian) <- hesnames
        d2KLMdt2 <- rep(0,3*ns)
        d2KLMdt2[i1] <- -D$d2Pb207U235dt2        # d2Kdt2
        d2KLMdt2[i2] <- -D$d2Pb206U238dt2        # d2Ldt2
        d2KLMdt2[i3] <- -D$d2Pb208Th232dt2       # d2Mdt2
        d2KLMda02 <- rep(0,3*ns)
        d2KLMda02[i2] <- -c0*W*exp(a0)           # d2Lda02
        d2KLMdb02 <- rep(0,3*ns)
        d2KLMdb02[i1] <- -c0*U*W*exp(b0)         # d2Kdb02
        d2KLMdc0da0 <- matrix(0,3*ns,ns)
        diag(d2KLMdc0da0[i2,i1]) <- -W*exp(a0)   # d2Ldc0da0
        d2KLMdc0db0 <- matrix(0,3*ns,ns)
        diag(d2KLMdc0db0[i1,i1]) <- -U*W*exp(b0) # d2Kdc0db0
        out$hessian[i1,i1] <- t(JKLM[,i1])%*%O%*%JKLM[,i1]                # d2dc02
        out$hessian['t','t'] <-
            t(JKLM[,'t'])%*%O%*%JKLM[,'t'] + KLM%*%O%*%d2KLMdt2           # d2dt2
        out$hessian['a0','a0'] <-
            t(JKLM[,'a0'])%*%O%*%JKLM[,'a0'] + KLM%*%O%*%d2KLMda02        # d2da02
        out$hessian['b0','b0'] <-
            t(JKLM[,'b0'])%*%O%*%JKLM[,'b0'] + KLM%*%O%*%d2KLMdb02        # d2db02
        out$hessian['t','a0'] <- t(JKLM[,'t'])%*%O%*%JKLM[,'a0']          # d2dtda0
        out$hessian['t','b0'] <- t(JKLM[,'t'])%*%O%*%JKLM[,'b0']          # d2dtdb0
        out$hessian['a0','b0'] <- t(JKLM[,'a0'])%*%O%*%JKLM[,'b0']        # d2da0db0
        out$hessian['t',i1] <- t(JKLM[,'t'])%*%O%*%JKLM[,i1]              # d2dc0dt
        out$hessian['a0',i1] <-
            t(JKLM[,'a0'])%*%O%*%JKLM[,i1] + KLM%*%O%*%d2KLMdc0da0        # d2dc0da0
        out$hessian['b0',i1] <-
            t(JKLM[,'b0'])%*%O%*%JKLM[,i1] + KLM%*%O%*%d2KLMdc0db0        # d2dc0db0
        out$hessian['a0','t'] <- out$hessian['t','a0']                    # d2da0dt
        out$hessian['b0','t'] <- out$hessian['t','b0']                    # d2db0dt
        out$hessian['b0','a0'] <- out$hessian['a0','b0']                  # d2db0da0
        out$hessian[i1,'t'] <- out$hessian['t',i1]                        # d2dtdc0
        out$hessian[i1,'a0'] <- out$hessian['a0',i1]                      # d2da0dc0
        out$hessian[i1,'b0'] <- out$hessian['b0',i1]                      # d2db0dc0
        if (np==4){
            d2EDdw2 <- get.Ew(w=w,format=x$format,ns=ns,tt=tt,D=D,deriv=2)
            d2lnDetEDdw2 <- trace(O%*%d2EDdw2) - trace(O%*%dEDdw%*%O%*%dEDdw) 
            d2Odw2 <- -(dOdw%*%dEDdw%*%O + O%*%d2EDdw2%*%O + O%*%dEDdw%*%dOdw)
            out$hessian['w','w'] <- (d2lnDetEDdw2 + KLM%*%d2Odw2%*%KLM)/2 # d2dw2
            out$hessian['w',i1] <- KLM%*%dOdw%*%JKLM[,i1]                 # d2dwdc0
            out$hessian['w','t'] <- KLM%*%dOdw%*%JKLM[,'t']               # d2dwdt
            out$hessian['w','a0'] <- KLM%*%dOdw%*%JKLM[,'a0']             # d2dwda0
            out$hessian['w','b0'] <- KLM%*%dOdw%*%JKLM[,'b0']             # d2dwdb0
            out$hessian[i1,'w'] <- out$hessian['w',i1]                    # d2dc0dw
            out$hessian['t','w'] <- out$hessian['w','t']                  # d2dtdw
            out$hessian['a0','w'] <- out$hessian['w','a0']                # d2da0dw
            out$hessian['b0','w'] <- out$hessian['w','b0']                # d2db0dw
        }
    }
    out
}

get.Ew <- function(w=0,format=1,ns=1,tt=0,D=mclean(),deriv=0){
    if (format<4) ndim <- 2
    else ndim <- 3
    J <- matrix(0,ndim,1)
    J[1,1] <- -D$dPb207U235dt # dKdt
    J[2,1] <- -D$dPb206U238dt # dLdt
    if (format>6) J[3,1] <- -D$dPb208Th232dt # dMdt
    if (deriv==1) {
        dEdw <- 2*exp(w)
        Ew <- J%*%dEdw%*%t(J)
    } else if (deriv==2) {
        d2Edw2 <- 2*exp(w)
        Ew <- J%*%d2Edw2%*%t(J)
    } else {
        E <- exp(w)^2
        Ew <- J%*%E%*%t(J)
    }
    out <- matrix(0,ndim*ns,ndim*ns)
    diag(out[1:ns,1:ns]) <- Ew[1,1]
    diag(out[(ns+1):(2*ns),(ns+1):(2*ns)]) <- Ew[2,2]
    diag(out[1:ns,(ns+1):(2*ns)]) <- Ew[1,2]
    diag(out[(ns+1):(2*ns),1:ns]) <- Ew[2,1]
    if (format>4){
        diag(out[(ns+1):(2*ns),(2*ns+1):(3*ns)]) <- Ew[3,3]
        diag(out[1:ns,(2*ns+1):(3*ns)]) <- Ew[1,3]
        diag(out[(2*ns+1):(3*ns),1:ns]) <- Ew[3,1]
        diag(out[(ns+1):(2*ns),(2*ns+1):(3*ns)]) <- Ew[2,3]
        diag(out[(2*ns+1):(3*ns),(ns+1):(2*ns)]) <- Ew[3,2]
    }
    out
}

fixit <- function(x,anchor=list(FALSE,NA),model=1,w=NA){
    if (x$format<4) np <- 2
    else np <- 3
    if (model==3) np <- np+1
    out <- rep(FALSE,np)
    if (model==3 & is.numeric(w)) out[np] <- TRUE # fix w
    if (anchor[[1]]){ # anchor t or a0(,b0)
        if (is.numeric(anchor[[2]])) out[1] <- TRUE # fix t
        else out[2:(np-1)] <- TRUE # fix a0(,b0)
    }
    out
}
