#' Leftmost and rightmost isochrons
#'
#' Modifies the Galbraith and Laslett's minimum age module to fit
#' overdispersed Pb-Pb and Th-U isochron data.
#' @param x an IsoplotR data object
#' @param left logical, switches between leftmost and rightmost
#'     isochrons
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the isochron age calculation.
#' @param inverse toggles between normal and inverse isochrons. See
#'     \code{\link{isochron}}.
#' @param ... unused optional arguments
#' @return a list with the intercept (\code{a}), slope (\code{b}) and
#'     other data that can be passed on to \code{scatterplot}.
#' @noRd
LRisochron <- function(x,...){ UseMethod("LRisochron",x) }
#' @noRd
LRisochron.default <- function(x,left=TRUE,hide=NULL,omit=NULL,...){
    x2calc <- clear(x,hide,omit)
    if (left){
        init_fit <- init_left(yd=x2calc)
    } else {
        init_fit <- init_right(yd=x2calc)
    }
    y0i <- init_fit$a
    gi <- -init_fit$b
    propi <- 0.5
    sigi <- stats::sd((y0i-x2calc[,'Y'])/x2calc[,'X'])
    init <- c(gi,propi,sigi,y0i)
    lower <- c(-Inf,0,.Machine$double.eps,-Inf)
    upper <- c(Inf,1,Inf,Inf)
    fun <- ifelse(left,yd2ratios_left,yd2ratios_right)
    fit <- contingencyfit(par=init,fn=get_LRisochron_L,
                          lower=lower,upper=upper,
                          yd=x2calc,fun=fun,hessian=TRUE)
    covmat <- inverthess(fit$hessian)
    a <- fit$par[4]
    sa <- sqrt(covmat[4,4])
    b <- ifelse(left,-1,1) * fit$par[1]
    sb <- sqrt(covmat[1,1])
    list(a=c('a'=unname(a),'s[a]'=unname(sa)),
         b=c('b'=unname(b),'s[b]'=unname(sb)),
         cov.ab=unname(covmat[1,4]),
         xyz=x,model=4,n=nrow(x2calc))
}
#' @noRd
LRisochron.PbPb <- function(x,inverse=TRUE,anchor=0,hide=NULL,omit=NULL,...){
    yd <- data2york(x,inverse=TRUE)
    yd2calc <- clear(yd,hide,omit)
    init_fit <- init_left(yd2calc)
    y0i <- init_fit$a
    gi <- -init_fit$b
    propi <- 0.5
    sigi <- stats::sd((y0i-yd2calc[,'Y'])/yd2calc[,'X'])
    x1 <- y1 <- y0 <- NULL
    eps <- .Machine$double.eps
    if (anchor[1]==1 & length(anchor)>1){
        y0 <- age2ratio(tt=anchor[2],ratio="Pb207Pb206")[1]
        init <- c(gi,propi,sigi)
        lower <- c(-Inf,0,eps)
        upper <- c(Inf,1,Inf)
    } else if (anchor[1]==2){
        Pb74 <- iratio('Pb207Pb204')[1]
        Pb64 <- iratio('Pb206Pb204')[1]
        y1 <- Pb74/Pb64
        x1 <- 1/Pb64
        init <- c(propi,sigi,y0i)
        lower <- c(0,eps,-Inf)
        upper <- c(1,Inf,Inf)
    } else {
        init <- c(gi,propi,sigi,y0i)
        lower <- c(-Inf,0,eps,-Inf)
        upper <- c(Inf,1,Inf,Inf)
    }
    fit <- contingencyfit(par=init,fn=get_LRisochron_L,
                          lower=lower,upper=upper,
                          fun=yd2ratios_left,
                          yd=yd2calc,y0=y0,x1=x1,y1=y1,
                          hessian=TRUE)
    covmat <- inverthess(fit$hessian)
    np <- length(fit$par)
    if (is.null(y0)){
        a <- fit$par[np]
        sa <- sqrt(covmat[np,np])
    } else {
        a <- y0
        sa <- 0
        cov.ab <- 0
    }
    if (is.null(x1) | is.null(y1)){
        b <- -fit$par[1]
        sb <- sqrt(covmat[1,1])
        if (is.null(y0)) cov.ab <- covmat[1,np]
    } else {
        b <- (y1-a)/x1
        sb <- sa/x1
        rho <- -1
        cov.ab <- rho*sa*sb
    }
    out <- list(a=c('a'=unname(a),'s[a]'=unname(sa)),
                b=c('b'=unname(b),'s[b]'=unname(sb)),
                cov.ab=unname(cov.ab),
                model=4,n=nrow(yd2calc))
    if (inverse){
        out$xyz <- yd
    } else {
        out <- invertfit(out,type="d")
        out$xyz <- normal2inverse(yd,type='d')
    }
    out
}
#' @noRd
LRisochron.ThU <- function(x,inverse=TRUE,anchor=0,hide=NULL,omit=NULL,...){
    if (x$format<3){
        stop("Rightmost isochrons are only available for ThU formats 3 and 4.")
    }
    yd <- data2york(x,inverse=FALSE)
    yd2calc <- clear(yd,hide,omit)
    init <- init_right(yd2calc)
    gi <- -init$b
    y0i <- ifelse(anchor[1]==1,1/x$U8Th2,init$a)
    propi <- 0.5
    sigi <- stats::sd((yd2calc[,'Y']-y0i)/(yd2calc[,'X']-y0i))
    init <- c(propi,sigi)
    lower <- c(0,.Machine$double.eps)
    upper <- c(1,Inf)
    y0 <- b <- NULL
    if (anchor[1]==1){
        y0 <- y0i
        lower <- c(-Inf,lower)
        upper <- c(Inf,upper)
        init <- c(gi,init)
    } else if (anchor[1]==2 & length(anchor)>1){
        b <- age2ratio(tt=anchor[2],ratio='Th230U238')[1]
        lower <- c(lower,-Inf)
        upper <- c(upper,Inf)
        init <- c(init,y0i)
    } else {
        lower <- c(-Inf,lower,-Inf)
        upper <- c(Inf,upper,Inf)
        init <- c(gi,init,y0i)
    }
    fit <- contingencyfit(par=init,fn=get_LRisochron_L,
                          lower=lower,upper=upper,
                          yd=yd2calc,y0=y0,gam=b,
                          fun=yd2ratios_ThU,
                          hessian=TRUE)
    covmat <- inverthess(fit$hessian)
    np <- length(fit$par)
    if (anchor[1]==1){
        sy0 <- 0
        b <- fit$par[1]
        sb <- sqrt(covmat[1,1])
        J <- rbind(c(-y0,1-b),
                   c(1,0))
        E <- J %*% rbind(c(covmat[1,1],0),c(0,0)) %*% t(J)
        a <- y0*(1-b)
        sa <- sqrt(E[1,1])
        cov.ab <- E[1,2]
    } else if (anchor[1]==2 & length(anchor)>1){
        y0 <- fit$par[np]
        sy0 <- sqrt(covmat[np,np])
        sb <- 0
        a <- y0*(1-b)
        sa <- sy0*(1-b)
        cov.ab <- 0
    } else {
        y0 <- fit$par[np]
        sy0 <- sqrt(covmat[np,np])
        b <- fit$par[1]
        sb <- sqrt(covmat[1,1])
        J <- rbind(c(-y0,1-b),
                   c(1,0))
        E <- J %*% covmat[c(1,np),c(1,np)] %*% t(J)
        a <- y0*(1-b)
        sa <- sqrt(E[1,1])
        cov.ab <- E[1,2]
    }
    out <- list(a=c('a'=unname(a),'s[a]'=unname(sa)),
                b=c('b'=unname(b),'s[b]'=unname(sb)),
                cov.ab=unname(cov.ab),
                model=4,n=nrow(yd2calc))
    if (inverse){
        out <- invertfit(out,type='d')
        out$xyz <- normal2inverse(yd,type='d')
    } else {
        out$xyz <- yd
    }
    out
}

init_left <- function(yd){
    X <- yd[,'X']
    Y <- yd[,'Y']
    leftmost <- which.min(X)
    topmost <- which.max(Y)
    if (leftmost==topmost){
        fit <- lm(Y ~ X)
        a <- fit$coefficients[1]
        b <- fit$coefficients[2]
    } else {
        a <- Y[leftmost]
        b <- (Y[topmost]-Y[leftmost])/(X[topmost]-X[leftmost])
    }
    list(a=a,b=b)
}
init_right <- function(yd){
    bottommost <- which.min(yd[,'Y'])
    list(a=yd[bottommost,'Y'],b=0)
}

yd2ratios_left <- function(yd,y0){
    r <- (y0-yd[,'Y'])/yd[,'X']
    sr <- sqrt(errorprop1x2(J1=-r/yd[,'X'],
                            J2=-1/yd[,'X'],
                            E11=yd[,'sX']^2,
                            E22=yd[,'sY']^2,
                            E12=yd[,'rXY']*yd[,'sX']*yd[,'sY']))
    cbind(r,sr)
}
yd2ratios_right <- function(yd,y0){
    r <- (yd[,'Y']-y0)/yd[,'X']
    sr <- sqrt(errorprop1x2(J1=-r/yd[,'X'],
                            J2=1/yd[,'X'],
                            E11=yd[,'sX']^2,
                            E22=yd[,'sY']^2,
                            E12=yd[,'rXY']*yd[,'sX']*yd[,'sY']))
    cbind(r,sr)    
}
yd2ratios_ThU <- function(yd,y0){
    r <- (yd[,'Y']-y0)/(yd[,'X']-y0)
    sr <- sqrt(errorprop1x2(J1=-r/(yd[,'X']-y0),
                            J2=-y0/(yd[,'X']-y0),
                            E11=yd[,'sX']^2,
                            E22=yd[,'sY']^2,
                            E12=yd[,'rXY']*yd[,'sX']*yd[,'sY']))
    cbind(r,sr)
}

get_LRisochron_L <- function(pars,yd,y0=NULL,
                             x1=NULL,y1=NULL,gam=NULL,
                             fun=yd2ratios_left){
    np <- length(pars)
    if (is.null(y0)){
        y0 <- pars[np]
    }
    if (is.null(gam)){
        if (is.null(x1) || is.null(y1)){
            gam <- pars[1]
            prop <- pars[2]
            sig <- pars[3]
        } else {
            gam <- (y0-y1)/x1
            prop <- pars[1]
            sig <- pars[2]
        }
    } else {
        prop <- pars[1]
        sig <- pars[2]
    }
    mu <- gam
    zs <- fun(yd=yd,y0=y0)
    z <- zs[,1]
    s <- zs[,2]
    AA  <- prop/sqrt(2*pi*s^2)
    BB <- -0.5*((z-gam)/s)^2
    CC <- (1-prop)/sqrt(2*pi*(sig^2+s^2))
    mu0 <- (mu/sig^2 + z/s^2)/(1/sig^2 + 1/s^2)
    s0 <- 1/sqrt(1/sig^2 + 1/s^2)
    DD <- 1-stats::pnorm((gam-mu0)/s0)
    EE <- 1-stats::pnorm((gam-mu)/sig)
    FF <- -0.5*((z-mu)^2)/(sig^2+s^2)
    fu <- AA*exp(BB) + CC*(DD/EE)*exp(FF)
    fu[fu<.Machine$double.xmin] <- .Machine$double.xmin
    fu[fu>.Machine$double.xmax] <- .Machine$double.xmax
    sum(-log(fu))
}

