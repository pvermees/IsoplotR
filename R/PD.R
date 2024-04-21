getDPratio <- function(tt,st=0,nuclide,exterr=FALSE,bratio=1){
    L <- lambda(nuclide)
    R <- bratio*(exp(L[1]*tt)-1)
    J1 <- bratio*L[1]*exp(L[1]*tt)
    if (exterr) J2 <- bratio*tt*exp(L[1]*tt)
    else J2 <- 0
    E11 <- st^2
    E22 <- L[2]^2
    vR <- errorprop1x2(J1,J2,E11,E22)
    out <- cbind(R,sqrt(vR))
    colnames(out) <- c('PD','s[PD]')
    out
}

getPDage <- function(DP,sDP=0,nuclide,exterr=FALSE,bratio=1){
    L <- lambda(nuclide)
    tt <- log(1 + DP/bratio)/L[1]
    E11 <- sDP^2
    if (exterr) E22 <- L[2]^2
    else E22 <- 0
    J1 <- 1/(L[1]*(bratio + DP)) # dt.dDP
    J2 <- -tt/L[1]               # dt.dL
    vt <- errorprop1x2(J1,J2,E11,E22)
    out <- cbind(tt,sqrt(vt))
    colnames(out) <- c('t','s[t]')
    out
}

# i2i = isochron to intercept
# bratio = branching ratio
# projerr = isochron projection error
PD_age <- function(x,nuclide,exterr=FALSE,i=NULL,i2i=TRUE,
                   bratio=1,omit4c=NULL,projerr=FALSE,...){
    ns <- length(x)
    out <- matrix(0,ns,2)
    if (ns<2) i2i <- FALSE
    if (i2i){
        y <- data2york(x,inverse=TRUE)
        fit <- regression(y,model=1,omit=omit4c)
        DP <- matrix(0,ns,2)
        DP[,1] <- 1/(y[,'X'] - y[,'Y']/fit$b[1])
        J1 <- -DP[,1]^2
        J2 <- (DP[,1]^2)/fit$b[1]
        J3 <- y[,'Y']/(y[,'Y']^2 - (fit$b[1]*y[,'X'])^2)
        E11 <- y[,'sX']^2
        E22 <- y[,'sY']^2
        E12 <- y[,'rXY']*y[,'sX']*y[,'sY']
        E33 <- fit$b[2]^2
        if (projerr) DP[,2] <- sqrt(errorprop1x3(J1,J2,J3,E11,E22,E33,E12))
        else DP[,2] <- sqrt(errorprop1x2(J1,J2,E11,E22,E12))
    } else {
        initial <- get_nominal_initials(x)
        dat <- data2york(x,exterr=exterr)
        dat[,'Y'] <- dat[,'Y'] - initial$y0
        if (projerr) dat[,'sY'] <- sqrt(dat[,'sY']^2 + initial$sy0^2)
        DP <- quotient(dat[,'X'],dat[,'sX'],dat[,'Y'],dat[,'sY'],dat[,'rXY'])
    }
    out <- getPDage(subset(DP,select=1),subset(DP,select=2),
                     nuclide,exterr=exterr,bratio=bratio,...)
    if (!is.null(i)) out <- out[i,,drop=FALSE]
    colnames(out) <- c('t','s[t]')
    out
}

get_nominal_initials <- function(x){
    if (is.RbSr(x)){
        out <- settings('iratio','Sr87Sr86')
    } else if (is.SmNd(x)){
        out <- settings('iratio','Nd143Nd144')
    } else if (is.ReOs(x)){
        Os72 <- settings('iratio','Os187Os192')
        Os82 <- settings('iratio','Os188Os192')
        out <- quotient(Os82[1],Os82[2],Os72[1],Os72[2],0)
    } else if (is.LuHf(x)){
        out <- settings('iratio','Hf176Hf177')
    } else if (is.KCa(x)){
        out <- settings('iratio','Ca40Ca44')
    } else if (is.ThPb(x)){
        out <- settings('iratio','Pb208Pb204')
    } else if (is.ArAr(x)){
        out <- settings('iratio','Ar40Ar36')
    } else {
        out <- c(0,0)
    }
    list(y0=out[1],sy0=out[2])
}

#' Convert isotope dilution derived concentrations to ratios
#' @param x an IsoplotR data object with \code{x$format=3}
#' @noRd
ppm2ratios <- function(x,...){ UseMethod("ppm2ratios",x) }
#' @noRd
ppm2ratios.default <- function(x,...){
    stop('Method ppm2ratios not available for this class.')
}
