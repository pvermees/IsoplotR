get.PD.ratio <- function(tt,st,nuclide,exterr=TRUE,bratio=1){
    L <- lambda(nuclide)
    R <- bratio*(exp(L[1]*tt)-1)
    Jac <- matrix(0,1,2)
    E <- matrix(0,2,2)
    Jac[1,1] <- bratio*L[1]*exp(L[1]*tt)
    if (exterr) Jac[1,2] <- bratio*tt*exp(L[1]*tt)
    E[1,1] <- st^2
    E[2,2] <- L[2]^2
    sR <- sqrt(Jac %*% E %*% t(Jac))
    out <- c(R,sR)
}

get.PD.age <- function(DP,sDP,nuclide,exterr=TRUE,bratio=1){
    L <- lambda(nuclide)
    tt <- log(1 + DP/bratio)/L[1]
    E <- matrix(0,2,2)
    J <- matrix(0,1,2)
    E11 <- sDP^2
    if (exterr) E22 <- L[2]^2
    else E22 <- 0
    E12 <- 0
    J1 <- 1/(L[1]*(bratio + DP)) # dt.dDP
    J2 <- -tt/L[1]               # dt.dL
    st <- errorprop1x2(J1,J2,E11,E22,E12)
    out <- cbind(tt,st)
    colnames(out) <- c('t','s[t]')
    out
}

# i2i = isochron to intercept
# bratio = branching ratio
PD.age <- function(x,nuclide,exterr=TRUE,i=NA,
                   sigdig=NA,i2i=TRUE,bratio=1,...){
    ns <- length(x)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    if (ns<2) i2i <- FALSE
    if (i2i){
        y <- data2york(x,inverse=TRUE)
        fit <- regression(y,model=1)
        DP <- matrix(0,ns,2)
        DP[,1] <- 1/(y[,'X'] - y[,'Y']/fit$b[1])
        J1 <- -DP[,1]^2
        J2 <- (DP[,1]^2)/fit$b[1]
        E11 <- y[,'sX']^2
        E22 <- y[,'sY']^2
        E12 <- y[,'rXY']*y[,'sX']*y[,'sY']
        DP[,2] <- errorprop1x2(J1,J2,E11,E22,E12)
    } else {
        initial <- get.nominal.initials(x)
        dat <- data2york(x,exterr=exterr)
        dat[,'Y'] <- dat[,'Y'] - initial$y0
        DP <- quotient(dat[,'X'],dat[,'sX'],dat[,'Y'],dat[,'sY'],dat[,'rXY'])
        if (exterr) dat[,'sY'] <- sqrt(dat[,'sY']^2 + initial$sy0^2)
    }
    tt <- get.PD.age(subset(DP,select=1),subset(DP,select=2),
                     nuclide,exterr=exterr,bratio=bratio)
    out <- roundit(subset(tt,select=1),subset(tt,select=2),sigdig=sigdig)
    if (!is.na(i)) out <- out[i,]
    out
}

get.nominal.initials <- function(x){
    out <- c(0,0)
    if (hasClass(x,'RbSr')){
        out <- settings('iratio','Sr87Sr86')
    } else if (hasClass(x,'SmNd')){
        out <- settings('iratio','Nd143Nd144')
    } else if (hasClass(x,'ReOs')){
        Os72 <- settings('iratio','Os187Os192')
        Os82 <- settings('iratio','Os188Os192')
        out <- quotient(Os82[1],Os82[2],Os72[1],Os72[2],0)
    } else if (hasClass(x,'LuHf')){
        out <- settings('iratio','Hf176Hf177')
    } else if (hasClass(x,'KCa')){
        out <- settings('iratio','Ca40Ca44')
    }
    list(y0=out[1],sy0=out[2])
}

ppm2ratios <- function(x,...){ UseMethod("ppm2ratios",x) }
ppm2ratios.default <- function(x,...){
    stop('Method ppm2ratios not available for this class.')
}

PD.inverse.ratios <- function(x,exterr=FALSE){
    X <- PD.normal.ratios(x,exterr=exterr)
    normal2inverse(X)
}
PD.normal.ratios <- function(x,exterr=FALSE){
    if (x$format==1)
        out <- x$x
    else if (x$format==2)
        out <- ppm2ratios(x,exterr=exterr)
    out
}
