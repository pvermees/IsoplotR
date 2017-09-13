get.PD.ratio <- function(tt,st,nuclide,exterr=TRUE){
    L <- lambda(nuclide)
    R <- exp(L[1]*tt)-1
    Jac <- matrix(0,1,2)
    E <- matrix(0,2,2)
    Jac[1,1] <- tt*exp(L[1]*tt)
    Jac[1,2] <- L[1]*exp(L[1]*tt)
    E[1,1] <- L[2]^2
    E[2,2] <- st^2
    sR <- sqrt(Jac %*% E %*% t(Jac))
    out <- c(R,sR)
}

get.PD.age <- function(DP,sDP,nuclide,exterr=TRUE){
    L <- lambda(nuclide)
    tt <- log(1 + DP)/L[1]
    E <- matrix(0,2,2)
    J <- matrix(0,1,2)
    E[1,1] <- sDP^2
    if (exterr) E[2,2] <- L[2]^2
    J[1,1] <- 1/(L[1]*(1 + DP))
    J[1,2] <- -tt/L[1]
    st <- sqrt(J %*% E %*% t(J))
    c(tt,st)
}

# i2i = isochron to intercept
PD.age <- function(x,nuclide,exterr=TRUE,i=NA,sigdig=NA,i2i=TRUE,...){
    ns <- length(x)
    dat <- data2york(x,exterr=exterr,common=FALSE)
    if (i2i){
        fit <- isochron(x,plot=FALSE,exterr=exterr)        
        dat[,'Y'] <- dat[,'Y'] - fit$a[1]
        if (exterr) dat[,'sY'] <- sqrt(dat[,'sY']^2 + fit$a[2]^2)
    }
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    E <- matrix(0,2,2)
    J <- matrix(0,1,2)
    for (j in 1:ns) {
        E[1,1] <- dat[j,'sX']^2
        E[2,2] <- dat[j,'sY']^2
        E[1,2] <- dat[j,'rXY']*dat[j,'sX']*dat[j,'sY']
        E[2,1] <- E[1,2]
        J[1,1] <- -dat[j,'Y']/dat[j,'X']^2
        J[1,2] <- 1/dat[j,'X']        
        DP <- dat[j,'Y']/dat[j,'X']
        sDP <- sqrt(J%*%E%*%t(J))
        tt <- get.PD.age(DP,sDP,nuclide,exterr=exterr)
        out[j,] <- roundit(tt[1],tt[2],sigdig=sigdig)
    }
    if (!is.na(i)) out <- out[i,]
    out
}

ppm2ratios <- function(x,...){ UseMethod("ppm2ratios",x) }
ppm2ratios.default <- function(x,...){
    stop('Method ppm2ratios not available for this class.')
}

length.PD <- function(x,...){ nrow(x$x) }
