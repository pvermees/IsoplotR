get.PD.ratio <- function(tt,st,nuclide,exterr=TRUE,bratio=1){
    L <- lambda(nuclide)
    R <- bratio*(exp(L[1]*tt)-1)
    Jac <- matrix(0,1,2)
    E <- matrix(0,2,2)
    Jac[1,1] <- tt*exp(L[1]*tt)
    Jac[1,2] <- L[1]*exp(L[1]*tt)
    E[1,1] <- L[2]^2
    E[2,2] <- st^2
    sR <- bratio*sqrt(Jac %*% E %*% t(Jac))
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
    if (ns<2) i2i <- FALSE
    dat <- data2york(x,exterr=exterr)
    if (i2i){
        fit <- isochron(x,plot=FALSE,exterr=exterr)
        initial <- fit$a
    } else {
        initial <- get.initial(x)
    }
    dat[,'Y'] <- dat[,'Y'] - initial[1]
    if (exterr) dat[,'sY'] <- sqrt(dat[,'sY']^2 + initial[2]^2)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    DP <- quotient(dat[,'X'],dat[,'sX'],dat[,'Y'],dat[,'sY'],dat[,'rXY'])
    tt <- get.PD.age(subset(DP,select=1),subset(DP,select=2),
                     nuclide,exterr=exterr,bratio=bratio)
    out <- roundit(subset(tt,select=1),subset(tt,select=2),sigdig=sigdig)
    if (!is.na(i)) out <- out[i,]
    out
}

get.initial <- function(x){
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
    out
}

ppm2ratios <- function(x,...){ UseMethod("ppm2ratios",x) }
ppm2ratios.default <- function(x,...){
    stop('Method ppm2ratios not available for this class.')
}

length.PD <- function(x,...){ nrow(x$x) }
