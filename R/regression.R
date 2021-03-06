regression <- function(xyz,model=1,type='york',omit=NULL){
    xyz2calc <- clear(xyz,omit)
    if (model==1) out <- model1regression(xyz2calc,type=type)
    else if (model==2) out <- model2regression(xyz2calc,type=type)
    else if (model==3) out <- model3regression(xyz2calc,type=type)
    else stop('invalid regression model')
    out$xyz <- xyz
    out$model <- model
    out$n <- nrow(xyz2calc)
    out$omit <- omit
    out
}

model1regression <- function(xyz,type='york'){
    if (identical(type,'york')){
        out <- york(xyz)
    } else if (identical(type,'titterington')){
        out <- titterington(xyz)
    } else {
        stop('invalid output type for model 1 regression')
    }
    out
}

model2regression <- function(xyz,type='york'){
    if (identical(type,'york')){
        out <- list()
        fit <- stats::lm(xyz[,'Y'] ~ xyz[,'X'])
        E <- stats::vcov(fit)
        out$df <- fit$df.residual
        out$a <- c(stats::coef(fit)[1],sqrt(E[1,1]))
        out$b <- c(stats::coef(fit)[2],sqrt(E[2,2]))
        names(out$a) <- c('a','s[a]')
        names(out$b) <- c('b','s[b]')
        out$cov.ab <- E[1,2]
    } else if (identical(type,'titterington')){
        out <- list()
        fit <- stats::lm(xyz[,c('Y','Z')] ~ xyz[,'X'])
        out$df <- fit$df.residual
        out$par <- c(stats::coef(fit))
        out$cov <- stats::vcov(fit)
        parnames <- c('a','b','A','B')
        names(out$par) <- parnames
        colnames(out$cov) <- parnames
        rownames(out$cov) <- parnames
    } else {
        stop('invalid output type for model 2 regression')
    }
    out
}

model3regression <- function(xyz,type='york'){
    out <- list()
    if (identical(type,'york')){
        out$w <- stats::optimize(LL.york,
                                 interval=c(0,stats::sd(xyz[,'Y'])),
                                 xy=xyz,maximum=TRUE)$maximum
        dd <- augment_york_errors(xyz,out$w)
        out <- c(out,york(dd))
    } else if (identical(type,'titterington')){
        out$w <- stats::optimize(LL.titterington,
                                 interval=c(0,stats::sd(xyz[,'Y'])),
                                 xyz=xyz,maximum=TRUE)$maximum
        dd <- augment_titterington_errors(xyz,out$w)
        out <- c(out,titterington(dd))
    } else {
        stop('invalid output type for model 3 regression')
    }
    out
}

LL.isochron <- function(w,xyz,type='york'){
    if (identical(type,'york'))
        out <- LL.york(w,xyz)
    else if (identical(type,'titterington'))
        out <- LL.titterington(w,xyz)
    out
}
LL.york <- function(w,xy){
    D <- augment_york_errors(xy,w)
    X <- matrix(0,1,2)
    fit <- york(D)
    P <- get.york.xy(D,fit$a[1],fit$b[1])
    x1 <- D[,'X']
    x2 <- D[,'Y']
    s1 <- D[,'sX']
    s2 <- D[,'sY']
    rho <- D[,'rXY']
    m1 <- P[,1]
    m2 <- P[,2]
    z <- ((x1-m1)/s1)^2 + ((x2-m2)/s2)^2 - 2*rho*(x1-m1)*(x2-m2)/(s1*s2)
    LL <- - log(s1) - log(s2) - 0.5*log(1-rho^2) - 0.5*z/(1-rho^2)
    sum(LL)
}
LL.titterington <- function(w,xyz){
    out <- 0
    D <- augment_titterington_errors(xyz,w)
    fit <- titterington(D)
    dat <- matrix2covlist(D)
    X <- matrix(0,1,3)
    ns <- nrow(D)
    for (i in 1:ns){
        XYZ <- dat$XYZ[[i]]
        O <- dat$omega[[i]]
        E <- solve(O)
        a <- fit$par[1]
        b <- fit$par[2]
        A <- fit$par[3]
        B <- fit$par[4]
        abg <- alpha.beta.gamma(a,b,A,B,XYZ,O)
        X[1,1] <- -abg[2]/abg[1]
        X[1,2] <- XYZ[2] - a - b*XYZ[1]
        X[1,3] <- XYZ[3] - A - B*XYZ[1]
        SS <- X %*% O %*% t(X)
        detE <- determinant(E,logarithm=TRUE)$modulus
        out <- out - 1.5*log(2*pi) - 0.5*detE - 0.5*SS
    }
    out
}
LL.ludwig <- function(w,xyz){
}

augment_york_errors <- function(xy,w){
    out <- xy
    out[,'sY'] <- sqrt(xy[,'sY']^2 + w^2)
    out[,'rXY'] <- xy[,'rXY']*xy[,'sY']/out[,'sY']
    out
}
augment_titterington_errors <- function(xyz,w){
    out <- xyz
    out[,'sY'] <- sqrt(xyz[,'sY']^2 + w^2)
    out[,'rXY'] <- xyz[,'rXY']*xyz[,'sY']/out[,'sY']
    out[,'rYZ'] <- xyz[,'rYZ']*xyz[,'sY']/out[,'sY']
    out
}
