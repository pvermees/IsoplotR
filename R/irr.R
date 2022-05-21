irr <- function(x,...){ UseMethod("irr",x) }
irr.default <- function(x,alpha=0.05,...){
    colnames(x) <- c('X','sX','Y','sY','rXY')
    yfit <- york(x)
    a <- log(yfit$a[1])
    b <- yfit$b[1]
    v <- log(get.york.xy(x,a=yfit$a[1],b=b)[,2])
    g <- log((v-a)/b)
    init <- c(a,b,g)
    XYZ <- york2clr(x)
    fit <- stats::optim(init,LL.irr.generic,XYZ=XYZ)
}
irr.RbSr <- function(x,alpha=0.05,model=1,inverse=FALSE,...){
    yfit <- isochron(x,plot=FALSE)
    y <- yfit$xyz
    ap <- log(yfit$a[1])
    bp <- log(yfit$age[1])
    ydat <- data2york(x)
    v <- log(get.york.xy(ydat,a=yfit$a[1],b=yfit$b[1])[,2])
    g <- log((v-ap)/yfit$b[1])
    init <- c(ap,bp,g)
    XYZ <- york2clr(ydat)
    fit <- stats::optim(init,LL.irr.a,XYZ=XYZ,
                        lam=lambda('Rb87'),inverse=inverse)
    fit
}
LL.irr.generic <- function(abg,XYZ){
    ns <- length(XYZ)
    a <- abg[1]
    b <- abg[2]
    g <- abg[-c(1,2)]
    v <- a + b*exp(g)
    u <- log((exp(v)-exp(a))/b)
    J <- cbind(c(2,-1,-1),c(-1,2,-1))/3
    out <- 0
    for (i in 1:ns){
        xyz <- J %*% c(u[i],v[i])
        D <- XYZ[[i]]$XYZ - xyz
        out <- out + t(D) %*% XYZ[[i]]$O %*% D
    }
    -out/2
}
LL.irr.a <- function(abgw,XYZ,lam,inverse=FALSE){
    wgba <- rev(abgw)
    abg <- rev(wgba[-1])
    LL.irr.a.w(abg,XYZ=XYZ,w=wgba[1],lam=lam,inverse=inverse)
}
LL.irr.a.w <- function(abg,XYZ,lam,w=-Inf,inverse=FALSE){
    if (inverse){
        a <- -abg[1]
        b <- -exp(abg[1])*(exp(lam[1]*exp(abg[2]))-1)
    } else {
        a <- abg[1]
        b <- exp(lam[1]*exp(abg[2]))-1
    }
    g <- abg[-c(1,2)]
    v <- a + b*exp(g)
    u <- log((exp(v)-exp(a))/b)
    J <- cbind(c(2,-1,-1),c(-1,2,-1))/3
    out <- 0
    if (is.finite(w)){
        Ei <- matrix(0,4,4)
        Ei[4,4] <- exp(w)
        du.da <- (1-exp(a))/(v-exp(a))
        da.dap <- ifelse(inverse,-1,1)
        Ji <- matrix(0,3,4)
        Ji[1:3,1:3] <- diag(3)
    }
    for (i in seq_along(XYZ)){
        if (is.finite(w)){
            Ei[1:3,1:3] <- XYZ[[i]]$E
            dxyz.da <- J %*% c(du.da[i],1)
            Ji[,4] <- dxyz.da*da.dap
            Eip <- Ji %*% Ei %*% t(Ji)
            O <- MASS::ginv(Eip)
        } else {
            O <- XYZ[[i]]$O
        }
        xyz <- J %*% c(u[i],v[i])
        D <- XYZ[[i]]$XYZ - xyz
        out <- out + t(D) %*% O %*% D
    }
    -out/2
}
LL.irr.b <- function(abgw,XYZ){
}

york2clr <- function(x){
    ns <- nrow(x)
    U <- log(x[,1])
    V <- log(x[,3])
    sU <- x[,2]/x[,1]
    sV <- x[,4]/x[,3]
    sUV <- x[,5]*sU*sV
    E <- matrix(0,2,2)
    O <- matrix(0,2,2)
    J <- cbind(c(2,-1,-1),c(-1,2,-1))/3
    out <- list()
    for (i in 1:ns){
        UV <- c(U[i],V[i])
        E[1,1] <- sU[i]^2
        E[2,2] <- sV[i]^2
        E[1,2] <- E[2,1] <- sUV[i]
        out[[i]] <- list()
        out[[i]]$XYZ <- J %*% UV
        out[[i]]$E <- J %*% E %*% t(J)
        out[[i]]$O <- MASS::ginv(out[[i]]$E)
    }
    out
}

alr2clr <- function(u,v){

}
