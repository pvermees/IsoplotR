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
    fit <- optim(init,LL.irr.generic,XYZ=XYZ)
}
irr.PD <- function(x,alpha=0.05,inverse=FALSE,...){
    
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
    Ei <- matrix(0,4,4)
    for (i in 1:ns){
        xyz <- J %*% c(u[i],v[i])
        D <- XYZ[[i]]$XYZ - xyz
        out <- out + t(D) %*% XYZ[[i]]$O %*% D
    }
    -out/2
}
LL.irr.a <- function(abg,XYZ,w=0){
}
LL.irr.b <- function(abg,XYZ,w=0){
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
