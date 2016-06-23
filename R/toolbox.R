roundit <- function(age,err){
    out <- list()
    out$err <- signif(err,2)
    nd <- log10(trunc(abs(age)/err))+2
    out$x <- signif(age,nd)
    out
}

# set minimum and maximum values of a dataset
getmM <- function(x,from=NA,to=NA,log=FALSE){
    if (is.na(from)) { from <- min(x); getm <- TRUE }
    else { getm <- FALSE }
    if (is.na(to)) { to <- max(x); getM <- TRUE }
    else { getM <- FALSE }
    if (getm) {
        if (log) { from <- from/2 }
        else {
            if (2*from-to<0) {from <- 0}
            else {from <- from-(to-from)/10}
        }
    }
    if (getM) {
        if (log) { to <- 2*to }
        else { to <- to+(to-from)/10 }
    }
    list(m=from,M=to)
}

emptyplot <- function(){
    graphics::plot(c(0,1),c(0,1),type='n',axes=FALSE,xlab="",ylab="")
}

cov2cor <- function(covmat){
    covmat[1,2]/sqrt(covmat[1,1]*covmat[2,2])
}

cor2cov <- function(sX,sY,rXY){
    covmat <- matrix(0,2,2)
    covmat[1,1] <- sX^2
    covmat[2,2] <- sY^2
    covmat[1,2] <- rXY*sX*sY
    covmat[2,1] <- covmat[1,2]
    covmat
}

get.covmat.xzyz <- function(xz,err.xz,yz,err.yz,err.xy){
    xy <- xz/yz
    0.5*xz*yz*((err.xz/xz)^2 + (err.yz/yz)^2 - (err.xy/xy)^2)
}

get.covmat.zxzy <- function(zx,err.zx,zy,err.zy,err.xy){
    xy <- zy/zx
    0.5*zx*zy*((err.zy/zy)^2+(err.zx/zx)^2+(err.xy/xy)^2)
}

get.covmat.xzzy <- function(xz,err.xz,zy,err.zy,err.xy){
    xy <- xz*zy
    0.5*xz*zy*((err.xy/xy)^2-(err.xz/xz)^2-(err.zy/zy)^2)
}
