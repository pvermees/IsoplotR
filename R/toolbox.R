numgrains <- function(x,...){ UseMethod("numgrains",x) }
numgrains.default <- function(x){ length(x[,1]) }
numgrains.UPb <- function(x){ length(x$x[,1]) }
numgrains.ArAr <- function(x){ length(x$x[,1]) }
numgrains.RbSr <- function(x){ length(x$x[,1]) }
numgrains.SmNd <- function(x){ length(x$x[,1]) }
numgrains.ReOs <- function(x){ length(x$x[,1]) }
numgrains.fissiontracks <- function(x){ length(x$x[,1]) }

select <- function(x,selection){
    out <- x
    i <- which(names(x$x) %in% selection)
    out$x <- x$x[i]
    out$covmat <- x$covmat[i,i]
    out
}

zip.matrix <- function(x){
    nc <- ncol(x)
    nr <- nrow(x)
    out <- rep(0,nr*nc)
    i <- seq(from=0,to=(nr-1)*nc,by=nc)
    for (j in 1:nc){
        out[i+j] <- x[,j]
    }
    out
}

unzip.vector <- function(x,nc=2){
    nr <- length(x)/nc
    out <- matrix(0,nr,nc)
    i <- seq(from=0,to=(nr-1)*nc,by=nc)
    for (j in 1:nc){
        out[,j] <- x[i+j]
    }
    out
}

roundit <- function(age,err,sigdig=2){
    out <- cbind(age,err)
    if (!is.na(sigdig)){
        out[,2] <- signif(err,sigdig)
        nd <- log10(trunc(abs(age)/err))+sigdig
        out[,1] <- signif(age,nd)
    }
    out
}

# set minimum and maximum values of a dataset
getmM <- function(x,from=NA,to=NA,log=FALSE){
    if (is.na(from)) { from <- min(x,na.rm=TRUE); getm <- TRUE }
    else { getm <- FALSE }
    if (is.na(to)) { to <- max(x,na.rm=TRUE); getM <- TRUE }
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

cor2cov <- function(sX,sY,rXY){
    covmat <- matrix(0,2,2)
    covmat[1,1] <- sX^2
    covmat[2,2] <- sY^2
    covmat[1,2] <- rXY*sX*sY
    covmat[2,1] <- covmat[1,2]
    covmat
}

get.cov.xzyz <- function(xz,err.xz,yz,err.yz,err.xy){
    xy <- xz/yz
    0.5*xz*yz*((err.xz/xz)^2 + (err.yz/yz)^2 - (err.xy/xy)^2)
}
get.cor.xzyz <- function(xz,err.xz,yz,err.yz,err.xy){
    get.cov.xzyz(xz,err.xz,yz,err.yz,err.xy)/(err.xz*err.yz)
}

get.cov.zxzy <- function(zx,err.zx,zy,err.zy,err.xy){
    xy <- zy/zx
    0.5*zx*zy*((err.zy/zy)^2 + (err.zx/zx)^2 - (err.xy/xy)^2)
}
get.cor.zxzy <- function(zx,err.zx,zy,err.zy,err.xy){
    get.cov.zxzy(zx,err.zx,zy,err.zy,err.xy)/(err.zx*err.zy)
}

get.cov.xzzy <- function(xz,err.xz,zy,err.zy,err.xy){
    xy <- xz*zy
    0.5*xz*zy*((err.xy/xy)^2 - (err.xz/xz)^2 - (err.zy/zy)^2)
}
get.cor.xzzy <- function(xz,err.xz,zy,err.zy,err.xy){
    get.cov.xzzy(xz,err.xz,zy,err.zy,err.xy)/(err.xz*err.zy)
}

# simultaneously performs error propagation for multiple samples
errorprop <- function(J11,J12,J21,J22,E11,E12,E22){
    out <- matrix(0,length(J11),3)
    colnames(out) <- c('varX','varY','cov')
    out[,'varX'] <- J11*J11*E11 + J11*J12*E12 + J11*J12*E12 + J12*J12*E22
    out[,'varY'] <- E11*J21*J21 + J21*J22*E12 + J21*J22*E12 + J22*J22*E22
    out[,'cov'] <- J11*J21*E11 + J12*J21*E12 + J11*J22*E12 + J12*J22*E22
    out
}

hasClass <- function(x,classname){
    classname %in% class(x)
}

get.covmat <- function(x,...){ UseMethod("get.covmat",x) }
get.covmat.default <- function(x,i,...){ stop('Invalid input into covmat() function') }

get.selection <- function(x,...){ UseMethod("get.selection",x) }
get.selection.default <- function(x,...){ x }

# negative multivariate log likelihood to be fed into R's optim function
LL.norm <- function(x,covmat){
    log(2*pi) + 0.5*determinant(covmat,logarithmic=TRUE)$modulus +
                                0.5*get.concordia.SS(x,covmat)
}
