length.UPb  <- function(x){ nrow(x$x) }
length.PbPb <- function(x){ nrow(x$x) }
length.ArAr <- function(x){ nrow(x$x) }
length.RbSr <- function(x){ nrow(x$x) }
length.SmNd <- function(x){ nrow(x$x) }
length.ReOs <- function(x){ nrow(x$x) }
length.LuHf <- function(x){ nrow(x$x) }
length.ThU <- function(x){ nrow(x$x) }
length.fissiontracks <- function(x){ length(x$Ns) }
length.UThHe <- function(x){ nrow(x) }

roundit <- function(age,err,sigdig=2){
    if (length(age)==1){
        min.err <- min(err,na.rm=TRUE)
        dat <- c(age,err)
    } else {
        min.err <- err
        dat <- cbind(age,err)
    }
    if (is.na(sigdig)) {
        out <- dat
    } else {
        nd <- log10(trunc(abs(dat)/min.err))+sigdig
        out <- signif(dat,nd)
    }
    out
}

# count the number of TRUEs in x
count <- function(x){ length(which(x)) }

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

cor2cov2 <- function(sX,sY,rXY){
    covmat <- matrix(0,2,2)
    covmat[1,1] <- sX^2
    covmat[2,2] <- sY^2
    covmat[1,2] <- rXY*sX*sY
    covmat[2,1] <- covmat[1,2]
    covmat
}
cor2cov3 <- function(sX,sY,sZ,rXY,rXZ,rYZ){
    covmat <- matrix(0,3,3)
    covmat[1,1] <- sX^2
    covmat[2,2] <- sY^2
    covmat[3,3] <- sZ^2
    covmat[1,2] <- rXY*sX*sY
    covmat[1,3] <- rXZ*sX*sZ
    covmat[2,3] <- rYZ*sY*sZ
    covmat[2,1] <- covmat[1,2]
    covmat[3,1] <- covmat[1,3]
    covmat[3,2] <- covmat[2,3]
    covmat
}

get.cov.div <- function(A,err.A,B,err.B,AB,err.AB){
    0.5*A*B*((err.A/A)^2+(err.B/B)^2-(err.AB/AB)^2)
}
get.cor.div <- function(A,err.A,B,err.B,AB,err.AB){
    get.cov.div(A,err.A,B,err.B,AB,err.AB)/(err.A*err.B)
}
get.cov.mult <- function(A,err.A,B,err.B,AB,err.AB){
    0.5*A*B*((err.AB/AB)^2 - (err.A/A)^2 - (err.B/B)^2)
}
get.cor.mult <- function(A,err.A,B,err.B,AB,err.AB){
    get.cov.mult(A,err.A,B,err.B,AB,err.AB)/(err.A*err.B)
}

# simultaneously performs error propagation for multiple samples
errorprop <- function(J11,J12,J21,J22,E11,E22,E12){
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

# negative multivariate log likelihood to be fed into R's optim function
LL.norm <- function(x,covmat){
    (log(2*pi) + determinant(covmat,logarithmic=TRUE)$modulus) +
        get.concordia.SS(x,covmat)/2
}

set.ellipse.colours <- function(ns=1,levels=NA,col=c('yellow','red')){
    if (any(!is.numeric(levels)) | any(is.na(levels))) levels <- NA
    nl <- length(levels)
    out <- NULL
    if (all(is.na(levels)) | (max(levels)==min(levels))){
        out <- rep(col[1],ns)
    } else if (nl<ns){
        out[1:nl] <- levels2colours(levels=levels,col=col)
    } else {
        out <- levels2colours(levels=levels,col=col)[1:ns]
    }
    out
}

levels2colours <- function(levels=c(0,1),col=c('yellow','red')){
    fn <- grDevices::colorRamp(colors=col,alpha=TRUE)
    normalised.levels <- (levels-min(levels))/(max(levels)-min(levels))
    col.matrix <- fn(normalised.levels)/255
    red <- col.matrix[,1]
    green <- col.matrix[,2]
    blue <- col.matrix[,3]
    alpha <- col.matrix[,4]
    grDevices::rgb(red,green,blue,alpha)
}

validLevels <- function(levels){
    (all(is.numeric(levels)) & !any(is.na(levels)))
}

colourbar <- function(z=c(0,1),col=c("#00FF0080","#FF000080"),
                      strip.width=0.02){
    if (!validLevels(z)) return()
    ucoord <- graphics::par()$usr
    plotwidth <- (ucoord[2]-ucoord[1])
    plotheight <- (ucoord[4]-ucoord[3])
    xe <- ucoord[2]
    xb <- xe - strip.width*plotwidth
    yb <- ucoord[3]
    ye <- ucoord[4]
    ndiv <- 50 # number of divisions
    dx <- (xe-xb)/ndiv
    dy <- (ye-yb)/ndiv
    zz <- seq(from=min(z),to=max(z),length.out=ndiv)
    cc <- levels2colours(levels=zz,col=col)
    for (i in 1:ndiv){
        graphics::rect(xb,yb+(i-1)*dy,xe,yb+i*dy,col=cc[i],border=NA)
    }
    graphics::rect(xb,yb,xe,ye)
    graphics::par(new=T)
    graphics::plot(rep(xe,length(z)),z,type='n',axes=F,xlab=NA,ylab=NA)
    graphics::axis(side=4)
}
