length.UPb  <- function(x){ nrow(x$x) }
length.PbPb <- function(x){ nrow(x$x) }
length.ArAr <- function(x){ nrow(x$x) }
length.RbSr <- function(x){ nrow(x$x) }
length.SmNd <- function(x){ nrow(x$x) }
length.ReOs <- function(x){ nrow(x$x) }
length.LuHf <- function(x){ nrow(x$x) }
length.ThU <- function(x){ nrow(x$x) }
length.UThHe <- function(x){ nrow(x) }
length.fissiontracks <- function(x){
    if (x$format==1) return(nrow(x$x))
    else return(length(x$Ns))
}

roundit <- function(age,err,sigdig=2){
    if (length(age)==1) dat <- c(age,err)
    else dat <- cbind(age,err)
    min.err <- min(abs(dat),na.rm=TRUE)
    if (is.na(sigdig)) {
        out <- dat
    } else if (any(dat<=0, na.rm=TRUE)){
        out <- signif(dat,sigdig)
    } else {
        nd <- ceiling(log10(abs(dat)))-ceiling(log10(min.err))+sigdig
        rounded <- signif(dat,nd)
        decimals <- trunc(log10(min.err))-sigdig
        if (decimals<1) nsmall <- abs(decimals)
        else nsmall <- 0
        out <- format(rounded,digits=max(0,nsmall-1),nsmall=nsmall,trim=TRUE)
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

# Implements Equations 6 & 7 of Ludwig (1998)
# x is an [n x 5] matrix with columns X, sX, Y, sY, rhoXY
wtdmean2D <- function(x){
    ns <- nrow(x)
    X <- x[,1]
    Y <- x[,3]
    O11 <- rep(0,ns)
    O22 <- rep(0,ns)
    O12 <- rep(0,ns)
    for (i in 1:ns){
        E <- cor2cov2(x[i,2],x[i,4],x[i,5])
        O <- solve(E)
        O11[i] <- O[1,1]
        O22[i] <- O[2,2]
        O12[i] <- O[1,2]
    }
    numx <- sum(O22)*sum(X*O11+Y*O12)-sum(O12)*sum(Y*O22+X*O12)
    den <- sum(O11)*sum(O22)-sum(O12)^2
    numy <- sum(O11)*sum(Y*O22+X*O12)-sum(O12)*sum(X*O11+Y*O12)
    xbar <- numx/den
    ybar <- numy/den
    out <- list()
    out$x <- c(xbar,ybar)
    Obar <- matrix(0,2,2)
    Obar[1,1] <- sum(O11)
    Obar[2,2] <- sum(O22)
    Obar[1,2] <- sum(O12)
    Obar[2,1] <- Obar[1,2]
    out$cov <- solve(Obar)
    out
}
# generalisation of wtdmean2D to 3 dimensions
wtdmean3D <- function(x){
    ns <- nrow(x)
    X <- x[,1]
    Y <- x[,3]
    Z <- x[,5]
    O11 <- rep(0,ns)
    O22 <- rep(0,ns)
    O33 <- rep(0,ns)
    O12 <- rep(0,ns)
    O13 <- rep(0,ns)
    O23 <- rep(0,ns)
    for (i in 1:ns){
        E <- cor2cov3(x[i,2],x[i,4],x[i,6],x[i,7],x[i,8],x[i,9])
        O <- solve(E)
        O11[i] <- O[1,1]
        O22[i] <- O[2,2]
        O33[i] <- O[3,3]
        O12[i] <- O[1,2]
        O13[i] <- O[1,3]
        O23[i] <- O[2,3]
    }
    AA <- sum(X*O11+Y*O12+Z*O13)
    BB <- sum(X*O12+Y*O22+Z*O23)
    CC <- sum(X*O13+Y*O23+Z*O33)
    DD <- sum(O22)*sum(O11)-sum(O12)^2
    EE <- BB*sum(O11)-AA*sum(O12)
    FF <- sum(O13)*sum(O12)-sum(O23)*sum(O11)
    GG <- EE/DD
    HH <- FF/DD
    II <- (AA-GG*sum(O12))/sum(O11)
    JJ <- (HH*sum(O12)+sum(O13))/sum(O11)
    KK <- CC-II*sum(O13)-GG*sum(O23)
    LL <- HH*sum(O23)-JJ*sum(O13)+sum(O33)
    xbar <- II-JJ*KK/LL
    ybar <- GG+HH*KK/LL
    zbar <- KK/LL
    out <- list()
    out$x <- c(xbar,ybar,zbar)
    Obar <- matrix(0,3,3)
    Obar[1,1] <- sum(O11)
    Obar[2,2] <- sum(O22)
    Obar[3,3] <- sum(O33)
    Obar[1,2] <- sum(O12)
    Obar[1,3] <- sum(O13)
    Obar[2,3] <- sum(O23)
    Obar[2,1] <- Obar[1,2]
    Obar[3,1] <- Obar[1,3]
    Obar[3,2] <- Obar[2,3]
    out$cov <- solve(Obar)
    out
}

# simultaneously performs error propagation for multiple samples
errorprop <- function(J11,J12,J21,J22,E11,E22,E12){
    out <- matrix(0,length(J11),3)
    colnames(out) <- c('varX','varY','cov')
    out[,'varX'] <- J11*J11*E11 + J11*J12*E12 + J11*J12*E12 + J12*J12*E22
    out[,'varY'] <- J21*J21*E11 + J21*J22*E12 + J21*J22*E12 + J22*J22*E22
    out[,'cov'] <- J11*J21*E11 + J12*J21*E12 + J11*J22*E12 + J12*J22*E22
    out
}
errorprop1x2 <- function(J1,J2,E11,E22,E12){
    sqrt(E11*J1^2 + 2*E12*J1*J2 + E22*J2^2)
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
                      strip.width=0.02,clabel=""){
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
    graphics::mtext(text=clabel,side=3,adj=1)
}

plot_points <- function(x,y,bg='white',show.numbers=FALSE,...){
    if ('pch' %in% graphics::par(list(...))) pch <- graphics::par()$pch
    else pch <- 21
    if ('cex' %in% graphics::par(list(...))) pch <- graphics::par()$cex
    else cex <- 1.5
    if (show.numbers) graphics::text(x,y,1:length(x),cex=cex,...)
    else graphics::points(x,y,bg=bg,pch=pch,cex=cex,...)
}

nfact <- function(alpha){
    stats::qnorm(1-alpha/2)
}
tfact <- function(alpha,df){
    stats::qt(1-alpha/2,df=df)
}
