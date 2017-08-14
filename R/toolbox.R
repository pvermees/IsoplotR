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
    log(2*pi) + 0.5*determinant(covmat,logarithmic=TRUE)$modulus +
                                0.5*get.concordia.SS(x,covmat)
}
