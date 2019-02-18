#' @export
length.UPb  <- function(x){ nrow(x$x) }
#' @export
length.PbPb <- function(x){ nrow(x$x) }
#' @export
length.ArAr <- function(x){ nrow(x$x) }
#' @export
length.PD <- function(x){ nrow(x$x) }
#' @export
length.KCa <- function(x){ nrow(x$x) }
#' @export
length.ThU <- function(x){ nrow(x$x) }
#' @export
length.UThHe <- function(x){ nrow(x) }
#' @export
length.KDE <- function(x){ length(x$ages) }
#' @export
length.KDEs <- function(x){ length(x$kdes) }
#' @export
length.fissiontracks <- function(x){
    if (x$format==1) return(nrow(x$x))
    else return(length(x$Ns))
}

#' @export
subset.UPb  <- function(x,subset,...){ subset_helper(x,subset,...) }
#' @export
subset.PbPb <- function(x,subset,...){ subset_helper(x,subset,...) }
#' @export
subset.ArAr <- function(x,subset,...){ subset_helper(x,subset,...) }
#' @export
subset.KCa <- function(x,subset,...){ subset_helper(x,subset,...) }
#' @export
subset.PD <- function(x,subset,...){ subset_helper(x,subset,...) }
#' @export
subset.ThU <- function(x,subset,...){ subset_helper(x,subset,...) }
#' @export
subset.detritals <- function(x,subset,...){
    out <- x[subset]
    class(out) <- class(x)
    out
}
#' @export
subset.UThHe <- function(x,subset,...){
    out <- x[subset,]
    class(out) <- class(x)
    out
}
#' @export
subset.fissiontracks <- function(x,subset,...){
    if (x$format==1){
        out <- subset_helper(x,subset,...)
    } else {
        out <- x
        out$Ns <- x$Ns[subset]
        out$A <- x$A[subset]
        out$U <- x$U[subset]
        out$sU <- x$sU[subset]
    }
    out
}
subset_helper <- function(x,subset){
    out <- x
    out$x <- subset(x$x,subset=subset)
    out
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
        nsmall <- max(0,-(trunc(log10(min.err))-sigdig))
        out <- format(dat,digits=sigdig,nsmall=nsmall,trim=TRUE)
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
# returns standard error
errorprop1x2 <- function(J1,J2,E11,E22,E12){
    sqrt(E11*J1^2 + 2*E12*J1*J2 + E22*J2^2)
}
errorprop1x3 <- function(J1,J2,J3,E11,E22,E33,E12=0,E13=0,E23=0){
    sqrt(E11*J1^2 + E22*J2^2 + E33*J3^2 + 2*J1*J2*E12 + 2*J1*J3*E13 + 2*J2*J3*E23)
}

quotient <- function(X,sX,Y,sY,rXY){
    YX <- Y/X
    sYX <- errorprop1x2(J1=(-Y/X^2),J2=(1/X),
                        E11=(sX^2),E22=(sY^2),E12=(rXY*sX*sY))
    cbind(YX,sYX)
}

hasClass <- function(x,classname){
    classname %in% class(x)
}

# negative multivariate log likelihood to be fed into R's optim function
LL.norm <- function(x,covmat){
    (log(2*pi) + determinant(covmat,logarithmic=TRUE)$modulus) +
        get.concordia.SS(x,covmat)/2
}

set.ellipse.colours <- function(ns=1,levels=NA,col=c('yellow','red'),
                                hide=NULL,omit=NULL,omit.col=NA){
    nl <- length(levels)
    if (nl > 1){
        levels[c(hide,omit)] <- NA
        levels[!is.numeric(levels)] <- NA
    }
    out <- NULL
    if (all(is.na(levels)) ||
        (min(levels,na.rm=TRUE)==max(levels,na.rm=TRUE))){
        out <- rep(col[1],ns)
    } else if (nl<ns){
        out[1:nl] <- levels2colours(levels=levels,col=col)
    } else {
        out <- levels2colours(levels=levels,col=col)[1:ns]
    }
    out[omit] <- omit.col
    out
}

levels2colours <- function(levels=c(0,1),col=c('yellow','red')){
    m <- min(levels,na.rm=TRUE)
    M <- max(levels,na.rm=TRUE)
    fn <- grDevices::colorRamp(colors=col,alpha=TRUE)
    normalised.levels <- (levels-m)/(M-m)
    nan <- is.na(normalised.levels)
    normalised.levels[nan] <- 0
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
                      strip.width=0.02,clabel="",...){
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
    graphics::plot(rep(xe,length(z)),z,type='n',
                   axes=F,xlab=NA,ylab=NA,...)
    graphics::axis(side=4)
    mymtext(text=clabel,side=3,adj=1)
}

# mybg, mypch etc are the default values that will be used if bg, pch are not in ...
# this function can be used to override parameters that are in the ...
plot_points <- function(x,y,mybg='yellow',mypch=21,mycol='black',mycex=1.5,mypos=1,
                        show.numbers=FALSE,hide=NULL,omit=NULL,...){
    ns <- length(x)
    sn <- clear(1:ns,hide)
    pargs <- get_points_pars(...,mybg=mybg,mypch=mypch,mycol=mycol,mycex=mycex)
    targs <- get_text_pars(...,mycol=mycol,mycex=mycex,mypos=mypos)
    pargs$x <- x[sn]
    pargs$y <- y[sn]
    targs$x <- x[sn]
    targs$y <- y[sn]
    hascol <- !all(is.na(mybg))
    show.points <- (hascol | !show.numbers)
    if (show.points){
        do.call(graphics::points,pargs)
    }
    if (show.numbers & show.points){
        targs$labels <- sn
        do.call(graphics::text,targs)
    } else if (show.numbers & !show.points){
        tcol <- rep('black',ns)
        tcol[omit] <- 'grey'
        targs$col <- tcol[sn]
        targs$labels <- sn
        targs$pos <- NULL
        do.call(graphics::text,targs)
    }
}
get_points_pars <- function(...,mybg='yellow',mypch=21,mycol='black',mycex=1.5){
    ellipsis <- list(...)
    pars <- names(ellipsis)
    out <- ellipsis
    if ('pch' %in% pars)
        out$pch <- ellipsis$pch
    else
        out$pch <- mypch
    if ('col' %in% pars)
        out$col <- ellipsis$col
    else
        out$col <- mycol
    if ('bg' %in% pars)
        out$bg <- ellipsis$bg
    else
        out$bg <- mybg
    if ('cex' %in% pars)
        out$cex <- ellipsis$cex
    else
        out$cex <- mycex
    out
}
get_text_pars <- function(...,mycol='black',mycex=1.5,mypos=NULL){
    ellipsis <- list(...)
    pars <- names(ellipsis)
    out <- ellipsis
    if ('col' %in% pars)
        out$col <- ellipsis$col
    else
        out$col <- mycol
    if ('cex' %in% pars)
        out$cex <- ellipsis$cex
    else
        out$cex <- mycex
    if ('pos' %in% pars)
        out$pos <- ellipsis$pos
    else
        out$pos <- mypos
    out
}

nfact <- function(alpha){
    stats::qnorm(1-alpha/2)
}
tfact <- function(alpha,df){
    stats::qt(1-alpha/2,df=df)
}

mymtext <- function(text,line=0,...){
    graphics::mtext(text,line=line,cex=graphics::par('cex'),...)
}

# if doall==FALSE, only returns the lower right submatrix
blockinverse <- function(AA,BB,CC,DD,doall=FALSE){
    invAA <- solve(AA)
    invDCAB <- solve(DD-CC%*%invAA%*%BB)
    if (doall){
        ul <- invAA + invAA %*% BB %*% invDCAB %*% CC %*% invAA
        ur <- - invAA %*% BB %*% invDCAB
        ll <- - invDCAB %*% CC %*% invAA
        lr <- invDCAB
        out <- rbind(cbind(ul,ur),cbind(ll,lr))
    } else {
        out <- invDCAB
    }
    out
}

# Optimise with some fixed parameters 
# Like optim, but with option to fix some parameters.
# parms: Parameters to potentially optimize in fn
# fixed: A vector of TRUE/FALSE values indicating which parameters in
# parms to hold constant (not optimize). If TRUE, the corresponding
# parameter in fn() is fixed. Otherwise it's variable and optimised over.
# Originally written by Barry Rowlingson, modified by PV
optifix <- function(parms, fixed, fn, gr = NULL, ...,
                    method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
                    lower = -Inf, upper = Inf, control = list(), hessian = FALSE){ 
    force(fn)
    force(fixed) 
    .npar=length(parms) 
    .fixValues = parms[fixed]
    names(.fixValues)=names(parms)[fixed]
    .parStart = parms[!fixed]
    names(.parStart)=names(parms)[!fixed]
  
    .fn <- function(par,pnames=names(parms),...){
        .par = rep(NA,sum(!fixed))
        .par[!fixed] = par
        .par[fixed] = .fixValues
        names(.par)=pnames
        fn(.par,...)
    }

    if(!is.null(gr)){
        .gr <- function(par,pnames=names(parms),...){
            .gpar = rep(NA,sum(!fixed))
            .gpar[!fixed] = par
            .gpar[fixed] = .fixValues
            names(.par)=pnames
            gr(.gpar,...)[!fixed]
        }
    } else {
        .gr <- NULL
    }

    .opt = stats::optim(.parStart,.fn,.gr,...,method=method,
                        lower=lower,upper=upper,
                        control=control,hessian=hessian) 
    
    .opt$fullpars = rep(NA,sum(!fixed)) 
    .opt$fullpars[fixed]=.fixValues 
    .opt$fullpars[!fixed]=.opt$par 
    names(.opt$fullpars)=names(parms)
    .opt$fixed = fixed

    # remove fullpars (PV)
    .opt$par <- .opt$fullpars
    .opt$fullpars <- NULL
    
    return(.opt)
}

'%ni%' <- function(x,y)!('%in%'(x,y))

clear <- function(x,...){
    i <- unlist(list(...))
    if (hasClass(x,'matrix')) sn <- 1:nrow(x)
    else sn <- 1:length(x)
    if (is.numeric(i)) out <- subset(x,subset=sn%ni%i)
    else out <- x
    out
}
