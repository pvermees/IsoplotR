roundit <- function(age,err,sigdig=2,oerr=3,text=FALSE,maxprecision=8){
    if (oerr>3){
        out <- roundit(age,err*age/100,sigdig=sigdig)
        out[-1] <- signif(err,sigdig)
    } else if (missing(err)){
        if (is.na(sigdig)) out <- age
        else out <- signif(age,digits=sigdig)
        nc <- 1
    } else if (any(age<0)){
        s <- sign(age)
        out <- roundit(age=abs(age),err=err,sigdig=sigdig)
        if (is.matrix(out)){
            out[,1] <- s*out[,1]
            nc <- ncol(out)
        } else {
            out[1] <- s*out[1]
            nc <- length(out)
        }
    } else {
        impossible <- (err>0 & log10(age/err)>maxprecision)
        err[impossible] <- 0 # impossibly good precision
        if (length(age)==1){
            dat <- c(age,err)
            nc <- length(dat)
        } else {
            dat <- cbind(age,err)
            nc <- ncol(dat)
        }
        min.err <- min(abs(dat),na.rm=TRUE)
        if (is.na(sigdig) | min.err==0) {
            out <- dat
        } else {
            scientific <- abs(log10(min.err))>10
            out <- format(dat,digits=sigdig,trim=TRUE,scientific=scientific)
        }
    }
    if (!text){
        suppressWarnings( # suppress NA warnings
            out <- matrix(as.numeric(out),ncol=nc)
        )
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
    dx <- to-from
    if (getm) {
        if (log) { from <- max(from/2,from-dx) }
        else {
            if ((2*from-to)<0) {from <- 0}
            else {from <- from-dx/10}
        }
    }
    if (getM) {
        if (log) { to <- min(2*to,to+dx) }
        else { to <- to+dx/10 }
    }
    list(m=from,M=to)
}

cor2cov2 <- function(sX,sY,rXY){
    covmat <- matrix(0,2,2)
    covmat[1,1] <- sX^2
    covmat[2,2] <- sY^2
    covmat[1,2] <- covmat[2,1] <- ifelse(rXY==0,0,rXY*sX*sY)
    covmat
}
cor2cov3 <- function(sX,sY,sZ,rXY,rXZ,rYZ){
    covmat <- matrix(0,3,3)
    covmat[1,1] <- sX^2
    covmat[2,2] <- sY^2
    covmat[3,3] <- sZ^2
    covmat[1,2] <- covmat[2,1] <- ifelse(rXY==0,0,rXY*sX*sY)
    covmat[1,3] <- covmat[3,1] <- ifelse(rXY==0,0,rXZ*sX*sZ)
    covmat[2,3] <- covmat[3,2] <- ifelse(rYZ==0,0,rYZ*sY*sZ)
    covmat
}
cor2cov4 <- function(sW,sX,sY,sZ,rWX,rWY,rWZ,rXY,rXZ,rYZ){
    covmat <- matrix(0,4,4)
    covmat[1,1] <- sW^2
    covmat[2,2] <- sX^2
    covmat[3,3] <- sY^2
    covmat[4,4] <- sZ^2
    covmat[1,2] <- covmat[2,1] <- rWX*sW*sX
    covmat[1,3] <- covmat[3,1] <- rWY*sW*sY
    covmat[1,4] <- covmat[4,1] <- rWZ*sW*sZ
    covmat[2,3] <- covmat[3,2] <- rXY*sX*sY
    covmat[2,4] <- covmat[4,2] <- rXZ*sX*sZ
    covmat[3,4] <- covmat[4,3] <- rYZ*sY*sZ
    covmat
}

get_cov_div <- function(A,err.A,B,err.B,AB,err.AB){
    0.5*A*B*((err.A/A)^2+(err.B/B)^2-(err.AB/AB)^2)
}
get_cor_div <- function(A,err.A,B,err.B,AB,err.AB){
    get_cov_div(A,err.A,B,err.B,AB,err.AB)/(err.A*err.B)
}
get_cov_mult <- function(A,err.A,B,err.B,AB,err.AB){
    0.5*A*B*((err.AB/AB)^2 - (err.A/A)^2 - (err.B/B)^2)
}
get_cor_mult <- function(A,err.A,B,err.B,AB,err.AB){
    get_cov_mult(A,err.A,B,err.B,AB,err.AB)/(err.A*err.B)
}

# Implements Equations 6 & 7 of Ludwig (1998)
# x is an [n x 5] matrix with columns X, sX, Y, sY, rXY
wtdmean2D <- function(x){
    ns <- nrow(x)
    X <- x[,1]
    Y <- x[,3]
    O11 <- rep(0,ns)
    O22 <- rep(0,ns)
    O12 <- rep(0,ns)
    for (i in 1:ns){
        E <- cor2cov2(x[i,2],x[i,4],x[i,5])
        DET <- E[1,1]*E[2,2]-E[1,2]*E[2,1]
        O11[i] <- E[2,2]/DET
        O22[i] <- E[1,1]/DET
        O12[i] <- -E[1,2]/DET
    }
    numx <- sum(O22)*sum(X*O11+Y*O12)-sum(O12)*sum(Y*O22+X*O12)
    numy <- sum(O11)*sum(Y*O22+X*O12)-sum(O12)*sum(X*O11+Y*O12)
    den <- sum(O11)*sum(O22)-sum(O12)^2
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

fixroundingerr <- function(v){
    v[(v<0) & (v>(-1e-10))] <- 0
    v
}

# simultaneously performs error propagation for multiple samples
errorprop <- function(J11,J12,J21,J22,E11,E22,E12=0){
    out <- matrix(0,length(J11),3)
    colnames(out) <- c('varX','varY','cov')
    out[,'varX'] <- fixroundingerr(J11*J11*E11 + J11*J12*E12 + J11*J12*E12 + J12*J12*E22)
    out[,'varY'] <- fixroundingerr(J21*J21*E11 + J21*J22*E12 + J21*J22*E12 + J22*J22*E22)
    out[,'cov'] <- J11*J21*E11 + J12*J21*E12 + J11*J22*E12 + J12*J22*E22
    out
}
errorprop1x2 <- function(J1,J2,E11,E22,E12=0){
    fixroundingerr(E11*J1^2 + 2*E12*J1*J2 + E22*J2^2)
}
errorprop1x3 <- function(J1,J2,J3,E11,E22,E33,E12=0,E13=0,E23=0){
    fixroundingerr(E11*J1^2 + E22*J2^2 + E33*J3^2 +
                   2*J1*J2*E12 + 2*J1*J3*E13 + 2*J2*J3*E23)
}

quotient <- function(X,sX,Y,sY,rXY=NULL,sXY=NULL){
    YX <- Y/X
    if (is.null(sXY)){
        if (is.null(rXY)) E12 <- 0
        else E12 <- rXY*sX*sY
    } else {
        E12 <- sXY
    }
    v <- errorprop1x2(J1=-Y/X^2,J2=1/X,E11=sX^2,E22=sY^2,E12=E12)
    cbind(YX,sqrt(v))
}

# negative multivariate log likelihood to be fed into R's optim function
LL_norm <- function(x,covmat){
    tryCatch({
        (log(2*pi) + determinant(covmat,logarithmic=TRUE)$modulus
            + stats::mahalanobis(x,center=FALSE,cov=covmat))/2
    },error=function(e) Inf)
}

set_ellipse_colours <- function(ns=1,levels=NA,col=c('yellow','red'),
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
# To be removed. Kept for backwards compatibility in provenance 4.2
set.ellipse.colours <- function(ns=1,levels=NA,col=c('yellow','red'),
                                hide=NULL,omit=NULL,omit.col=NA){
    set_ellipse_colours(ns=ns,levels=levels,col=col,
                        hide=hide,omit=omit,omit.col=omit.col)
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

colourbar <- function(z=c(0,1),fill=c("#00FF0080","#FF000080"),
                      stroke='black',strip.width=0.02,clabel="",...){
    if (!validLevels(z)) return()
    if (all(is.na(fill)) | length(unique(fill))==1) col <- stroke
    else col <- fill
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

plot_points <- function(x,y,bg='yellow',pch=21,cex=1.5,pos,col,
                        show.numbers=FALSE,hide=NULL,omit=NULL,...){
    ns <- length(x)
    sn <- clear(1:ns,hide)
    X <- x[sn]
    Y <- y[sn]
    hascol <- !all(is.na(bg))
    show.points <- (hascol | !show.numbers)
    if (show.points){
        if (missing(pos)) pos <- 1
        if (missing(col)) col <- 'black'
        graphics::points(X,Y,bg=bg,pch=pch,col=col,cex=cex,...)
    } else {
        if (missing(pos)) pos <- NULL
        if (missing(col)){
            tcol <- rep('black',ns)
            tcol[omit] <- 'grey'
            col <- tcol[sn]
        }
    }
    if (show.numbers){
        graphics::text(X,Y,col=col,cex=cex,pos=pos,labels=sn,...)
    }
}

mymtext <- function(text,line=0,...){
    graphics::mtext(text,line=line,cex=graphics::par('cex'),...)
}

# if doall==FALSE, only returns the lower right submatrix
blockinverse <- function(AA,BB,CC,DD,invAA=NULL,doall=FALSE){
    if (is.null(invAA)) invAA <- solve(AA)
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
blockinverse3x3 <- function(AA,BB,CC,DD,EE,FF,GG,HH,II){
    ABDE <- rbind(cbind(AA,BB),cbind(DD,EE))
    invABDE <- blockinverse(AA=AA,BB=BB,CC=DD,DD=EE,doall=TRUE)
    CF <- rbind(CC,FF)
    GH <- cbind(GG,HH)
    blockinverse(AA=ABDE,BB=CF,CC=GH,DD=II,invAA=invABDE,doall=TRUE)
}

'%ni%' <- function(x,y)!('%in%'(x,y))

clear <- function(x,...,OGLS=FALSE){
    i <- unlist(list(...))
    if (is.matrix(x)) ns <- ifelse(OGLS,nrow(x)/2,nrow(x))
    else ns <- length(x)
    keep <- (1:ns)%ni%i
    if (length(i)>0){
        if (is.matrix(x) && OGLS) out <- subset_ogls(x,subset=keep)
        else out <- subset(x,subset=keep)
    } else {
        out <- x
    }
    out
}

logit <- function(x,m=0,M=1,inverse=FALSE){
    out <- x
    if (inverse){
        toobig <- x>1000
        toosmall <- x<(-1000)
        easy <- !(toobig|toosmall)
        out[easy] <- m + M*exp(x[easy])/(1+exp(x[easy]))
        out[toobig] <- m + M
        out[toosmall] <- m
    } else {
        easy <- (x>m & x<M)
        out[easy] <- log((x[easy]-m)/(M-x[easy]))
        out[x>=M] <- Inf
        out[x<=m] <- -Inf
    }
    out
}

invertible <- function(hess){
    tryCatch({
        E <- solve(hess)
        return(all(diag(E)>0))
    }, error = function(e){
        return(FALSE)
    })
}

inverthess <- function(hess){
    if (invertible(hess)){
        return(solve(hess))
    } else {
        H <- nearPD(hess)
        return(solve(H))
    }
}

det3x3 <- function(vx,vy,vz,sxy,sxz,syz){
    vx*vy*vz + 2*sxy*syz*sxz - vy*sxz^2 - vz*sxy^2 - vx*syz^2
}

invertcovmat <- function(vx,vy,vz,sxy=0,sxz=0,syz=0){
    if (missing(vz)){
        den <- vx*vy - sxy^2
        out <- cbind('xx'=vy/den,'yy'=vx/den,'xy'=-sxy/den)
    } else {
        aa <- vx
        ee <- vy
        ii <- vz
        bb <- dd <- sxy
        cc <- gg <- sxz
        ff <- hh <- syz
        den <- det3x3(vx,vy,vz,sxy,sxz,syz)
        xx <- (ee*ii-ff*hh)/den
        yy <- (aa*ii-cc*gg)/den
        zz <- (aa*ee-bb*dd)/den
        xy <- (cc*hh-bb*ii)/den
        xz <- (bb*ff-cc*ee)/den
        yz <- (cc*dd-aa*ff)/den
        out <- cbind('xx'=xx,'yy'=yy,'zz'=zz,'xy'=xy,'xz'=xz,'yz'=yz)
    }
    out
}

# recursive function that returns log(exp(u)+exp(v))
log_sum_exp <- function(u,v){
    if (missing(v)){
        if (length(u)>1){
            v <- u[-1]
            u <- u[1]
        } else {
            return(u)
        }
    }
    if (length(v)>1){
        v <- log_sum_exp(v[1],v[-1])
    } 
    max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))
}

contingencyfit <- function(par,fn,lower,upper,hessian=TRUE,control=NULL,...){
    fit <- stats::optim(par=par,fn=fn,method='L-BFGS-B',lower=lower,
                        upper=upper,hessian=hessian,control=control,...)
    failed <- fit$convergence>0 | (hessian & !invertible(fit$hessian)) |
        any((fit$par-lower)==0) | any((upper-fit$par)==0)
    if (failed){
        NMfit <- stats::optim(par=par,fn=fn,hessian=hessian,control=control,...)
        if (NMfit$convergence>0){
            warning('Optimisation did not converge.')
            if (NMfit$value<fit$value) {
                fit <- NMfit
            }
        } else {
            fit <- NMfit
        }
    }
    if (hessian & !invertible(fit$hessian)){
        warning('Ill-conditioned Hessian matrix')
    }
    fit
}

getMSWD <- function(X2,df){
    out <- list()
    if (df>0){
        out$mswd <- as.numeric(X2/df)
        out$p.value <- as.numeric(1-stats::pchisq(X2,df))
    } else {
        out$mswd <- 1
        out$p.value <- 1
    }
    out
}
