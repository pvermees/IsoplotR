#' Visualise U-Th-He data on a logratio plot or ternary diagram
#'
#' Plot U-Th(-Sm)-He data on a (log[He/Th] vs. log[U/He]) logratio
#' plot or U-Th-He ternary diagram
#'
#' @param x an object of class \code{UThHe}
#' @param logratio Boolean flag indicating whether the data should be
#'     shown on bivariate log[He/Th] vs. log[U/He] diagramme, or a
#'     U-Th-He ternary diagramme.
#' @param show.central.comp show the geometric mean composition as a
#'     white ellipse?
#' @param show.numbers show the grain numbers inside the error
#'     ellipses?
#' @param alpha confidence cutoff for the error ellipses
#' @param contour.col two-element vector with the fill colours to be
#'     assigned to the minimum and maximum age contour
#' @param ellipse.col background colour of the error ellipses
#' @param sigdig number of significant digits for the central age
#' @param xlim optional limits of the x-axis (log[U/He]) of the
#'     logratio plot. If \code{xlim=NA}, the axis limits are
#'     determined automatically.
#' @param ylim optional limits of the y-axis (log[Th/He]) of the
#'     logratio plot. If \code{ylim=NA}, the axis limits are
#'     determined automatically.
#' @param fact three-element vector with the scaling factors of the
#'     ternary diagram if \code{fact=NA}, these will be determined
#'     automatically
#' @param ... optional arguments to the generic \code{plot} function
#' @references
#' Vermeesch, P., 2010. HelioPlot, and the treatment of overdispersed
#' (U-Th-Sm)/He data. Chemical Geology, 271(3), pp.108-111.
#' @examples
#' data(examples)
#' helioplot(examples$UThHe)
#' dev.new()
#' helioplot(examples$UThHe,logratio=FALSE)
#' @export
helioplot <- function(x,logratio=TRUE,show.central.comp=TRUE,
                      show.numbers=FALSE,alpha=0.05,
                      contour.col=c('white','red'),
                      ellipse.col=rgb(0,1,0,0.5),
                      sigdig=2,xlim=NA,
                      ylim=NA,fact=NA,...){
    fit <- central.UThHe(x)
    if (logratio) {
        plot.logratio.contours(x,contour.col=contour.col,...)
        plot.logratio.ellipses(x,alpha=alpha,ellipse.col=ellipse.col,
                               show.numbers=show.numbers)
    } else {
        if (all(is.na(fact))) fact <- getfact(x,fit)
        plot.helioplot.contours(x,fact=fact,contour.col=contour.col,
                                xlim=xlim,ylim=ylim)
        plot.helioplot.ellipses(x,fact=fact,alpha=alpha,
                                ellipse.col=ellipse.col,
                                show.numbers=show.numbers)
    }
    if (show.central.comp){
        plot.central.ellipse(fit,fact=fact,logratio=logratio,
                             alpha=alpha,doSm=doSm(x))
        title(helioplot.title(fit,sigdig=sigdig))
    }
}

plot.logratio.frame <- function(lims,...){
    graphics::plot(lims[1:2],lims[3:4],type='n',bty='n',
                   xlab='log[U/He]',ylab='log[Th/He]',...)
}

plot.helioplot.frame <- function(lims,fact=c(1,1,1),fill.col=NA,...){
    graphics::plot(c(0,1),c(0,1),type='n',xaxt='n',yaxt='n',
                   xlab='',ylab='',asp=1,bty='n',...)
    corners <- xyz2xy(matrix(c(1,0,0,1,0,1,0,0,0,0,1,0),ncol=3))
    graphics::polygon(corners,col=fill.col)
    HeLabel <- paste0(fact[1],' x He')
    ULabel <- paste0(fact[3],' x U')
    ThLabel <- paste0(fact[2],' x Th')
    labels <- c(HeLabel,ULabel,ThLabel)
    graphics::text(corners[1:3,],labels=labels,pos=c(3,1,1))
}

plot.logratio.ellipses <- function(x,alpha=0.05,show.numbers=FALSE,
                                   ellipse.col=rgb(0,1,0,0.5)){
    valid <- is.finite(rowSums(x))
    X <- x[valid,]
    ns <- nrow(X)
    doSm <- doSm(X)
    for (i in 1:ns){
        if (doSm){
            uvwc <- UThHe2uvw.covmat(X,i)
            x0 <- uvwc$uvw[1]
            y0 <- uvwc$uvw[2]
            ell <- ellipse(x=x0,y=y0,covmat=uvwc$covmat,alpha=alpha)
        } else {
            uvc <- UThHe2uv.covmat(X,i)
            x0 <- uvc$uv[1]
            y0 <- uvc$uv[2]
            ell <- ellipse(x=x0,y=y0,covmat=uvc$covmat,alpha=alpha)
        }
        graphics::polygon(ell,col=ellipse.col)
        graphics::points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) graphics::text(x0,y0,i)
    }
   
}
plot.helioplot.ellipses <- function(x,fact=c(1,1,1),alpha=0.05,
                                    show.numbers=FALSE,
                                    ellipse.col=rgb(0,1,0,0.5)){
    valid <- is.finite(rowSums(x))
    X <- x[valid,]
    ns <- nrow(X)
    doSm <- doSm(X)
    for (i in 1:ns){
        if (doSm){
            uvwc <- UThHe2uvw.covmat(X,i)
            x0 <- uvwc$uvw[1]
            y0 <- uvwc$uvw[2]
            ell <- ellipse(x=x0,y=y0,covmat=uvwc$covmat,alpha=alpha)
            HeUTh0 <- uv2HeUTh(uvwc$uvw[1:2])
        } else {
            uvc <- UThHe2uv.covmat(X,i)
            x0 <- uvc$uv[1]
            y0 <- uvc$uv[2]
            ell <- ellipse(x=x0,y=y0,covmat=uvc$covmat,alpha=alpha)
            HeUTh0 <- uv2HeUTh(uvc$uv)
        }
        HeUTh <- uv2HeUTh(ell)
        xyz <- renormalise(HeUTh,fact=fact)
        xy <- xyz2xy(xyz)
        graphics::polygon(xy,col=ellipse.col)
        xyz0 <- renormalise(HeUTh0,fact=fact)
        x0y0 <- xyz2xy(xyz0)
        graphics::points(x0y0[1],x0y0[2],pch=19,cex=0.25)
        if (show.numbers) graphics::text(x0y0[1],x0y0[2],i)
    }
}

plot.central.ellipse <- function(fit,fact=c(1,1,1),logratio=TRUE,
                                 alpha=0.05,doSm=TRUE){
    if (doSm){
        ell <- ellipse(x=fit$uvw[1],y=fit$uvw[2],
                       covmat=fit$covmat[1:2,1:2],alpha=alpha)
    } else {
        ell <- ellipse(x=fit$uv[1],y=fit$uv[2],
                       covmat=fit$covmat,alpha=alpha)
    }
    if (logratio){
        graphics::polygon(ell,col='white')
    } else {
        HeUTh <- uv2HeUTh(ell)
        xyz <- renormalise(HeUTh,fact=fact)
        xy <- xyz2xy(xyz)
        graphics::polygon(xy,col='white')            
    }
}

plot.logratio.contours <- function(x,contour.col=c('white','red'),
                                   xlim=NA,ylim=NA,...){
    cntrs <- get.logratio.contours(x,xlim=xlim,ylim=ylim)
    plot.logratio.frame(cntrs$lims,...)
    tticks <- cntrs$tticks
    nt <- length(tticks)
    crp <- grDevices::colorRampPalette(contour.col)(nt)
    for (i in 1:nt){
        uv.plot <- cntrs$uv[[i]]
        uv.first <- uv.plot[1,]
        uv.last <- uv.plot[nrow(uv.plot),]
        uv.plot <- rbind(uv.plot,c(uv.last[1],cntrs$lims[3]))
        uv.plot <- rbind(uv.plot,cntrs$lims[c(1,3)])
        uv.plot <- rbind(uv.plot,c(cntrs$lims[1],uv.first[2]))
        uv.plot <- rbind(uv.plot,uv.first)
        graphics::polygon(uv.plot,col=crp[i])
        graphics::text(uv.first[1],uv.first[2],
                       labels=paste(tticks[i],'Ma'),pos=4)
    }
}

plot.helioplot.contours <- function(x,fact=c(1,1,1),
                                    contour.col=c('white','red'),
                                    xlim=NA,ylim=NA,...){
    cntrs <- get.helioplot.contours(x,fact=fact)
    plot.helioplot.frame(cntrs$lims,fact=fact,
                         fill.col=contour.col[2],...)
    tticks <- cntrs$tticks
    nt <- length(tticks)
    crp <- grDevices::colorRampPalette(contour.col)(nt+1)
    for (i in nt:1){
        xyz <- cntrs$xyz[[i]]
        xyz.first <- xyz[1,]
        xyz <- rbind(xyz,c(0,0,1),c(0,1,0),xyz.first)
        xyz <- renormalise(xyz,fact=fact)
        xy <- xyz2xy(xyz)
        graphics::polygon(xy,col=crp[i])
        label <- paste(signif(cntrs$tticks[[i]],2),'Ma')
        graphics::text(xy[1,1],xy[1,2],labels=label,pos=2)
    }
}

helioplot.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age[1],fit$age[2],sigdig=sigdig)
    line1 <- substitute('central age ='~a%+-%b~'[Ma] (1'~sigma~')',
                        list(a=rounded.age[1], b=rounded.age[2]))
    line2 <- substitute('MSWD (concordance) ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,2),
                             b=signif(fit$p.value,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}

# UVW = central composition
SS.UThHe.uvw <- function(UVW,x){
    ns <- nrow(x)
    SS <- 0
    for (i in 1:ns){
        uvwc <- UThHe2uvw.covmat(x,i)
        X <- UVW-uvwc$uvw
        Ei <- solve(uvwc$covmat)
        SSi <- X %*% Ei %*% t(X)
        if (is.finite(SSi)) SS <- SS + SSi
    }
    SS
}
SS.UThHe.uv <- function(UV,x){
    ns <- nrow(x)
    SS <- 0
    for (i in 1:ns){
        uvc <- UThHe2uv.covmat(x,i)
        X <- UV-uvc$uv
        Ei <- solve(uvc$covmat)
        SSi <- X %*% Ei %*% t(X)
        if (is.finite(SSi)) SS <- SS + SSi
    }
    SS
}

get.logratio.contours <- function(x,xlim=NA,ylim=NA,res=500){
    out <- list()
    R <- iratio('U238U235')[1]
    L8 <- lambda('U238')[1]
    L5 <- lambda('U235')[1]
    L2 <- lambda('Th232')[1]
    L7 <- lambda('Sm147')[1]
    f147 <- f147Sm()[1]
    doSm <- doSm(x)
    if (doSm){
        uvw <- UThHe2uvw(x)
        out$lims <- get.logratioplot.limits(uvw[,c('u','v')])
        w <- mean(uvw[,'w'],na.rm=TRUE)
        Sm <- geomean.Sm(x)
    } else {
        uv <- UThHe2uv(x)
        out$lims <- get.logratioplot.limits(uv)
        w <- 0
        Sm <- 0
    }
    if (all(is.finite(xlim))) out$lims[1:2] <- xlim
    if (all(is.finite(ylim))) out$lims[3:4] <- ylim
    du <- out$lims[2]-out$lims[1]
    dv <- out$lims[4]-out$lims[3]
    tticks <- get.logratio.tticks(out$lims,Sm=Sm)
    out$uv <- list()
    out$tticks <- NULL
    nt <- 0
    for (i in 1:length(tticks)){
        tt <- tticks[i]
        aa <- 8*(exp(L8*tt)-1)*R/(1+R) +
             7*(exp(L5*tt)-1)/(1+R)
        bb <- 6*(exp(L2*tt)-1)
        if (doSm)
            cc <- f147*(exp(L7*tt)-1)
        else
            cc <- 0
        uv.plot <- NULL
        for (j in 1:res){
            u <- out$lims[1]+du*j/res
            exp.v <- (1-aa*exp(u)-cc*exp(w))/bb
            if (exp.v > exp(out$lims[3]) & exp.v < exp(out$lims[4])){
                v <- log(exp.v)
                uv.plot <- rbind(uv.plot,c(u,v))
            }
        }
        if (!is.null(uv.plot)){
            nt <- nt + 1
            out$tticks[nt] <- tt
            out$uv[[nt]] <- uv.plot
        }
    }
    out
}
get.helioplot.contours <- function(x,fact=c(1,1,1),res=100){
    out <- list()
    doSm <- doSm(x)
    if (doSm){
        uvw <- UThHe2uvw(x)
        SmU <- exp(mean(uvw[,'w']-uvw[,'u'],na.rm=TRUE))
    }
    out$tticks <- get.helioplot.tticks(fact)
    out$xyz <- list()
    nt <- length(out$tticks)
    for (i in 1:nt){
        tt <- out$tticks[i]
        U <- seq(1,0,length.out=res)
        Th <- seq(0,1,length.out=res)
        if (doSm) Sm <- SmU*U
        else Sm <- 0
        He <- get.He(tt,U,Th,Sm)
        out$xyz[[i]] <- cbind(He,U,Th)
    }
    out
}

# f = the distance from the X- and Y-margins for
# the maximum and minimum age contours to be plotted
get.logratio.tticks <- function(lims,f=0.05,Sm=0){
    uv.lims <- rep(0,4)
    uv.lims[1] <- lims[1]+f*(lims[3]-lims[1])
    uv.lims[2] <- lims[2]+f*(lims[4]-lims[2])
    uv.lims[3] <- lims[3]-f*(lims[3]-lims[1])
    uv.lims[4] <- lims[4]-f*(lims[4]-lims[2])
    Umax <- exp(uv.lims[1])/(exp(uv.lims[1])+exp(uv.lims[3])+1)
    Thmax <- exp(uv.lims[3])/(exp(uv.lims[1])+exp(uv.lims[3])+1)
    Hemax <- 1/(exp(uv.lims[1])+exp(uv.lims[3])+1)
    Umin <- exp(uv.lims[2])/(exp(uv.lims[2])+exp(uv.lims[4])+1)
    Thmin <- exp(uv.lims[4])/(exp(uv.lims[2])+exp(uv.lims[4])+1)
    Hemin <- 1/(exp(uv.lims[2])+exp(uv.lims[4])+1)
    mint <- get.UThHe.age(Umin,0,Thmin,0,Hemin,0,Sm,0)[1]
    maxt <- get.UThHe.age(Umax,0,Thmax,0,Hemax,0,Sm,0)[1]
    get.tticks(mint,maxt)
}
# f = the distance from the ternary corners for
# the maximum and minimum age contours to be plotted
get.helioplot.tticks <- function(fact,f=0.05,Sm=0){
    m <- c(f,1-f,0)/fact
    M <- c(1-f,f,0)/fact
    mint <- get.UThHe.age(m[2],0,m[3],0,m[1],0,Sm,0)[1]
    maxt <- get.UThHe.age(M[2],0,M[3],0,M[1],0,Sm,0)[1]
    get.tticks(mint,maxt)
}
get.tticks <- function(mint,maxt){
    m <- log10(mint)
    M <- log10(maxt)
    grDevices::axisTicks(usr=c(m,M),log=TRUE)
}

get.logratioplot.limits <- function(uv,f=1){
    ru <- range(uv[,1],na.rm=TRUE)
    rv <- range(uv[,2],na.rm=TRUE)
    du <- diff(ru)
    dv <- diff(rv)
    minu <- ru[1] - f*du
    maxu <- ru[2] + f*du
    minv <- rv[1] - f*dv
    maxv <- rv[2] + f*dv
    c(minu,maxu,minv,maxv)
}

# x is an object of class UThHe
UThHe2uvw <- function(x){
    if (hasClass(x,'UThHe'))
        out <- log(x[,c('U','Th','Sm')])-log(x[,'He'])
    else
        out <- matrix(log(x[c('U','Th','Sm')])-log(x['He']),1,3)
    colnames(out) <- c('u','v','w')
    out
}
UThHe2uv <- function(x){
    if (hasClass(x,'UThHe'))
        out <- log(x[,c('U','Th')])-log(x[,'He'])
    else
        out <- matrix(log(x[c('U','Th')])-log(x['He']),1,2)
    colnames(out) <- c('u','v')
    out
}

uvw2UThHe <- function(uvw,covmat=matrix(0,3,3)){
    u <- uvw[1]
    v <- uvw[2]
    w <- uvw[3]
    D <- exp(u)+exp(v)+exp(w)+1
    U <- exp(u)/D
    Th <- exp(v)/D
    Sm <- exp(w)/D
    He <- 1/D
    J <- matrix(0,4,3)
    J[1,1] <- exp(u)*(exp(v)+exp(w)+1)/D^2
    J[1,2] <- -exp(u)*exp(v)/D^2
    J[1,3] <- -exp(u)*exp(w)/D^2
    J[2,1] <- -exp(v)*exp(u)/D^2
    J[2,2] <- exp(v)*(exp(u)+exp(w)+1)/D^2
    J[2,3] <- -exp(v)*exp(w)/D^2
    J[3,1] <- -exp(w)*exp(u)/D^2
    J[3,2] <- -exp(w)*exp(v)/D^2
    J[3,3] <- exp(w)*(exp(u)+exp(v)+1)/D^2
    J[4,1] <- -exp(u)/D^2
    J[4,2] <- -exp(v)/D^2
    J[4,3] <- -exp(w)/D^2
    E <- J %*% covmat %*% t(J)
    sU <- sqrt(E[1,1])
    sTh <- sqrt(E[2,2])
    sSm <- sqrt(E[3,3])
    sHe <- sqrt(E[4,4])
    out <- c(U,sU,Th,sTh,Sm,sSm,He,sHe)
    names(out) <- c('U','sU','Th','sTh','Sm','sSm','He','sHe')
    out
}
uv2UThHe <- function(uv,covmat=matrix(0,2,2)){
    u <- uv[1]
    v <- uv[2]
    D <- exp(u)+exp(v)+1
    U <- exp(u)/D
    Th <- exp(v)/D
    He <- 1/D
    J <- matrix(0,3,2)
    J[1,1] <- exp(u)*(exp(v)+1)/D^2
    J[1,2] <- -exp(u)*exp(v)/D^2
    J[2,1] <- -exp(v)*exp(u)/D^2
    J[2,2] <- exp(v)*(exp(u)+1)/D^2
    J[3,1] <- -exp(u)/D^2
    J[3,2] <- -exp(v)/D^2
    E <- J %*% covmat %*% t(J)
    sU <- sqrt(E[1,1])
    sTh <- sqrt(E[2,2])
    sHe <- sqrt(E[3,3])
    out <- c(U,sU,Th,sTh,He,sHe)
    names(out) <- c('U','sU','Th','sTh','He','sHe')
    out
}
uv2HeUTh <- function(uv){
    if (hasClass(uv,"matrix")){
        u <- uv[,1]
        v <- uv[,2]
    } else {
        u <- uv[1]
        v <- uv[2]
    }
    D <- exp(u)+exp(v)+1
    U <- exp(u)/D
    Th <- exp(v)/D
    He <- 1/D
    cbind(He,U,Th)
}

UThHe2uvw.covmat <- function(x,i){
    out <- list()
    U <- x[i,'U']
    sU <- x[i,'errU']
    Th <- x[i,'Th']
    sTh <- x[i,'errTh']
    Sm <- x[i,'Sm']
    sSm <- x[i,'errSm']
    He <- x[i,'He']
    sHe <- x[i,'errHe']
    out$uvw <- UThHe2uvw(x[i,])
    out$covmat <- matrix(0,3,3)
    names(out$uvw) <- c("u","v","w")
    rownames(out$covmat) <- c("u","v","w")
    colnames(out$covmat) <- c("u","v","w")
    out$covmat[1,1] <- (sU/U)^2 + (sHe/(U*He))^2
    out$covmat[2,2] <- (sTh/Th)^2 + (sHe/(Th*He))^2
    out$covmat[3,3] <- (sSm/Sm)^2 + (sHe/(Sm*He))^2
    out$covmat[1,2] <- (sHe^2)/(U*Th*He^2)
    out$covmat[1,3] <- (sHe^2)/(U*Sm*He^2)
    out$covmat[2,3] <- (sHe^2)/(Th*Sm*He^2)
    out$covmat[2,1] <- out$covmat[1,2]
    out$covmat[3,1] <- out$covmat[1,3]
    out$covmat[3,2] <- out$covmat[2,3]
    out
}
UThHe2uv.covmat <- function(x,i){
    out <- list()
    U <- x[i,'U']
    sU <- x[i,'errU']
    Th <- x[i,'Th']
    sTh <- x[i,'errTh']
    He <- x[i,'He']
    sHe <- x[i,'errHe']
    out$uv <- UThHe2uv(x[i,])
    out$covmat <- matrix(0,2,2)
    names(out$uv) <- c("u","v")
    rownames(out$covmat) <- c("u","v")
    colnames(out$covmat) <- c("u","v")
    out$covmat[1,1] <- (sU/U)^2 + (sHe/(U*He))^2
    out$covmat[2,2] <- (sTh/Th)^2 + (sHe/(Th*He))^2
    out$covmat[1,2] <- (sHe^2)/(U*Th*He^2)
    out$covmat[2,1] <- out$covmat[1,2]
    out
}

# ternary compositions to plot coordinates
xyz2xy <- function(xyz){
    if (hasClass(xyz,"matrix")){
        n <- nrow(xyz)
        x <- xyz[,1]
        y <- xyz[,2]
        z <- xyz[,3]
    } else {
        n <- 1
        x <- xyz[1]
        y <- xyz[2]
        z <- xyz[3]
    }
    xy <- matrix(0,n,2)
    xy[,1] <- 0.5*(x+2*z)/(x+y+z)
    xy[,2] <- sin(pi/3)*x/(x+y+z)
    xy
}

renormalise <- function(xyz,fact=c(1,1,1)){
    if (hasClass(xyz,"matrix")){
        nr <- nrow(xyz)
        FACT <- matrix(rep(fact,nr),nrow=nr,byrow=TRUE)
        out <- xyz*FACT
        NORM <- rowSums(out)
        out <- out/NORM
    } else {
        out <- xyz*fact/sum(xyz*fact)
    }
    out
}

geomean.Sm <- function(x){
    UVW <- colMeans(UThHe2uvw(x),na.rm=TRUE)
    exp(UVW[3])/(exp(UVW[1])+exp(UVW[2])+exp(UVW[3])+1)
}

getfact <- function(x,fit){
    if (doSm(x)) HeUTh <- uv2HeUTh(fit$uvw[1:2])
    else HeUTh <- uv2HeUTh(fit$uv[1:2])
    fact <- signif(1/HeUTh,1)
}
