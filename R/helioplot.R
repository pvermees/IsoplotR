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
#' @param levels a vector with additional values to be displayed as
#'     different background colours within the error ellipses.
#' @param ellipse.col a vector of two background colours for the error
#'     ellipses. If \code{levels=NA}, then only the first colour will
#'     be used. If \code{levels} is a vector of numbers, then
#'     \code{ellipse.col} is used to construct a colour ramp.
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
                      contour.col=c('white','red'),levels=NA,
                      ellipse.col=c("#00FF0080","#0000FF80"),
                      sigdig=2,xlim=NA,ylim=NA,fact=NA,...){
    fit <- central.UThHe(x)
    if (logratio) {
        plot_logratio_contours(x,contour.col=contour.col,
                               xlim=xlim,ylim=ylim,...)
        plot_logratio_ellipses(x,alpha=alpha,levels=levels,
                               ellipse.col=ellipse.col,
                               show.numbers=show.numbers)
    } else {
        if (all(is.na(fact))) fact <- getfact(x,fit)
        plot_helioplot_contours(x,fact=fact,contour.col=contour.col,
                                xlim=xlim,ylim=ylim)
        plot_helioplot_ellipses(x,fact=fact,alpha=alpha,
                                levels=levels,ellipse.col=ellipse.col,
                                show.numbers=show.numbers)
    }
    if (show.central.comp){
        plot_central_ellipse(fit,fact=fact,logratio=logratio,
                             alpha=alpha,doSm=doSm(x))
        graphics::title(helioplot_title(fit,sigdig=sigdig))
    }
    colourbar(z=levels,col=ellipse.col)
}

plot_logratio_frame <- function(lims,...){
    graphics::plot(lims[1:2],lims[3:4],type='n',bty='n',
                   xlab='log[U/He]',ylab='log[Th/He]',...)
}

plot_helioplot_frame <- function(lims,fact=c(1,1,1),fill.col=NA,...){
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

plot_logratio_ellipses <- function(x,alpha=0.05,show.numbers=FALSE,levels=NA,
                                   ellipse.col=c("#00FF0080","#0000FF80")){
    ns <- nrow(x)
    ellipse.cols <- set.ellipse.colours(ns=ns,levels=levels,col=ellipse.col)
    for (i in 1:ns){
        uvc <- UThHe2uv.covmat(x,i)
        x0 <- uvc$uv[1]
        y0 <- uvc$uv[2]
        ell <- ellipse(x=x0,y=y0,covmat=uvc$covmat,alpha=alpha)
        graphics::polygon(ell,col=ellipse.cols[i])
        if (show.numbers) graphics::text(x0,y0,i)
        else graphics::points(x0,y0,pch=19,cex=0.25)
    }
}
plot_helioplot_ellipses <- function(x,fact=c(1,1,1),alpha=0.05,
                                    show.numbers=FALSE,levels=NA,
                                    ellipse.col=c("#00FF0080","#0000FF80")){
    ns <- nrow(x)
    ellipse.cols <- set.ellipse.colours(ns=ns,levels=levels,col=ellipse.col)
    for (i in 1:ns){
        uvc <- UThHe2uv.covmat(x,i)
        x0 <- uvc$uv[1]
        y0 <- uvc$uv[2]
        ell <- ellipse(x=x0,y=y0,covmat=uvc$covmat,alpha=alpha)
        HeUTh0 <- uv2HeUTh(uvc$uv)
        HeUTh <- uv2HeUTh(ell)
        xyz <- renormalise(HeUTh,fact=fact)
        xy <- xyz2xy(xyz)
        graphics::polygon(xy,col=ellipse.cols[i])
        xyz0 <- renormalise(HeUTh0,fact=fact)
        x0y0 <- xyz2xy(xyz0)
        if (show.numbers) graphics::text(x0y0[1],x0y0[2],i)
        else graphics::points(x0y0[1],x0y0[2],pch=19,cex=0.25)
    }
}

plot_central_ellipse <- function(fit,fact=c(1,1,1),logratio=TRUE,
                                 alpha=0.05,doSm=TRUE,...){
    ell <- ellipse(x=fit$uvw[1],y=fit$uvw[2],
                   covmat=fit$covmat[1:2,1:2],alpha=alpha)
    if (logratio){
        graphics::polygon(ell,col='white')
    } else {
        HeUTh <- uv2HeUTh(ell)
        xyz <- renormalise(HeUTh,fact=fact)
        xy <- xyz2xy(xyz)
        graphics::polygon(xy,col='white')            
    }
}

plot_logratio_contours <- function(x,contour.col=c('white','red'),
                                   xlim=NA,ylim=NA,...){
    cntrs <- get.logratio.contours(x,xlim=xlim,ylim=ylim)
    plot_logratio_frame(cntrs$lims,...)
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

plot_helioplot_contours <- function(x,fact=c(1,1,1),
                                    contour.col=c('white','red'),
                                    xlim=NA,ylim=NA,...){
    cntrs <- get.helioplot.contours(x,fact=fact)
    plot_helioplot_frame(cntrs$lims,fact=fact,
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

helioplot_title <- function(fit,sigdig=2){
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

get.logratio.contours <- function(x,xlim=NA,ylim=NA,res=50){
    out <- list()
    R <- iratio('U238U235')[1]
    L8 <- lambda('U238')[1]
    L5 <- lambda('U235')[1]
    L2 <- lambda('Th232')[1]
    L7 <- lambda('Sm147')[1]
    f147 <- f147Sm()[1]
    doSm <- doSm(x)
    out$lims <- get.logratioplot.limits(x)
    if (doSm){
        uvw <- UThHe2uvw(x)        
        w <- mean(uvw[,'w'],na.rm=TRUE)
        Sm <- geomean.Sm(x)
    } else {
        uv <- UThHe2uv(x)
        w <- 0
        Sm <- 0
        cc <- 0
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
        aa <- 8*(exp(L8*tt)-1)*R/(1+R) + 7*(exp(L5*tt)-1)/(1+R)
        bb <- 6*(exp(L2*tt)-1)
        if (doSm) cc <- f147*(exp(L7*tt)-1)
        uv.plot <- NULL
        # evaluate the maximum v value
        pred.exp.u <- 1-bb*exp(out$lims[4])-cc*exp(w)
        if (pred.exp.u > 0){
            u4maxv <- log(pred.exp.u) - log(aa)
            if (u4maxv > out$lims[1] && u4maxv < out$lims[2])
                uv.plot <- rbind(uv.plot,c(u4maxv,out$lims[4]))
        }
        # evaluate all the whole range of u values
        for (j in 0:res){
            u <- out$lims[1]+du*j/res
            exp.v <- (1-aa*exp(u)-cc*exp(w))/bb
            if (exp.v > exp(out$lims[3]) & exp.v < exp(out$lims[4])){
                v <- log(exp.v)
                uv.plot <- rbind(uv.plot,c(u,v))
            }
        }
        # evaluate the minimum v value
        pred.exp.u <- 1-bb*exp(out$lims[3])-cc*exp(w)
        if (pred.exp.u > 0){
            u4minv <- log(pred.exp.u) - log(aa)
            if (u4minv > out$lims[1] && u4minv < out$lims[2])
                uv.plot <- rbind(uv.plot,c(u4minv,out$lims[3]))
        }
        # add to the list if any solutions were found
        if (!is.null(uv.plot)){
            nt <- nt + 1
            out$tticks[nt] <- tt
            out$uv[[nt]] <- uv.plot
        }
    }
    out
}
get.helioplot.contours <- function(x,fact=c(1,1,1),res=50){
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

get.logratioplot.limits <- function(x,nse=3){
    ns <- length(x)
    doSm <- doSm(x)
    minu <- Inf
    maxu <- -Inf
    minv <- Inf
    maxv <- -Inf
    for (i in 1:ns){
        d <- UThHe2uv.covmat(x,i)
        uv <- d$uv
        uv.err <- sqrt(diag(d$covmat)[c('u','v')])
        umin <- uv['u'] - nse*uv.err['u']
        umax <- uv['u'] + nse*uv.err['u']
        vmin <- uv['v'] - nse*uv.err['v']
        vmax <- uv['v'] + nse*uv.err['v']
        if (umax>maxu) maxu <- umax
        if (umin<minu) minu <- umin
        if (vmax>maxv) maxv <- vmax
        if (vmin<minv) minv <- vmin
    }
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
    U <- x[i,'U']
    sU <- x[i,'errU']
    Th <- x[i,'Th']
    sTh <- x[i,'errTh']
    Sm <- x[i,'Sm']
    sSm <- x[i,'errSm']
    He <- x[i,'He']
    sHe <- x[i,'errHe']
    out <- list()
    out$uvw <- UThHe2uvw(x[i,])
    out$covmat <- matrix(0,3,3)
    J <- matrix(0,3,4)
    E <- matrix(0,4,4)
    diag(E) <- c(sU,sTh,sSm,sHe)^2
    J[1,1] <- 1/U   # du.dU
    J[1,4] <- -1/He # du.dHe
    J[2,2] <- 1/Th  # dv.dTh
    J[2,4] <- -1/He # dv.dHe
    J[3,3] <- 1/Sm  # dw.dSm
    J[3,4] <- -1/He # dw.dHe
    out$covmat <- J %*% E %*% t(J)
    names(out$uvw) <- c("u","v","w")
    rownames(out$covmat) <- c("u","v","w")
    colnames(out$covmat) <- c("u","v","w")
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
    J <- matrix(0,2,3)
    E <- matrix(0,3,3)
    diag(E) <- c(sU,sTh,sHe)^2
    J[1,1] <- 1/U   # du.dU
    J[1,3] <- -1/He # du.dHe
    J[2,2] <- 1/Th  # dv.dTh
    J[2,3] <- -1/He # dv.dHe
    out$covmat <- J %*% E %*% t(J)
    names(out$uv) <- c("u","v")
    rownames(out$covmat) <- c("u","v")
    colnames(out$covmat) <- c("u","v")
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
    HeUTh <- uv2HeUTh(fit$uvw[1:2])
    fact <- signif(1/HeUTh,1)
}
