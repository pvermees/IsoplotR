#' Visualise U-Th-He data on a logratio plot or ternary diagram
#'
#' Plot U-Th(-Sm)-He data on a (log[He/Th] vs. log[U/He]) logratio
#' plot or U-Th-He ternary diagram
#'
#' @param x an object of class \code{UThHe}
#' @param logratio Boolean flag indicating whether the data should be
#'     shown on bivariate log[He/Th] vs. log[U/He] diagramme, or a
#'     U-Th-He ternary diagramme.
#' @param show.age calculate the U-Th-He central age?
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
#' @param fact three element vector with the scaling factors of the
#'     ternary diagram if \code{fact=NA}, these will be determined
#'     automatically
#' @param ... optional arguments to the generic \code{plot} function
#' @examples
#' data(examples)
#' helioplot(examples$UThHe)
#' dev.new()
#' helioplot(examples$UThHe,logratio=FALSE)
#' @export
helioplot <- function(x,logratio=TRUE,show.age=TRUE,
                      show.numbers=FALSE,alpha=0.05,
                      contour.col=c('white','red'),
                      ellipse.col=rgb(0,1,0,0.5),sigdig=2, xlim=NA,
                      ylim=NA,fact=NA,...){
    fit <- centralage.UThHe(x)
    if (logratio) {
        plot.logratio.contours(x,contour.col=contour.col,...)
        plot.logratio.ellipses(x,alpha=alpha,ellipse.col=ellipse.col,
                               show.numbers=show.numbers)
    } else {
        if (all(is.na(fact))){
            HeUTh <- uv2HeUTh(fit$uvw[1:2])
            fact <- signif(1/HeUTh,1)
        }
        plot.helioplot.contours(x,fact=fact,contour.col=contour.col,
                                xlim=xlim,ylim=ylim)
        plot.helioplot.ellipses(x,fact=fact,alpha=alpha,
                                ellipse.col=ellipse.col,
                                show.numbers=show.numbers)
    }
    if (show.age){
        plot.central.ellipse(fit,fact=fact,logratio=logratio,
                             alpha=alpha)
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
    ns <- nrow(x)
    for (i in 1:ns){
        uvwc <- UThHe2uvw.covmat(x,i)
        x0 <- uvwc$uvw[1]
        y0 <- uvwc$uvw[2]
        ell <- ellipse(x=x0,y=y0,covmat=uvwc$covmat,alpha=alpha)
        graphics::polygon(ell,col=ellipse.col)
        graphics::points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) graphics::text(x0,y0,i)
    }
}

plot.helioplot.ellipses <- function(x,fact=c(1,1,1),alpha=0.05,
                                    show.numbers=FALSE,
                                    ellipse.col=rgb(0,1,0,0.5)){
    ns <- nrow(x)
    for (i in 1:ns){
        uvwc <- UThHe2uvw.covmat(x,i)
        x0 <- uvwc$uvw[1]
        y0 <- uvwc$uvw[2]
        ell <- ellipse(x=x0,y=y0,covmat=uvwc$covmat,alpha=alpha)
        HeUTh <- uv2HeUTh(ell)
        xyz <- renormalise(HeUTh,fact=fact)
        xy <- xyz2xy(xyz)
        graphics::polygon(xy,col=ellipse.col)
        HeUTh0 <- uv2HeUTh(uvwc$uvw[1:2])
        xyz0 <- renormalise(HeUTh0,fact=fact)
        x0y0 <- xyz2xy(xyz0)
        graphics::points(x0y0[1],x0y0[2],pch=19,cex=0.25)
        if (show.numbers) graphics::text(x0y0[1],x0y0[2],i)
    }
}

plot.central.ellipse <- function(fit,fact=c(1,1,1),logratio=TRUE,
                                 alpha=0.05){
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
                        list(a=rounded.age$x, b=rounded.age$err))
    line2 <- substitute('MSWD (concordance) ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,2),
                             b=signif(fit$p.value,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}

centralage.UThHe <- function(x){
    uvw <- UThHe2uvw(x)
    fit <- stats::optim(c(0,0,0),SS.UThHe,method='BFGS',
                        hessian=TRUE,x=x)
    out <- list()
    ns <- nrow(x)
    df <- 2*(ns-1)
    out$uvw <- fit$par
    out$covmat <- solve(fit$hessian)
    nms <- c('u','v','w')
    names(out$uvw) <- nms
    colnames(out$covmat) <- nms
    rownames(out$covmat) <- nms
    SS <- SS.UThHe(out$uvw,x,Sm=FALSE)
    out$mswd <- SS/df
    out$p.value <- 1-pchisq(SS,df)
    cc <- uvw2UThHe(out$uvw,out$covmat)
    out$age <- get.UThHe.age(cc['U'],cc['sU'],cc['Th'],cc['sTh'],
                             cc['He'],cc['sHe'],cc['Sm'],cc['sSm'])
    out
}

# UVW = central composition
SS.UThHe <- function(UVW,x,Sm=TRUE){
    ns <- nrow(x)
    SS <- 0
    for (i in 1:ns){
        uvwc <- UThHe2uvw.covmat(x,i)
        X <- UVW-uvwc$uvw
        Ei <- solve(uvwc$covmat)
        if (!Sm) {
            X <- matrix(X[1:2],1,2)
            Ei <- Ei[1:2,1:2]
        }
        SS <- SS + X %*% Ei %*% t(X)
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
    uvw <- UThHe2uvw(x)
    f <- 0.5
    out$lims <- get.logratioplot.limits(uvw[,c('u','v')],f=f)
    if (all(!is.na(xlim))) out$lims[1:2] <- xlim
    if (all(!is.na(ylim))) out$lims[3:4] <- ylim
    du <- out$lims[2]-out$lims[1]
    dv <- out$lims[4]-out$lims[3]
    out$uv <- list()
    tticks <- get.logratio.tticks(out$lims)
    nt <- length(tticks)
    out$tticks <- rep(0,nt)
    for (i in 1:nt){
        tt <- tticks[i]
        a <- 8*R*(exp(L8*tt)-1)/(1+R) +
             7*(exp(L5*tt)-1)/(1+R)
        b <- 6*(exp(L2*tt)-1)
        uv.plot <- NULL
        for (j in 1:res){
            u <- out$lims[1]+du*j/res
            exp.v <- (1-a*exp(u))/b
            if (exp.v > exp(out$lims[3]) & exp.v < exp(out$lims[4])){
                v <- log(exp.v)
                uv.plot <- rbind(uv.plot,c(u,v))
            }
        }
        out$tticks[i] <- tt
        out$uv[[i]] <- uv.plot
    }
    out
}

get.helioplot.contours <- function(x,fact=c(1,1,1),res=100){
    out <- list()
    out$tticks <- get.helioplot.tticks(fact)
    out$xyz <- list()
    nt <- length(out$tticks)
    for (i in 1:nt){
        tt <- out$tticks[i]
        U <- seq(1,0,length.out=res)
        Th <- seq(0,1,length.out=res)
        He <- get.He(tt,U,Th)
        out$xyz[[i]] <- cbind(He,U,Th)
    }
    out
}

# f = the distance from the ternary corners for
# the maximum and minimum age contours to be plotted
get.helioplot.tticks <- function(fact,f=0.05){
    m <- c(f,1-f,0)/fact
    M <- c(1-f,f,0)/fact
    mint <- get.UThHe.age(m[2],0,m[3],0,m[1],0)[1]
    maxt <- get.UThHe.age(M[2],0,M[3],0,M[1],0)[1]
    get.tticks(mint,maxt)
}
# f = the distance from the X- and Y-margins for
# the maximum and minimum age contours to be plotted
get.logratio.tticks <- function(lims,f=0.05){
    tlims <- rep(0,4)
    tlims[1] <- lims[1]+f*(lims[3]-lims[1])
    tlims[2] <- lims[2]+f*(lims[4]-lims[2])
    tlims[3] <- lims[3]-f*(lims[3]-lims[1])
    tlims[4] <- lims[4]-f*(lims[4]-lims[2])
    Umax <- exp(tlims[1])/(exp(tlims[1])+exp(tlims[2])+1)
    Thmax <- exp(tlims[2])/(exp(tlims[1])+exp(tlims[2])+1)
    Hemax <- 1/(exp(tlims[1])+exp(tlims[2])+1)
    Umin <- exp(tlims[3])/(exp(tlims[3])+exp(tlims[4])+1)
    Thmin <- exp(tlims[4])/(exp(tlims[3])+exp(tlims[4])+1)
    Hemin <- 1/(exp(tlims[3])+exp(tlims[4])+1)
    mint <- get.UThHe.age(Umin,0,Thmin,0,Hemin,0)[1]
    maxt <- get.UThHe.age(Umax,0,Thmax,0,Hemax,0)[1]
    get.tticks(mint,maxt)
}
get.tticks <- function(mint,maxt){
    m <- log10(mint)
    M <- log10(maxt)
    grDevices::axisTicks(usr=c(m,M),log=TRUE)
}

get.logratioplot.limits <- function(uv,f=1){
    ru <- range(uv[,1])
    rv <- range(uv[,2])
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
    if (methods::is(x,'UThHe'))
        out <- log(x[,c('U','Th','Sm')])-log(x[,'He'])
    else
        out <- matrix(log(x[c('U','Th','Sm')])-log(x['He']),1,3)
    colnames(out) <- c('u','v','w')
    out
}

uv2HeUTh <- function(uv){
    if (methods::is(uv,"matrix")){
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

# ternary compositions to plot coordinates
xyz2xy <- function(xyz){
    if (methods::is(xyz,"matrix")){
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
    if (methods::is(xyz,"matrix")){
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
