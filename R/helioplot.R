# x is an object of class UThHe
helioplot <- function(x,f=0.5,...){
    graphics::plot(c(0,1),c(0,1),type='n',xaxt='n',yaxt='n',
                   xlab='',ylab='',asp=1,bty='n',...)
    corners <- xyz2xy(matrix(c(1,0,0,1,0,1,0,0,0,0,1,0),ncol=3))
    graphics::lines(corners)
    graphics::text(corners[1:3,],labels=c('He','U','Th'),pos=c(3,1,1))
}

get.UThHe.age.limits <- function(lims){
    Umax <- exp(lims[1])/(exp(lims[1])+exp(lims[2])+1)
    Thmax <- exp(lims[2])/(exp(lims[1])+exp(lims[2])+1)
    Hemax <- 1/(exp(lims[1])+exp(lims[2])+1)
    Umin <- exp(lims[3])/(exp(lims[3])+exp(lims[4])+1)
    Thmin <- exp(lims[4])/(exp(lims[3])+exp(lims[4])+1)
    Hemin <- 1/(exp(lims[3])+exp(lims[4])+1)
    mint <- get.UThHe.age(Umin,0,Thmin,0,Hemin,0)[1]
    maxt <- get.UThHe.age(Umax,0,Thmax,0,Hemax,0)[1]
    c(mint,maxt)
}

logratioplot <- function(x,show.age=TRUE,alpha=0.05,
                         contour.col=c('white','red'),
                         ellipse.col=rgb(0,1,0,0.5),
                         xlim=NA,ylim=NA,sigdig=2,...){
    plot.logratio.contours(x,contour.col=contour.col,
                           xlim=xlim,ylim=ylim)
    plot.logratio.ellipses(x,alpha=alpha,
                           ellipse.col=ellipse.col)
    if (show.age){
        fit <- central.age.UThHe(x)
        ell <- ellipse(x=fit$uvw[1],y=fit$uvw[2],
                       covmat=fit$covmat[1:2,1:2],alpha=alpha)
        graphics::polygon(ell,col='white')
        title(helioplot.title(fit,sigdig=sigdig))
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

central.age.UThHe <- function(x){
    uvw <- UThHe2uvw(x)
    fit <- optim(c(0,0,0),SS.UThHe,method='BFGS',hessian=TRUE,x=x)
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

plot.logratio.ellipses <- function(x,alpha=0.05,
                                   ellipse.col=rgb(0,1,0,0.5)){
    ns <- nrow(x)
    for (i in 1:ns){
        uvwc <- UThHe2uvw.covmat(x,i)
        x0 <- uvwc$uvw[1]
        y0 <- uvwc$uvw[2]
        ell <- ellipse(x=x0,y=y0,covmat=uvwc$covmat,alpha=alpha)
        graphics::polygon(ell,col=ellipse.col)
        graphics::points(x0,y0,pch=19,cex=0.25)
    }
}

plot.logratio.contours <- function(x,contour.col=c('white','red'),
                                   xlim=NA,ylim=NA,...){
    R <- iratio('U238U235')[1]
    L8 <- lambda('U238')[1]
    L5 <- lambda('U235')[1]
    L2 <- lambda('Th232')[1]
    L7 <- lambda('Sm147')[1]
    uvw <- UThHe2uvw(x)
    f <- 0.5
    lims <- get.logratioplot.limits(uvw[,c('u','v')],f=f)
    if (all(!is.na(xlim))) lims[1:2] <- xlim
    if (all(!is.na(ylim))) lims[3:4] <- ylim
    graphics::plot(lims[1:2],lims[3:4],type='n',bty='n',
                   xlab='log[U/He]',ylab='log[Th/He]',...)
    tticks <- get.tticks(lims)
    res <- 500 # resolution
    du <- lims[2]-lims[1]
    dv <- lims[4]-lims[3]
    nt <- length(tticks)
    crp <- colorRampPalette(contour.col)(nt)
    for (i in 1:nt){
        tt <- tticks[i]
        a <- 8*R*(exp(L8*tt)-1)/(1+R) +
             7*(exp(L5*tt)-1)/(1+R)
        b <- 6*(exp(L2*tt)-1)
        uv.plot <- NULL
        for (j in 1:res){
            u <- lims[1]+du*j/res
            exp.v <- (1-a*exp(u))/b
            if (exp.v > exp(lims[3]) & exp.v < exp(lims[4])){
                v <- log(exp.v)
                uv.plot <- rbind(uv.plot,c(u,v))
            }
        }
        if (!is.null(uv.plot)){
            uv.first <- uv.plot[1,]
            uv.last <- uv.plot[nrow(uv.plot),]
            uv.plot <- rbind(uv.plot,c(uv.last[1],lims[3]))
            uv.plot <- rbind(uv.plot,lims[c(1,3)])
            uv.plot <- rbind(uv.plot,c(lims[1],uv.first[2]))
            uv.plot <- rbind(uv.plot,uv.first)
            graphics::polygon(uv.plot,col=crp[i])
            graphics::text(uv.first[1],uv.first[2],
                           labels=paste(tt,'Ma'),pos=4)
        }
    }
}

get.tticks <- function(lims){
    out <- pretty(get.UThHe.age.limits(lims))
    if (out[1]==0) {
        tticks <- pretty(out[1:2]*10)/10
        out[1] <- tticks[2]
    }
    out
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
    if (methods::is(x,'UThHe')) out <- log(x[,c('U','Th','Sm')])-log(x[,'He'])
    else out <- matrix(log(x[c('U','Th','Sm')])-log(x['He']),1,3)
    colnames(out) <- c('u','v','w')
    out
}

uvw2UThHe <- function(uvw,covmat){
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
    xy <- matrix(0,nrow=n,ncol=2)
    xy[,1] <- 0.5*(x+2*z)/(x+y+z)
    xy[,2] <- sin(pi/3)*x/(x+y+z)
    return(xy)
}
