# finds the (log of) the lower concordia intercept
# and the (log of) the common Pb intercept(s)
# plus the search ranges for ludwig regression
york2ludwig <- function(x,anchor=0,buffer=2){
    if (x$format<4){
        out <- york2ludwigTW(x=x,anchor=anchor,buffer=buffer)
    } else if (x$format<7){
        out <- york2ludwig204(x=x,anchor=anchor,buffer=buffer)
    } else if (x$format<9){
        out <- york2ludwig208(x=x,anchor=anchor,buffer=buffer)
    } else {
        stop("Invalid U-Pb format")
    }
    out
}

york2ludwigTW <- function(x,anchor=0,buffer=2){
    par <- lower <- upper <- vector()
    yd <- data2york(x,option=1)
    if (anchor[1]==1){
        Pb76c <- iratio('Pb207Pb206')[1]
        yfit <- MLyork(yd,anchor=c(1,Pb76c))
        tm <- WconcordiaIntersection(yfit=yfit,d=x$d)
        par['t'] <- log(tm[1])
        upper['t'] <- log(tm[2])
        if (iratio('Pb207Pb206')[2]>0){
            par['a0'] <- log(Pb76c)
            lower['a0'] <- log(age_to_Pb207Pb206_ratio(tt=tm[2],d=x$d)[1])
        }
    } else if (anchor[1]==2 & length(anchor)>1){
        tt <- anchor[2]
        if ((length(anchor)>2 && anchor[3]>0)){
            par['t'] <- log(tt)
            upper['t'] <- par['t'] + buffer
        }
        McL <- mclean(tt,d=x$d)
        Xt <- McL$Pb207U235
        Yt <- McL$Pb207U235
        YD <- yd
        YD[,'X'] <- yd[,'X'] - Xt # shift left
        yfit <- MLyork(YD,anchor=c(1,Yt))
        Pb76c <- 1/unname(yfit$b[1]*iratio('U238U235')[1])
        par['a0'] <- log(Pb76c)
        lower['a0'] <- log(age_to_Pb207Pb206_ratio(tt=tt,d=x$d)[1])
    } else { # no anchor
        yfit <- york(yd)
        tm <- WconcordiaIntersection(yfit=yfit,d=x$d)
        par['t'] <- log(tm[1])
        upper['t'] <- log(tm[2])
        Pb76c <- 1/unname(yfit$b[1]*iratio('U238U235')[1])
        par['a0'] <- log(Pb76c)
        lower['a0'] <- log(age_to_Pb207Pb206_ratio(tt=tm[2],d=x$d)[1])
    }
    pnames <- names(par)
    if ('t' %in% pnames) lower['t'] <- par['t'] - buffer
    if ('a0' %in% pnames) upper['a0'] <- par['a0'] + buffer
    list(par=par,lower=lower,upper=upper)
}

york2ludwig204 <- function(x,anchor=0,type=0,buffer=2){
    par <- lower <- upper <- vector()
    yda <- data2york(x,option=3)
    ydb <- data2york(x,option=4)
    if (type==1) yd <- yda
    if (type==2) yd <- ydb
    if (anchor[1]==1){
        if (type%in%c('joint',0,1)){
            Pb64c <- iratio('Pb206Pb204')[1]
            abxa <- inithelper(yd=yda,y0=1/Pb64c)
            tt <- get.Pb206U238.age(x=abxa['x0inv'],d=x$d)[1]
            par['t'] <- log(tt)
            if (iratio('Pb206Pb204')[2]>0) par['a0'] <- log(Pb64c)
        }
        if (type%in%c('joint',0,2)){
            Pb74c <- iratio('Pb207Pb204')[1]
            abxb <- inithelper(yd=ydb,y0=1/Pb74c)
        }
        if (type==2){
            tt <- get.Pb206U238.age(x=abxb['x0inv'],d=x$d)[1]
            par['t'] <- log(tt)
        }
        if (type%in%c('joint',0,2)){
            if (iratio('Pb207Pb204')[2]>0) par['b0'] <- log(Pb74c)
        }
    } else if (anchor[1]==2 && length(anchor)>1){
        tt <- anchor[2]
        if ((length(anchor)>2 && anchor[3]>0)) par['t'] <- log(tt)
        if (type%in%c('joint',0,1)){
            Pb6U8r <- mclean(tt=tt)$Pb206U238
            abxa <- inithelper(yd=yda,x0=1/Pb6U8r)
            par['a0'] <- log(1/abxa['a'])
        }
        if (type%in%c('joint',0,2)){
            Pb7U5r <- mclean(tt=tt)$Pb207U235
            abxb <- inithelper(yd=ydb,x0=1/Pb7U5r)
            par['b0'] <- log(1/abxb['a'])
        }
    } else {
        if (type%in%c('joint',0,1)){
            abxa <- inithelper(yd=yda)
            tt <- get.Pb206U238.age(x=abxa['x0inv'],d=x$d)[1]
        }
        if (type%in%c('joint',0,2)){
            abxb <- inithelper(yd=ydb)
        }
        if (type==2){
            tt <- get.Pb207U235.age(x=abxb['x0inv'],d=x$d)[1]
        }
        par['t'] <- log(tt)
        if (type%in%c('joint',0,1)){
            par['a0'] <- log(1/abxa['a'])
        }
        if (type%in%c('joint',0,2)){
            par['b0'] <- log(1/abxb['a'])
        }
    }
    list(par=par,lower=par-buffer,upper=par+buffer)
}

york2ludwig208 <- function(x,anchor=0,type=0,buffer=2){
    par <- lower <- upper <- vector()
    if (anchor[1]==1){
        yd <- data2york(x,option=2)
        abx <- inithelper(yd=yd)
        pilott <- get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1]
        if (type==1){ # 0806 vs 38/06
            yd <- data2york(x,option=6,tt=pilott)
            y0 <- iratio('Pb208Pb206')[1]
        } else if (type==2){ # 0807 vs 35/07
            yd <- data2york(x,option=7,tt=pilott)
            y0 <- iratio('Pb208Pb207')[1]
        } else if (type==3){ # 0608 vs 32/08
            yd <- data2york(x,option=8,tt=pilott)
            y0 <- 1/iratio('Pb208Pb206')[1]
        } else if (type==4){ # 0708 vs 32/08
            yd <- data2york(x,option=9,tt=pilott)
            y0 <- 1/iratio('Pb208Pb207')[1]
        } else { # joint, 0 or 1
            y0 <- iratio('Pb207Pb206')[1]
        }
        abx <- inithelper(yd=yd,y0=y0)
        if (type==1){
            par['t'] <- log(get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1])
            if (iratio('Pb208Pb206')[2]>0) par['a0'] <- log(y0)
        } else if (type==2){
            par['t'] <- log(get.Pb207U235.age(x=abx['x0inv'],d=x$d)[1])
            if (iratio('Pb208Pb207')[2]>0) par['b0'] <- log(y0)
        } else if (type==3){
            par['t'] <- log(get.Pb208Th232.age(x=abx['x0inv'],d=x$d)[1])
            if (iratio('Pb208Pb206')[2]>0) par['a0'] <- log(1/y0)
        } else if (type==4){
            par['t'] <- log(get.Pb208Th232.age(x=abx['x0inv'],d=x$d)[1])
            if (iratio('Pb208Pb207')[2]>0) par['b0'] <- log(1/y0)
        } else { # joint, 0 or 1
            par['t'] <- log(get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1])
            if (iratio('Pb208Pb206')[2]>0){
                par['a0'] <- log(iratio('Pb208Pb206')[1])
            }
            if (iratio('Pb208Pb207')[2]>0){
                par['b0'] <- log(iratio('Pb208Pb207')[1])
            }
        }
    } else if (anchor[1]==2 && length(anchor)>1){
        tt <- anchor[2]
        if (length(anchor)>2 && anchor[3]>0) par['t'] <- log(tt)
        if (type==2){ # 0807 vs 35/07
            yd <- data2york(x,option=7,tt=tt)
            x0 <- age_to_U235Pb207_ratio(tt)[1]
        } else if (type==3){ # 0608 vs 32/08
            yd <- data2york(x,option=8,tt=tt)
            x0 <- 1/age_to_Pb208Th232_ratio(tt)[1]
        } else if (type==4){ # 0708 vs 32/08
            yd <- data2york(x,option=9,tt=tt)
            x0 <- 1/age_to_Pb208Th232_ratio(tt)[1]
        } else { # joint, 0 or 1: 0806 vs 38/06
            yd <- data2york(x,option=6,tt=tt)
            x0 <- age_to_U238Pb206_ratio(tt)[1]
        }
        abx <- inithelper(yd=yd,x0=x0)
        if (type==1){
            par['a0'] <- log(1/abx['a'])
        } else if (type==2){
            par['b0'] <- log(1/abx['a'])
        } else if (type==3){
            par['a0'] <- log(abx['a'])
        } else if (type==4){
            par['b0'] <- log(abx['a'])
        } else {
            par['a0'] <- log(1/abx['a'])
            ydb <- data2york(x,option=7)
            abxb <- inithelper(yd=ydb,x0=x0)
            par['b0'] <- log(1/abxb['a'])
        }  
    } else {
        yd <- data2york(x,option=2)
        abx <- inithelper(yd=yd)
        pilott <- get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1]
        if (type==1){ # 0806 vs 38/06
            yd <- data2york(x,option=6,tt=pilott)
        } else if (type==2){ # 0807 vs 35/07
            yd <- data2york(x,option=7,tt=pilott)
        } else if (type==3){ # 0608 vs 32/08
            yd <- data2york(x,option=8,tt=pilott)
        } else if (type==4){ # 0708 vs 32/08
            yd <- data2york(x,option=9,tt=pilott)
        } else { # joint, 0 or 1
                                        # keep yd
        }
        abx <- inithelper(yd=yd)
        par['t'] <- log(pilott)
        if (type==1){ # 0806 vs 38/06
            par['a0'] <- log(1/abx['a'])
        } else if (type==2){ # 0807 vs 35/07
            par['b0'] <- log(1/abx['a'])
        } else if (type==3){ # 0608 vs 32/08
            par['a0'] <- log(abx['a'])
        } else if (type==4){ # 0708 vs 32/08
            par['b0'] <- log(abx['a'])
        } else { # joint, 0 or 1
            par['t'] <- log(get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1])
            yda <- data2york(x,option=6)
            ydb <- data2york(x,option=7)
            abxa <- inithelper(yd=yda)
            abxb <- inithelper(yd=ydb)
            par['a0'] <- log(1/abxa['a'])
            par['b0'] <- log(1/abxb['a'])
        }            
    }
    list(par=par,lower=par-buffer,upper=par+buffer)
}

WconcordiaIntersection <- function(yfit,d=diseq()){
    misfit <- function(tt,a,b,d,gradient=FALSE){
        McL <- mclean(tt=tt,d=d)
        if (gradient){
            dXydt <- McL$dPb207U235dt
            dYwdt <- McL$dPb206U238dt
            out <- dYwdt - b*dXydt
        } else {
            Xy <- McL$Pb207U235 # York coordinate
            Yw <- McL$Pb206U238 # Wetherill coordinate
            out <- Yw - a - b*Xy
        }
        out
    }
    a <- unname(abs(yfit$a[1]))
    b <- unname(yfit$b[1])
    midpoint <- stats::uniroot(misfit,lower=0,upper=4600,
                               a=a,b=b,d=d,gradient=TRUE)$root
    if (misfit(tt=midpoint,a=a,b=b,d=d) > 0){
        tt <- stats::uniroot(misfit,lower=0,upper=midpoint,a=a,b=b,d=d)$root
    } else {
        tt <- stats::optimise(misfit,lower=0,upper=4600,a=a,b=b,d=d)$minimum
    }
    c(tt,midpoint)
}

inithelper <- function(yd,x0=NULL,y0=NULL){
    out <- c('a'=NA,'b'=NA,'x0inv'=NA)
    if (is.null(x0) & is.null(y0)){
        fit <- york(yd)
        out['a'] <- max(sqrt(.Machine$double.eps),fit$a[1])
        out['b'] <- min(-sqrt(.Machine$double.eps),fit$b[1])
    } else if (is.null(y0)){ # anchor x
        i <- which.min(yd[,'X'])
        if (x0>yd[i,'X']){
            out['b'] <- yd[i,'Y']/(yd[i,'X']-x0)
        } else {
            out['b'] <- yd[i,'Y']/(0-x0)
        }
        if (is.finite(x0)){
            out['a'] <- abs(out['b']*x0)
        } else {
            out['a'] <- mean(yd[,'Y'])
        }
    } else { # anchor y
        i <- which.min(yd[,'Y'])
        if (y0>yd[i,'Y']){
            out['b'] <- (yd[i,'Y']-y0)/yd[i,'X']
        } else {
            out['b'] <- y0/(0-yd[i,'X'])
        }
        out['a'] <- y0
    }
    out['x0inv'] <- -out['b']/out['a']
    out
}
