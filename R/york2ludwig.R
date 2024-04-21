# helper function for init_ludwig
york2ludwig <- function(x,anchor=0,buffer=2,type=0,model=1){
    if (anchor[1]==3){ # stacey-kramers
        init <- york2ludwig(x=x,anchor=0,buffer=buffer,type=type,model=model)
        out <- list(par=init$par['t'],
                    lower=init$lower['t'],
                    upper=init$upper['t'])
    } else if (x$format<4){
        out <- york2ludwigTW(x=x,anchor=anchor,buffer=buffer,model=model)
    } else if (x$format%in%c(4,5,6,9,10,85,119,1210)){
        out <- york2ludwig20x(x=x,anchor=anchor,buffer=buffer,type=type,model=model)
    } else if (x$format%in%c(7,8,11,12)){
        out <- york2ludwig208(x=x,anchor=anchor,buffer=buffer,type=type,model=model)
    } else {
        stop("Invalid U-Pb format")
    }
    out
}

york2ludwigTW <- function(x,anchor=0,buffer=2,model=1){
    par <- lower <- upper <- c()
    yd <- data2york(x,option=1)
    if (anchor[1]==1){
        Pb76c <- iratio('Pb207Pb206')
        U85 <- iratio('U238U235')[1]
        yfit <- MLyork(yd,anchor=c(2,1/(Pb76c[1]*U85)))
        tm <- WconcordiaIntersection(yfit=yfit,d=x$d)
        par['t'] <- log(tm[1])
        lower['t'] <- par['t'] - buffer
        upper['t'] <- log(tm[2])
        if (model==1 & Pb76c[2]>0){
            par['a0'] <- log(Pb76c[1])
            lower['a0'] <- log(age_to_Pb207Pb206_ratio(tt=tm[2],d=x$d)[1])
            upper['a0'] <- par['a0'] + buffer
        } else if (model==3 & Pb76c[2]<=0){
            Pb76err <- data2york(x,option=2)[,'sY']
            par['w'] <- log(stats::median(Pb76err))
            lower['w'] <- par['w'] - max(buffer,10)
            upper['w'] <- log(Pb76c[1])
        }
    } else if (anchor[1]==2 & length(anchor)>1){
        tt <- anchor[2]
        if (model==1 & length(anchor)>2 && anchor[3]>0){
            par['t'] <- log(tt)
            lower['t'] <- par['t'] - buffer
            upper['t'] <- par['t'] + buffer
        }
        McL <- mclean(tt,d=x$d)
        Xt <- McL$Pb207U235
        Yt <- McL$Pb206U238
        YD <- yd
        YD[,'X'] <- yd[,'X'] - Xt # shift left
        yfit <- MLyork(YD,anchor=c(1,Yt))
        Pb76c <- 1/unname(yfit$b[1]*iratio('U238U235')[1])
        par['a0'] <- log(Pb76c)
        lower['a0'] <- log(age_to_Pb207Pb206_ratio(tt=tt,d=x$d)[1])
        upper['a0'] <- par['a0'] + buffer
        if (model==3 & length(anchor)>2 && anchor[3]<=0){
            stPb68 <- get_Pb206U238_age(x)[,2]
            par['w'] <- log(stats::median(stPb68))
            lower['w'] <- par['w'] - max(buffer,10)
            upper['w'] <- log(tt) + buffer
        }
    } else { # no anchor
        yfit <- york(yd)
        tm <- WconcordiaIntersection(yfit=yfit,d=x$d)
        par['t'] <- log(tm[1])
        lower['t'] <- par['t'] - buffer
        upper['t'] <- log(tm[2])
        Pb76c <- 1/unname(yfit$b[1]*iratio('U238U235')[1])
        par['a0'] <- log(Pb76c)
        lower['a0'] <- log(age_to_Pb207Pb206_ratio(tt=tm[2],d=x$d)[1])
        upper['a0'] <- par['a0'] + buffer
        if (model==3){
            stPb68 <- get_Pb206U238_age(x)[,2]
            par['w'] <- log(stats::median(stPb68))
            lower['w'] <- par['w'] - max(buffer,10)
            upper['w'] <- par['t'] + buffer
        }
    }
    list(par=par,lower=lower,upper=upper)
}

york2ludwig20x <- function(x,anchor=0,type=0,buffer=2,model=1){
    par <- lower <- upper <- vector()
    if (x$format%in%c(4,5,6,9,85,119)) yda <- data2york(x,option=3)
    if (x$format%in%c(4,5,6,10,85,1210)) ydb <- data2york(x,option=4)
    if (type==1) yd <- yda
    if (type==2) yd <- ydb
    if (anchor[1]==1){
        if (type%in%c('joint',0,1)){
            if (x$format%in%c(85,119)){
                Pb6xc <- iratio('Pb206Pb208')
            } else {
                Pb6xc <- iratio('Pb206Pb204')
            }
            abxa <- inithelper(yd=yda,y0=1/Pb6xc[1])
            tt <- get_Pb206U238_age(x=abxa['x0inv'],d=x$d)[1]
            par['t'] <- log(tt)
            if (model==1 & Pb6xc[2]>0){
                par['a0'] <- log(Pb6xc[1])
            }
        }
        if (type%in%c('joint',0,2)){
            if (x$format%in%c(85,1210)){
                Pb7xc <- iratio('Pb207Pb208')
            } else {
                Pb7xc <- iratio('Pb207Pb204')
            }
            abxb <- inithelper(yd=ydb,y0=1/Pb7xc[1])
        }
        if (type==2){
            tt <- get_Pb207U235_age(x=abxb['x0inv'],d=x$d)[1]
            par['t'] <- log(tt)
        }
        if (model==1 & type%in%c('joint',0,2) && Pb7xc[2]>0){
            par['b0'] <- log(Pb7xc[1])
        }
        if (model==3){
            if (type==1 && Pb6xc[2]<=0){
                par['w'] <- log(stats::median(yda[,'sY']))
                upper['w'] <- log(Pb6xc[1])
            } else if (type==2 && Pb7xc[2]<=0){
                par['w'] <- log(stats::median(ydb[,'sY']))
                upper['w'] <- log(Pb7xc[1])
            }
        }
    } else if (anchor[1]==2 & length(anchor)>1){
        tt <- anchor[2]
        if (model==1 & length(anchor)>2 & anchor[3]>0){
            par['t'] <- log(tt)
        }
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
        if (model==3 & length(anchor)>2 && anchor[3]<=0){
            if (type==1){
                stPb6U8 <- get_Pb206U238_age(x)[,2]
                par['w'] <- log(stats::median(stPb6U8))
            } else if (type==2){
                stPb7U5 <- get_Pb207U235_age(x)[,2]
                par['w'] <- log(stats::median(stPb7U5))
            }
            upper['w'] <- log(tt)
        }
    } else {
        if (type%in%c('joint',0,1)){
            abxa <- inithelper(yd=yda)
            tt <- get_Pb206U238_age(x=abxa['x0inv'],d=x$d)[1]
        }
        if (type%in%c('joint',0,2)){
            abxb <- inithelper(yd=ydb)
        }
        if (type==2){
            tt <- get_Pb207U235_age(x=abxb['x0inv'],d=x$d)[1]
        }
        par['t'] <- log(tt)
        if (type%in%c('joint',0,1)){
            par['a0'] <- log(1/abxa['a'])
        }
        if (type%in%c('joint',0,2)){
            par['b0'] <- log(1/abxb['a'])
        }
        if (model==3){
            if (type%in%c('joint',0,1)){
                par['w'] <- log(stats::median(get_Pb206U238_age(x=x)[,2]))
            } else {
                par['w'] <- log(stats::median(get_Pb207U235_age(x=x)[,2]))
            }
            upper['w']  <- log(tt)
        }
    }
    for (pname in c('t','a0','b0')){
        if (pname%in%names(par)){
            lower[pname] <- par[pname] - buffer
            upper[pname] <- par[pname] + buffer
        }
    }
    if ('w'%in%names(par)){
        lower['w'] <- par['w'] - max(buffer,10)
    }
    list(par=par,lower=lower,upper=upper)
}

york2ludwig208 <- function(x,anchor=0,type=0,buffer=2,model=1){
    par <- lower <- upper <- vector()
    if (anchor[1]==1){
        if (type==1){ # 0806 vs 38/06
            Pb68c <- iratio('Pb206Pb208')
            pilott <- min(get_Pb206U238_age(x=x)[,1])
            yd <- data2york(x,option=6,tt=pilott)
            y0 <- 1/Pb68c[1]
        } else if (type==2){ # 0807 vs 35/07
            Pb78c <- iratio('Pb207Pb208')
            pilott <- min(get_Pb207U235_age(x=x)[,1])
            yd <- data2york(x,option=7,tt=pilott)
            y0 <- 1/Pb78c[1]
        } else if (type==3){ # 0608 vs 32/08
            Pb68c <- iratio('Pb206Pb208')
            pilott <- min(get_Pb208Th232_age(x=x)[,1])
            yd <- data2york(x,option=8,tt=pilott)
            y0 <- Pb68c[1]
        } else if (type==4){ # 0708 vs 32/08
            Pb78c <- iratio('Pb207Pb208')
            pilott <- min(get_Pb208Th232_age(x=x)[,1])
            yd <- data2york(x,option=9,tt=pilott)
            y0 <- Pb78c[1]
        } else { # joint, 0 or 1
            yd <- data2york(x,option=2)
            pilott <- min(get_Pb206U238_age(x=x)[,1])
            y0 <- iratio('Pb207Pb206')[1]
        }
        abx <- inithelper(yd=yd,y0=y0)
        if (type==1){
            par['t'] <- log(get_Pb206U238_age(x=abx['x0inv'],d=x$d)[1])
            if (model==1 & Pb68c[2]>0){
                par['a0'] <- -log(y0)
            } else if (model==3 & Pb68c[2]<=0){
                par['w'] <- stats::median(yd[,'sY']/yd[,'Y']^2)
                upper['w'] <- log(Pb68c[1])
            }
        } else if (type==2){
            par['t'] <- log(get_Pb207U235_age(x=abx['x0inv'],d=x$d)[1])
            if (model==1 & Pb78c[2]>0){
                par['b0'] <- -log(y0)
            } else if (model==3 & Pb78c[2]<=0){
                par['w'] <- stats::median(yd[,'sY']/yd[,'Y']^2)
                upper['w'] <- log(Pb78c[1])
            }
        } else if (type==3){
            par['t'] <- log(get_Pb208Th232_age(x=abx['x0inv'],d=x$d)[1])
            if (model==1 & Pb68c[2]>0){
                par['a0'] <- log(y0)
            } else if (model==3 & Pb68c[2]<=0){
                par['w'] <- stats::median(yd[,'sY'])
                upper['w'] <- log(Pb68c[1])
            }
        } else if (type==4){
            par['t'] <- log(get_Pb208Th232_age(x=abx['x0inv'],d=x$d)[1])
            if (model==1 & Pb78c[2]>0){
                par['b0'] <- log(y0)
            } else if (model==3 & Pb78c[2]<=0){
                par['w'] <- stats::median(yd[,'sY'])
                upper['w'] <- log(Pb78c[1])
            }
        } else { # joint, 0 or 1
            par['t'] <- log(get_Pb206U238_age(x=abx['x0inv'],d=x$d)[1])
            Pb68c <- iratio('Pb206Pb208')
            if (model==1 & Pb68c[2]>0){
                par['a0'] <- -log(Pb68c[1])
            }
            Pb78c <- iratio('Pb207Pb208')
            if (model==1 & Pb78c[2]>0){
                par['b0'] <- -log(Pb78c[1])
            }
        }
    } else if (anchor[1]==2 & length(anchor)>1){
        tt <- anchor[2]
        if (model==1 & length(anchor)>2 & anchor[3]>0){
            par['t'] <- log(tt)
        }
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
        if (model==3 & length(anchor)>2 && anchor[3]<=0){
            if (type%in%c('joint',0,1)){
                stPb6U8 <- get_Pb206U238_age(x)[,2]
                par['w'] <- log(stats::median(stPb6U8))
            } else if (type==2){
                stPb7U5 <- get_Pb207U235_age(x)[,2]
                par['w'] <- log(stats::median(stPb7U5))
            } else if (type%in%c(3,4)){
                stPb8Th2 <- get_Pb208Th232_age(x)[,2]
                par['w'] <- log(stats::median(stPb8Th2))
            } else {
                stop('Invalid isochron type.')
            }
            upper['w'] <- log(tt)
        }        
    } else {
        if (type==2){ # 0807 vs 35/07
            pilott <- min(get_Pb207U235_age(x)[,1])
            yd <- data2york(x,option=7,tt=pilott)
        } else if (type==3){ # 0608 vs 32/08
            pilott <- min(get_Pb208Th232_age(x)[,1])
            yd <- data2york(x,option=8,tt=pilott)
        } else if (type==4){ # 0708 vs 32/08
            pilott <- min(get_Pb208Th232_age(x)[,1])
            yd <- data2york(x,option=9,tt=pilott)
        } else { # joint, 0 or 1 (0806 vs 38/06)
            pilott <- min(get_Pb206U238_age(x)[,1])
            yd <- data2york(x,option=6,tt=pilott)
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
        } else { # joint or 0
            par['t'] <- log(get_Pb206U238_age(x=abx['x0inv'],d=x$d)[1])
            yda <- data2york(x,option=6)
            ydb <- data2york(x,option=7)
            abxa <- inithelper(yd=yda)
            abxb <- inithelper(yd=ydb)
            par['a0'] <- log(1/abxa['a'])
            par['b0'] <- log(1/abxb['a'])
        }            
        if (model==3){
            if (type%in%c('joint',0,1)){
                par['w'] <- log(stats::median(get_Pb206U238_age(x=x)[,2]))
            } else if (type==2){
                par['w'] <- log(stats::median(get_Pb207U235_age(x=x)[,2]))
            } else if (type%in%c(3,4)){
                par['w'] <- log(stats::median(get_Pb208Th232_age(x=x)[,2]))
            } else {
                stop('Invalid isochron type.')
            }
            upper['w']  <- log(pilott)
        }
    }
    for (pname in c('t','a0','b0')){
        if (pname%in%names(par)){
            lower[pname] <- par[pname] - buffer
            upper[pname] <- par[pname] + buffer
        }
    }
    if ('w'%in%names(par)){
        lower['w'] <- par['w'] - max(buffer,10)
    }
    list(par=par,lower=lower,upper=upper)
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
