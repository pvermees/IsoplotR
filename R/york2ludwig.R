# helper function for init_ludwig
york2ludwig <- function(x,type=2){
    if (x$format<4){
        yd <- data2york(x=x,option=2)
    } else if (x$format%in%c(9,10,119,1210)){
        yd <- x$x[,1:5]
    } else {
        stop("Invalid format or type for york2ludwig")
    }
    yfit <- york(yd)
    out <- ab2t(x=x,yfit=yfit)
    out$y0 <- c(1/yfit$a[1],yfit$a[2]/yfit$a[1]^2)
    names(out$y0) <- c('y','s[y]')
    out
}

ab2t <- function(x,yfit){
    if (x$format<4){
        out <- TWconcordiaIntersection(yfit=yfit,d=x$d)
    } else if (x$format%in%c(9,10,119,1210)){
        DP <- -fit$b[1]/fit$a[1]
        vDP <- errorprop1x2(J1=fit$b[1]/fit$a[1]^2,
                            J2=-1/fit$a[1],
                            E11=fit$a[2]^2,
                            E22=fit$b[2]^2,
                            E12=fit$cov.ab)
        sDP <- sqrt(vDP)
        out <- yfit
        if (x$format %in% c(9,119)){
            out$t1 <- age(x=c(DP,sDP),method="U238-Pb206",d=x$d)
        } else {
            out$t1 <- age(x=c(DP,sDP),method="U235-Pb207",d=x$d)
        }
    } else {
        stop("york2ludwig regression not yet implemented for this format")
    }
    out
}

TWconcordiaIntersection <- function(yfit,d=diseq()){
    out <- yfit
    misfit <- function(tt,yfit,d=diseq()){
        McL <- mclean(tt,d=d)
        yfit$a[1] + yfit$b[1]/McL$Pb206U238 - McL$Pb207Pb206
    }
    error <- function(tt,yfit,d=diseq){
        McL <- mclean(tt,d=d)
        dfda <- 1
        dfdb <- 1/McL$Pb206U238
        dfdt <- -yfit$b[1]*McL$dPb206U238dt/McL$Pb206U238^2 - McL$dPb207Pb206dt
        dtda <- -dfda/dfdt
        dtdb <- -dfdb/dfdt
        vt <- errorprop1x2(J1=dtda,J2=dtdb,E11=yfit$a[2]^2,E22=yfit$b[2]^2,E12=yfit$cov.ab)
        sqrt(vt)
    }
    if (yfit$b[1] > 0){ # one root
        ll <- age(x=yfit$a[1],method="Pb207-Pb206",d=d)
        t1 <- uniroot(f=misfit,lower=ll,upper=2*ll,yfit=yfit,d=d,extendInt="yes")$root
        st1 <- error(tt=t1,yfit=yfit,d=d)
        t2 <- st2 <- NULL
    } else { # two or no roots
        init <- find_concordia_age_parallel_to_line(yfit=yfit,d=d)
        McL <- mclean(tt=init$midpoint,d=d)
        x_midpoint <- 1/McL$Pb206U238
        y_midpoint <- 1/McL$Pb207Pb206
        if (y_midpoint > yfit$a[1] + yfit$b[1]*x_midpoint){ # two roots
            t1 <- uniroot(f=misfit,lower=init$ll,upper=init$midpoint,yfit=yfit,d=d)$root
            t2 <- uniroot(f=misfit,lower=init$midpoint,upper=init$ul,yfit=yfit,d=d)$root
            st1 <- error(tt=t1,yfit=yfit,d=d)
            st2 <- error(tt=t2,yfit=yfit,d=d)
        } else { # no roots
            t1 <- t2 <- NULL
        }
    }
    out$t1 <- c(t1,st1)
    out$t2 <- c(t2,st2)
    out
}

find_concordia_age_parallel_to_line <- function(yfit,d=diseq()){
    misfit <- function(tt,b,d=diseq()){
        McL <- mclean(tt,d=d)
        dU238Pb206dt <- -McL$dPb206U238dt/McL$Pb206U238^2
        slope <- McL$dPb207Pb206dt/dU238Pb206dt
        slope - b
    }
    ll <- age(x=-yfit$b[1]/yfit$a[1],method="U238-Pb206",d=d)[1]
    ul <- age(x=yfit$a[1],method="Pb207-Pb206")[1]
    midpoint <- uniroot(f=misfit,lower=ll,upper=ul,b=yfit$b[1],d=d)$root
    list(ll=ll,midpoint=midpoint,ul=ul)
}
