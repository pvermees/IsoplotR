# helper function for init_ludwig
york2ludwig <- function(x){
    if (x$format<4){
        yd <- data2york(x=x,option=2)
        out <- TWconcordiaIntersection(yfit=york(yd),d=x$d)
        type <- 0
    } else if (x$format%in%c(9,119)){
        yd <- x$x[,1:5]
        out <- Pb6U8isochron2ludwig(yfit=york(yd),d=x$d)
        type <- 1
    } else if (x$format%in%c(10,1210)){
        yd <- x$x[,1:5]
        out <- Pb7U5isochron2ludwig(yfit=york(yd),d=x$d)
        type <- 2
    } else {
        stop("Invalid format or type for york2ludwig")
    }
    LLpar <- c(log(out$par[1:2]),out$par[-(1:2)])
    out$value <- LL_ludwig(par=LLpar,x=x,type=type)
    out$n <- nrow(yd)
    out
}

TWconcordiaIntersection <- function(yfit,d=diseq()){
    misfit <- function(tt,yfit,d=diseq()){
        McL <- mclean(tt,d=d)
        yfit$a[1] + yfit$b[1]/McL$Pb206U238 - McL$Pb207Pb206
    }
    parcov <- function(tt,yfit,d=diseq){
        out <- list()
        McL <- mclean(tt,d=d)
        dfda <- 1
        dfdb <- 1/McL$Pb206U238
        dfdt <- -yfit$b[1]*McL$dPb206U238dt/McL$Pb206U238^2 - McL$dPb207Pb206dt
        dtda <- -dfda/dfdt
        dtdb <- -dfdb/dfdt
        out$par <- c(tt,yfit$a[1])
        pnames <- c('t','a0')
        J <- matrix(0,2,2)
        J[1,1] <- dtda
        J[1,2] <- dtdb
        J[2,1] <- 1 # da0da
        if (d$U48$option==2){
            out$par <- append(out$par,McL$U48i)
            pnames <- append(pnames,'U48i')
            dU48ida <- McL$dU48idt*dtda
            dU48idb <- McL$dU48idt*dtdb
            J <- rbind(J,c(dU48ida,dU48idb))
        }
        if (d$ThU$option==2){
            out$par <- append(out$par,McL$ThUi)
            pnames <- append(pnames,'ThUi')
            dThUida <- McL$dThUidt*dtda
            dThUidb <- McL$dThUidt*dtdb
            J <- rbind(J,c(dThUida,dThUidb))
        }
        E <- rbind(c(yfit$a[2]^2,yfit$cov.ab),
                   c(yfit$cov.ab,yfit$b[2]^2))
        out$cov <- J %*% E %*% t(J)
        names(out$par) <- rownames(out$cov) <- colnames(out$cov) <- pnames
        out
    }
    out <- yfit
    if (yfit$b[1] > 0){ # one root
        ll <- age(x=yfit$a,method="Pb207-Pb206")[1]
        t1 <- uniroot(f=misfit,lower=ll,upper=2*ll,yfit=yfit,d=d,extendInt="yes")$root
        t2 <- Inf
        out <- append(out,parcov(t1,yfit))
    } else { # two or no roots
        init <- find_concordia_age_parallel_to_line(yfit=yfit,d=d)
        McL <- mclean(tt=init$midpoint)
        x_midpoint <- 1/McL$Pb206U238
        y_midpoint <- 1/McL$Pb207Pb206
        if (y_midpoint > yfit$a[1] + yfit$b[1]*x_midpoint){ # two roots
            t1 <- uniroot(f=misfit,lower=init$ll,upper=init$midpoint,yfit=yfit,d=d)$root
            t2 <- uniroot(f=misfit,lower=init$midpoint,upper=init$ul,yfit=yfit)$root
        } else { # no roots
            stop('No concordia intersections')
        }
        out <- append(out,parcov(t1,yfit,d))
    }
    out$t <- setNames(c(t1,t2),c('t[l]','t[u]'))
    out
}

find_concordia_age_parallel_to_line <- function(yfit,d=diseq()){
    misfit <- function(tt,b,d=diseq()){
        McL <- mclean(tt,d=d)
        dU238Pb206dt <- -McL$dPb206U238dt/McL$Pb206U238^2
        slope <- McL$dPb207Pb206dt/dU238Pb206dt
        (slope - b)^2
    }
    ll <- age(x=-yfit$b[1]/yfit$a[1],method="U238-Pb206",d=d)[1]
    ul <- age(x=yfit$a[1],method="Pb207-Pb206")[1]
    midpoint <- stats::optim(ul,fn=misfit,method="BFGS",b=yfit$b[1])$par
    list(ll=ll,midpoint=midpoint,ul=ul)
}

Pb6U8isochron2ludwig <- function(yfit,d=diseq()){
    out <- yfit
    DP <- -yfit$b[1]/yfit$a[1]
    t1 <- age(x=DP,method="U238-Pb206",d=d)[1]
    pnames <- c('t','a0')
    McL <- mclean(t1,d=d)
    np <- 2
    out$par <- c(t1,1/yfit$a[1])
    if (d$U48$option==2) np <- np + 1
    if (d$ThU$option==2) np <- np + 1
    J <- matrix(0,np,2)
    E <- rbind(c(yfit$a[2]^2,yfit$cov.ab),
               c(yfit$cov.ab,yfit$b[2]^2))
    dtdDP <- 1/McL$dPb206U238dt
    dDPda <- yfit$b[1]/yfit$a[1]^2
    dDPdb <- -1/yfit$a[1]
    dtda <- dtdDP*dDPda
    dtdb <- dtdDP*dDPdb
    J[1,1] <- dtda
    J[1,2] <- dtdb
    J[2,1] <- -1/yfit$a[1]^2 # da0da
    if (d$U48$option==2){
        J[3,1] <- McL$dU48idt*dtda # dU48ida
        J[3,2] <- McL$dU48idt*dtdb # dU48idb
        out$par <- append(out$par,McL$U48i)
        pnames <- append(pnames,'U48i')
    }
    if (d$ThU$option==2){
        J[np,1] <- McL$dThUidt*dtda # dThUida
        J[np,2] <- McL$dThUidt*dtdb # dThUidb
        out$par <- append(out$par,McL$ThUi)
        pnames <- append(pnames,'ThUi')
    }
    out$cov <- J %*% E %*% t(J)
    names(out$par) <- rownames(out$cov) <- colnames(out$cov) <- pnames
    out
}

Pb7U5isochron2ludwig <- function(yfit,d=diseq()){
    out <- yfit
    DP <- -yfit$b[1]/yfit$a[1]
    t1 <- age(x=DP,method="U235-Pb207",d=d)[1]
    pnames <- c('t','b0')
    McL <- mclean(t1,d=d)
    out$par <- c(t1,1/yfit$a[1])
    J <- matrix(0,2,2)
    E <- rbind(c(yfit$a[2]^2,yfit$cov.ab),
               c(yfit$cov.ab,yfit$b[2]^2))
    dtdDP <- 1/McL$dPb207U235dt
    dDPda <- yfit$b[1]/yfit$a[1]^2
    dDPdb <- -1/yfit$a[1]
    dtda <- dtdDP*dDPda
    dtdb <- dtdDP*dDPdb
    J[1,1] <- dtda
    J[1,2] <- dtdb
    J[2,1] <- -1/yfit$a[1]^2
    out$cov <- J %*% E %*% t(J)
    names(out$par) <- rownames(out$cov) <- colnames(out$cov) <- pnames
    out
}
