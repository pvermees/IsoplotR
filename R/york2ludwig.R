# helper function for init_ludwig
york2ludwig <- function(x,anchor=0,type=0,model=1){
    if (x$format %in% c(1,2,3,9,10,119,1210)){
        yd <- x$x[,1:5]
    } else if (x$format %in% 4:8){
        stop("york2ludwig regression not yet implemented for this format")
    } else {
        stop("Invalid format or isochron type for york2ludwig regression")
    }
    out <- yfit <- york(yd)
    out$age <- ab2t(x=x,fit=yfit)[1,]
    names(out$age) <- c('t','s[t]')
    out$y0 <- c(1/yfit$a[1],yfit$a[2]/yfit$a[1]^2)
    names(out$y0) <- c('y','s[y]')
    out
}

ab2t <- function(x,fit){
    if (x$format%in%c(9,10,119,1210)){
        DP <- -fit$b[1]/fit$a[1]
        vDP <- errorprop1x2(J1=fit$b[1]/fit$a[1]^2,
                            J2=-1/fit$a[1],
                            E11=fit$a[2]^2,
                            E22=fit$b[2]^2,
                            E12=fit$cov.ab)
        sDP <- sqrt(vDP)
        if (x$format %in% c(9,119)){
            tst <- age(x=c(DP,sDP),method="U238-Pb206",d=x$d)
        } else {
            tst <- age(x=c(DP,sDP),method="U235-Pb207",d=x$d)
        }
    } else {
        stop("york2ludwig regression not yet implemented for this format")
    }
    tst
}
