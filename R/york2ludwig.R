york2ludwig <- function(x,anchor=0){
    yd <- data2york(x,option=1)
    if (anchor[1]==1){
        Pb76i <- iratio('Pb207Pb206')[1]
        yfit <- MLyork(yd,anchor=c(1,y0))
        tt <- concordiaIntersection(yfit=yfit,d=x$d)
    } else if (anchor[1]==2 & length(anchor)>1){
        tt <- anchor[2]
        Pb7U5 <- 
        yfit <- york(yd)
        Pb76i <- unname(yfit$b[1])
    } else{ # no anchor
        yfit <- york(yd)
        Pb76i <- unname(yfit$b[1])
        tt <- concordiaIntersection(yfit=yfit,d=x$d)
    }
    c(tt,Pb76i)
}

concordiaIntersection <- function(yfit,d=diseq()){
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
    a <- unname(yfit$a[1])
    b <- unname(yfit$b[1])
    midpoint <- uniroot(misfit,lower=0,upper=4600,
                        a=a,b=b,d=d,gradient=TRUE)$root
    if (misfit(tt=midpoint,a=a,b=b,d=d) > 0){
        tt <- uniroot(misfit,lower=0,upper=midpoint,a=a,b=b,d=d)$root
    } else {
        tt <- optimise(misfit,lower=0,upper=4600,a=a,b=b,d=d)$minimum
    }
    tt
}
