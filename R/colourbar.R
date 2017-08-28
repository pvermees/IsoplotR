colourbar <- function(z=c(0,1),col=c("#00FF0080","#FF000080"),corner=2,
                      buffer=0,strip.width=0.02,strip.length=1,...){
    ucoord <- par()$usr
    plotwidth <- (ucoord[2]-ucoord[1])
    plotheight <- (ucoord[4]-ucoord[3])
    if (corner==1){
        xb <- ucoord[1] + buffer*plotwidth
        xe <- xb + strip.length*plotwidth
        yb <- ucoord[3] + buffer*plotheight
        ye <- yb + strip.width*plotheight
    } else if (corner==2){
        xe <- ucoord[2] - buffer*plotwidth
        xb <- xe - strip.width*plotwidth
        yb <- ucoord[3] + buffer*plotwidth
        ye <- yb + strip.length*plotheight
    }
    ndiv <- 50 # number of divisions
    dx <- (xe-xb)/ndiv
    dy <- (ye-yb)/ndiv
    zz <- seq(from=min(z),to=max(z),length.out=ndiv)
    cc <- levels2colours(levels=zz,colours=col)
    for (i in 1:ndiv){
        if (corner==1)
            rect(xb+(i-1)*dx,yb,xb+i*dx,ye,col=cc[i],border=NA)
        else if (corner==2)
            rect(xb,yb+(i-1)*dy,xe,yb+i*dy,col=cc[i],border=NA)
    }
    rect(xb,yb,xe,ye)
}
