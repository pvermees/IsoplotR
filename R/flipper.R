# helper function for isochron() that flips X- and Y- axis and back
# if necessary to facilitate anchored regression and model-3 fits
flipper <- function(x,inverse=FALSE,hide=NULL,omit=NULL,
                    model=1,wtype=0,anchor=1,type,y0rat,t2rfun,...){
    xyz <- data2york(x,inverse=inverse)
    if (anchor[1]<1 & model<3){
        fitinverse <- inverse
        d2calc <- clear(xyz,hide,omit)
        fit <- regression(d2calc,model=model)
    } else if (anchor[1]==1){
        fitinverse <- FALSE
        yd <- data2york(x,inverse=fitinverse)
        d2calc <- clear(yd,hide,omit)
        anchor[2:3] <- iratio(y0rat)
        fit <- MLyork(d2calc,anchor=anchor,model=model)
    } else if (anchor[1]==2){
        fitinverse <- TRUE
        yd <- data2york(x,inverse=fitinverse)
        d2calc <- clear(yd,hide,omit)
        st <- ifelse(length(anchor)<3,0,anchor[3])
        anchor <- c(1,do.call(t2rfun,args=list(t=anchor[2],st=st,...)))
        fit <- MLyork(d2calc,anchor=anchor,model=model)
    } else if (wtype==1){
        fitinverse <- FALSE
        yd <- data2york(x,inverse=fitinverse)
        d2calc <- clear(yd,hide,omit)
        fit <- MLyork(d2calc,model=model,wtype=1)
    } else if (wtype==2){
        fitinverse <- TRUE
        yd <- data2york(x,inverse=fitinverse)
        d2calc <- clear(yd,hide,omit)
        fit <- MLyork(d2calc,model=model,wtype=1)
    } else {
        stop("Invalid anchor and/or wtype value.")
    }
    if (inverse == fitinverse) out <- fit
    else out <- invertfit(fit,type=type)
    out$xyz <- xyz
    out
}

invertfit <- function(fit,type="p"){
    out <- fit
    if (type%in%c(1,"p")){
        out$a <- quotient(X=fit$a[1],sX=fit$a[2],Y=1,sY=0,rXY=0)
        out$b <- quotient(X=fit$a[1],sX=fit$a[2],
                          Y=-fit$b[1],sY=fit$b[2],
                          rXY=fit$cov.ab/(fit$a[2]*fit$b[2]))
    } else if (type%in%c(2,"d")){
        out$a <- fit$b
        out$b <- fit$a
    } else {
        stop("Invalid isochron type.")
    }
    out
}
