# helper function for isochron() that flips X- and Y- axis and back
# if necessary to facilitate anchored regression and model-3 fits
flipper <- function(x,inverse=FALSE,hide=NULL,omit=NULL,
                    model=1,wtype=0,anchor=0,type='p',...){
    y0rat <- gety0rat(x)
    DPrat <- getDPrat(x)
    yd <- data2york(x,inverse=inverse)
    if (model<3 & anchor[1]<1){
        ifi <- rep(FALSE,3)
        d2calc <- clear(yd,hide,omit)
        fit <- regression(d2calc,model=model)
    } else if (anchor[1]==1){
        wtype <- 1 # override
        ifi <- get_ifi(wtype=wtype,type=type,inverse=inverse)
        d2calc <- flipinvert(yd=yd,ifi=ifi,type=type,hide=hide,omit=omit)
        anchor[2:3] <- iratio(y0rat)
        if (model<2) fit <- anchoredYork(d2calc,y0=anchor[2],sy0=anchor[3])
        else fit <- MLyork(d2calc,anchor=anchor,model=model)
    } else if (anchor[1]==2){
        wtype <- 2 # override
        ifi <- get_ifi(wtype=wtype,type=type,inverse=inverse)
        d2calc <- flipinvert(yd=yd,ifi=ifi,type=type,hide=hide,omit=omit)
        if (is.null(DPrat)){
            DP <- anchor[2:3]
        } else {
            st <- ifelse(length(anchor)<3,0,anchor[3])
            DP <- age2ratio(tt=anchor[2],st=st,ratio=DPrat,...)
        }
        if (model<2){
            fit <- anchoredYork(d2calc,y0=DP[1],sy0=DP[2])
        } else {
            anchor <- c(1,DP)
            fit <- MLyork(d2calc,anchor=anchor,model=model)
        }
    } else if (wtype==1){
        ifi <- get_ifi(wtype=wtype,type=type,inverse=inverse)
        d2calc <- flipinvert(yd=yd,ifi=ifi,type=type,hide=hide,omit=omit)
        fit <- MLyork(d2calc,model=model,wtype='a')
    } else if (wtype==2){
        ifi <- get_ifi(wtype=wtype,type=type,inverse=inverse)
        d2calc <- flipinvert(yd=yd,ifi=ifi,type=type,hide=hide,omit=omit)
        fit <- MLyork(d2calc,model=model,wtype='a')
    } else {
        stop("Invalid anchor and/or wtype value.")
    }
    fit$anchor <- anchor
    out <- list()
    out$flippedfit <- fit
    if (ifi[3]){
        fit <- invertfit(fit,type=type,wtype=wtype)
    }
    if (ifi[2]){
        fit <- unflipfit(fit)
    }
    if (ifi[1]){
        fit <- invertfit(fit,type=type,wtype=wtype)
    }
    out$wtype <- wtype
    out <- append(out,fit)
    out$xyz <- yd
    out
}

# ifi = invert, flip, invert
get_ifi <- function(wtype,type,inverse){
    if (wtype==1){
        if (inverse){
            out <- c(TRUE,FALSE,FALSE)
        } else {
            out <- c(FALSE,FALSE,FALSE)
        }
    } else if (wtype==2){
        if (type=='p'){
            if (inverse){
                out <- c(FALSE,TRUE,TRUE)
            } else {
                out <- c(TRUE,TRUE,TRUE)
            }
        } else if (type=='d'){
            if (inverse){
                out <- c(FALSE,FALSE,FALSE)
            } else {
                out <- c(TRUE,FALSE,FALSE)
            }
        } else {
            stop('Invalid type')
        }
    } else {
        out <- c(FALSE,FALSE,FALSE)
    }
    out
}

flipinvert <- function(yd,ifi=rep(FALSE,3),type='p',hide=NULL,omit=NULL){
    if (ifi[1]) yd <- normal2inverse(yd,type=type)
    if (ifi[2]) yd[,c('X','sX','Y','sY','rXY')] <- yd[,c(3,4,1,2,5)]
    if (ifi[3]) yd <- normal2inverse(yd)
    clear(yd,hide,omit)
}

# the purpose of flipping and inverting is to always use the intercept
# to improve the fit, therefore the inverse operations always attribute
# any overdispersion to the incoming intercept
unflipfit <- function(fit){
    out <- fit
    a <- -fit$a[1]/fit$b[1]
    b <- 1/fit$b[1]
    J11 <- -1/fit$b[1]
    J12 <- -a/fit$b[1]
    J21 <- 0
    J22 <- -b/fit$b[1]
    E11 <- fit$a[2]^2
    E22 <- fit$b[2]^2
    E12 <- fit$cov.ab
    vcovab <- errorprop(J11,J12,J21,J22,E11,E22,E12)
    out$a <- c(a,sqrt(vcovab[1]))
    out$b <- c(b,sqrt(vcovab[2]))
    out$cov.ab <- unname(vcovab[3])
    out
}
invertfit <- function(fit,type="p",wtype=0){
    out <- fit
    if (type%in%c(1,"p")){
        a <- 1/fit$a[1]
        b <- -fit$b[1]/fit$a[1]
        J11 <- -a/fit$a[1]
        J12 <- 0
        J21 <- -b/fit$a[1]
        J22 <- -1/fit$a[1]
        E11 <- fit$a[2]^2
        E22 <- fit$b[2]^2
        E12 <- fit$cov.ab
        vcovab <- errorprop(J11,J12,J21,J22,E11,E22,E12)
        out$a <- c(a,sqrt(vcovab[1]))
        out$b <- c(b,sqrt(vcovab[2]))
        out$cov.ab <- vcovab[3]
    } else if (type%in%c(2,"d")){
        out$a <- fit$b
        out$b <- fit$a
    } else {
        stop("Invalid isochron type.")
    }
    names(out$a) <- c('a','s[a]')
    names(out$b) <- c('b','s[b]')
    out
}

anchoredYork <- function(x,y0=0,sy0=0){
    eps <- .Machine$double.eps
    X <- rbind(x,c(0,eps,y0,max(sy0,eps),0))
    out <- yorkhelper(X,np=1)
    if (y0==0) out$a[1] <- 0
    if (sy0==0) out$a[2] <- 0
    out$model <- 1
    out$n <- nrow(x)
    out
}

#' Returns the dependent variable of a conventional isochron
#'
#' Helper function that returns a string with the Dd ratio of a given
#' chronometer
#' @param x an IsoplotR data object
#' @param ... optional arguments
#' @noRd
gety0rat <- function(x,...){ UseMethod("gety0rat",x) }
#' @noRd
gety0rat.default <- function(x,...){ NULL }
#' @noRd
gety0rat.ArAr <- function(x,...){ 'Ar40Ar36' }
#' @noRd
gety0rat.PbPb <- function(x,...){ 'Pb207Pb204' }
#' @noRd
gety0rat.ThPb <- function(x,...){ 'Pb208Pb204' }
#' @noRd
gety0rat.KCa <- function(x,...){ paste0('Ca40Ca',x$sister) }
#' @noRd
gety0rat.RbSr <- function(x,...){ 'Sr87Sr86' }
#' @noRd
gety0rat.ReOs <- function(x,...){ 'Os187Os188' }
#' @noRd
gety0rat.SmNd <- function(x,...){ 'Nd143Nd144' }
#' @noRd
gety0rat.LuHf <- function(x,...){ 'Hf176Hf177' }

#' Returns the Daughter-Parent ratio
#'
#' Helper function that returns a string with the DP-ratio of a given
#' chronometer
#' 
#' @param x an IsoplotR data object
#' @param ... optional arguments
#' @noRd
getDPrat <- function(x,...){ UseMethod("getDPrat",x) }
#' @noRd
getDPrat.default <- function(x,...){ NULL }
#' @noRd
getDPrat.ArAr <- function(x,...){ 'Ar40Ar39' }
#' @noRd
getDPrat.PbPb <- function(x,...){ 'Pb207Pb206' }
#' @noRd
getDPrat.ThPb <- function(x,...){ 'Pb208Th232' }
#' @noRd
getDPrat.KCa <- function(x,...){ 'Ca40K40' }
#' @noRd
getDPrat.RbSr <- function(x,...){ 'Sr87Rb87' }
#' @noRd
getDPrat.ReOs <- function(x,...){ 'Os187Re187' }
#' @noRd
getDPrat.SmNd <- function(x,...){ 'Nd143Sm147' }
#' @noRd
getDPrat.LuHf <- function(x,...){ 'Hf176Lu176' }

#' Returns the radioactive parent
#'
#' Helper function that returns a string with the parent nuclide
#' of a given chronometer
#' 
#' @param x an IsoplotR data object
#' @param ... optional arguments
#' @noRd
getParent <- function(x,...){ UseMethod("getParent",x) }
#' @noRd
getParent.default <- function(x,...){ NULL }
#' @noRd
getParent.ThPb <- function(x,...){ 'Th232' }
#' @noRd
getParent.KCa <- function(x,...){ 'K40' }
#' @noRd
getParent.RbSr <- function(x,...){ 'Rb87' }
#' @noRd
getParent.ReOs <- function(x,...){ 'Re187' }
#' @noRd
getParent.SmNd <- function(x,...){ 'Sm147' }
#' @noRd
getParent.LuHf <- function(x,...){ 'Lu176' }
