# functions for flipped isochron regression, which is useful for
# INVERSE isochrons anchored by age or using model-3 with wtype=2

flippedregression <- function(x,model=1,inverse=FALSE,wtype=1,
                              anchor=0,hide=NULL,omit=NULL){
    flippedmodel3 <- (model==3 & inverse & wtype==2)
    flippedinverse <- (inverse & anchor[1]==2)
    flip <- flippedmodel3 | flippedinverse
    ANCHOR <- anchor2anchor(x,anchor=anchor,inverse=inverse)
    WTYPE <- wtypecheck(wtype=wtype,anchor=anchor)
    YD <- yd <- data2york(x,inverse=inverse)
    if (flip) YD[,c('X','sX','Y','sY','rXY')] <- yd[,c('Y','sY','X','sX','rXY')]
    d2calc <- clear(YD,hide,omit)
    out <- regression(d2calc,model=model,wtype=WTYPE,anchor=ANCHOR)
    out$inverse <- inverse
    out$flipped <- flip
    out$yd <- yd
    out
}

wtypecheck <- function(wtype,anchor){
    if (anchor[1]==1) out <- 2
    else if (anchor[1]==2) out <- 1
    else out <- wtype
    out
}

age2anchor <- function(x,...){ UseMethod("age2anchor",x) }
age2anchor.default <- function(x,anchor,inverse,...){ anchor }
age2anchor.PbPb <- function(x,tt,inverse,...){
    Pb76 <- age_to_Pb207Pb206_ratio(tt)
    if (inverse) out <- c(1,Pb76)
    else out <- c(2,Pb76)
    out
}
age2anchor.ArAr <- function(x,tt,inverse,...){
    Ar09 <- (exp(lambda("K40")[1]*tt)-1)/x$J[1]
    if (inverse) out <- c(1,1/Ar09)
    else out <- c(2,Ar09)
    out
}

# converts common ratio and age anchors to intercept and slope anchors
anchor2anchor <- function(x,anchor=0,inverse=TRUE){
    if (anchor[1]==1){
        out <- common2anchor(x,inverse=inverse)
    } else if (anchor[1]==2 && length(anchor)>1){
        out <- age2anchor(x,tt=anchor[2],inverse=inverse)
    } else {
        out <- anchor
    }
    out
}

common2anchor <- function(x,...){ UseMethod("common2anchor",x) }
common2anchor.default <- function(x,...){ stop("Not implemented") }
common2anchor.PbPb <- function(x,inverse,...){
    Pb74 <- settings('iratio','Pb207Pb204')[1]
    if (inverse) out <- c(2,Pb74)
    else out <- c(1,Pb74)
    out
}
common2anchor.ArAr <- function(x,inverse,...){
    Ar06 <- settings('iratio','Ar40Ar36')[1]
    if (inverse) out <- c(1,1/Ar06)
    else out <- c(1,Ar06)
    out
}

# converts generic regression parameters a, b and lw
# to geochronologically meaningful parameters DP, Dd and lw
# fit = the output of flippedregression
generic2DPDd <- function(fit){
    a <- fit$a['a']
    b <- fit$b['b']
    if (fit$inverse){
        J <- matrix(0,3,3)
        if (fit$flipped){
            DP <- 1/a
            Dd <- -b/a
            J[1,1] <- -1/a^2
            J[2,1] <- b/a^2
            J[2,2] <- -1/a
        } else {
            DP <- -b/a
            Dd <- 1/a
            J[1,1] <- b/a^2
            J[1,2] <- -1/a
            J[2,1] <- -1/a^2
        }
        lw <- fit$par['lw'] - 2*log(fit$par['a'])
        J[3,1] <- -2/fit$par['a']
        J[3,3] <- 1
    } else {
        DP <- fit$b[1]
        sDP <- fit$b[2]
        Dd <- fit$a[1]
        sDd <- fit$a[2]
        lw <- fit$par['lw']
        J <- diag(1,3,3)
    }
    out <- list()
    out$par <- unname(c(DP,Dd,lw))
    out$cov <- unname(J %*% fit$cov %*% t(J))
    out$DP <- c('DP'=out$par[1],'s[DP]'=sqrt(out$cov[1,1]))
    out$Dd <- c('Dd'=out$par[2],'s[Dd]'=sqrt(out$cov[2,2]))
    names(out$par) <- colnames(out$cov) <-
        rownames(out$cov) <- c('DP','Dd','lw')
    out
}
