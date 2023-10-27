flippedregression <- function(x,model=1,inverse=FALSE,wtype=1,
                              anchor=0,hide=NULL,omit=NULL){
    abanchor <- 0
    abwtype <- wtype
    if (inverse){
        if (model==3){
            if (anchor[1]==1) { # anchor by inherited composition
                flip <- FALSE
                abanchor <- common2anchor(x,inverse=inverse)
                abwtype <- 2 # overdispersion of the slope
            } else if (anchor[1]==2) { # anchor by age
                flip <- TRUE
                abanchor <- age2anchor(x,tt=anchor[2],inverse=inverse)
                abwtype <- 2 # overdispersion of the slope
            } else { # do not anchor
                if (wtype==2){ # age dispersion
                    flip <- TRUE
                    abwtype <- 1 # overdispersion of the intercept
                } else { # dispersion of the inherited composition
                    flip <- FALSE
                }
            }
        } else { # model 1 or 2
            if (anchor[1]==2){ # anchor by age
                flip <- TRUE
                abanchor <- age2anchor(x,tt=anchor[2],inverse=inverse)
            } else if (anchor[1]==1){ # anchor by inherited composition
                flip <- FALSE
                abanchor <- common2anchor(x,inverse=inverse)
            } else { # no anchor
                flip <- FALSE
            }
        }
    } else {
        flip <- FALSE
        if (anchor[1]==2){ # anchor by age
            abanchor <- age2anchor(x,tt=anchor[2],inverse=inverse)
        } else if (anchor[1]==1){ # anchor by inherited composition
            abanchor <- common2anchor(x,inverse=inverse)
        }
    }
    YD <- yd <- data2york(x,inverse=inverse)
    if (flip) YD[,c('X','sX','Y','sY','rXY')] <- yd[,c('Y','sY','X','sX','rXY')]
    d2calc <- clear(YD,hide,omit)
    out <- regression(d2calc,model=model,abwtype=abwtype,abanchor=abanchor)
    out$inverse <- inverse
    out$flipped <- flip
    out$yd <- yd
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

flipback <- function(fit,model,wtype){
    ao <- ai <- unname(fit$a['a'])
    bo <- bi <- unname(fit$b['b'])
    lw <- unname(fit$par['lw'])
    w <- NA
    JDPDd <- matrix(0,3,3)
    Jab <- cbind(diag(1,2,2),0)
    if (fit$inverse){
        if (fit$flipped){
            DP <- 1/ai
            Dd <- -bi/ai
            ao <- -ai/bi
            bo <- 1/bi
            JDPDd[1,1] <- -1/ai^2
            JDPDd[2,1] <- bi/ai^2
            JDPDd[2,2] <- -1/ai
            Jab[1,1] <- -1/bi
            Jab[1,2] <- ai/bi^2
            Jab[2,2] <- -1/bi^2
        } else {
            DP <- -bi/ai
            Dd <- 1/ai
            JDPDd[1,1] <- bi/ai^2
            JDPDd[1,2] <- -1/ai
            JDPDd[2,1] <- -1/ai^2
        }
    } else {
        DP <- bi
        Dd <- ai
        JDPDd[1,2] <- 1
        JDPDd[2,1] <- 1
    }
    if (model==3 & ( (wtype==2 & fit$flipped) | (wtype==1 & fit$inverse))){
        w <- exp(lw)/ai^2
        JDPDd[3,1] <- -w/ai
    } else {
        w <- exp(lw)
        JDPDd[3,3] <- w
    }
    EDPDd <- unname(JDPDd %*% fit$cov %*% t(JDPDd))
    Eab <- unname(Jab %*% fit$cov %*% t(Jab))
    out <- list()
    out$yd <- fit$yd
    out$n <- fit$n
    out$model <- fit$model
    if (fit$model==1){
        out$p.value <- fit$p.value
        out$mswd <- fit$mswd
        out$df <- fit$df
    }
    out$DP <- c('DP'=DP,'s[DP]'=sqrt(EDPDd[1,1]))
    out$Dd <- c('Dd'=Dd,'s[Dd]'=sqrt(EDPDd[2,2]))
    out$disp <- c('w'=w,'s[w]'=sqrt(EDPDd[3,3]))
    out$a <- c('a'=ao,'s[a]'=sqrt(Eab[1,1]))
    out$b <- c('b'=bo,'s[b]'=sqrt(Eab[2,2]))
    out$cov.ab <- Eab[1,2]
    out
}
