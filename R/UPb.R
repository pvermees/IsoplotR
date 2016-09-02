get.selection.UPb <- function(x,wetherill=TRUE,...){
    if (wetherill) selection <- c('Pb207U235','Pb206U238')
    else selection <- c('U238Pb206','Pb207Pb206')
    selection
}

get.covmat.UPb <- function(x,i,...){
    out <- matrix(rep(0,16),nrow=4)
    rownames(out) <- c('Pb207U235','Pb206U238','U238Pb206','Pb207Pb206')
    colnames(out) <- rownames(out)
    if (x$format == 1){
        out['Pb207U235','Pb207U235'] <- x$x[i,'errPb207U235']^2
        out['Pb206U238','Pb206U238'] <- x$x[i,'errPb206U238']^2
        out['U238Pb206','U238Pb206'] <- x$x[i,'errU238Pb206']^2
        out['Pb207Pb206','Pb207Pb206'] <- x$x[i,'errPb207Pb206']^2

        Pb207U235 <- x$x[i,'Pb207U235']
        errPb207U235 <- x$x[i,'errPb207U235']
        Pb206U238 <- x$x[i,'Pb206U238']
        errPb206U238 <- x$x[i,'errPb206U238']
        U238Pb206 <- x$x[i,'U238Pb206']
        errU238Pb206 <- x$x[i,'errU238Pb206']
        Pb207Pb206 <- x$x[i,'Pb207Pb206']
        errPb207Pb206 <- x$x[i,'errPb207Pb206']

        R <- iratio('U238U235')[1]
        Pb206U235 <- Pb206U238*R
        errPb206U235 <- errPb206U238*R
        U235Pb206 <- U238Pb206/R
        errU235Pb206 <- errU238Pb206/R
        Pb207U238 <- Pb207U235/R
        errPb207U238 <- errPb207U235/R

        out.Pb207U235.Pb206U235 <- get.cov.xzyz(
            Pb207U235,errPb207U235,Pb206U235,errPb206U235,errPb207Pb206
        )
        out.Pb207U235.U235Pb206 <- get.cov.xzzy(
            Pb207U235,errPb207U235,U235Pb206,errU235Pb206,errPb207Pb206
        )

        out['Pb207U235','Pb206U238'] <- out.Pb207U235.Pb206U235/R
        out['Pb207U235','U238Pb206'] <- out.Pb207U235.U235Pb206*R
        out['Pb207U235','Pb207Pb206'] <- get.cov.zxzy(
            Pb207U235,errPb207U235,Pb207Pb206,errPb207Pb206,errU235Pb206
        )
        out['Pb206U238','U238Pb206'] <- -0
        out['Pb206U238','Pb207Pb206'] <- get.cov.xzzy(
            Pb207Pb206,errPb207Pb206,Pb206U238,errPb206U238,errPb207U238
        )
        out['U238Pb206','Pb207Pb206'] <- get.cov.xzyz(
            Pb207Pb206,errPb207Pb206,U238Pb206,errU238Pb206,errPb207U238
        )
        out['Pb206U238','Pb207U235'] <- out['Pb207U235','Pb206U238']
        out['U238Pb206','Pb207U235'] <- out['Pb207U235','U238Pb206']
        out['Pb207Pb206','Pb207U235'] <- out['Pb207U235','Pb207Pb206']
        out['U238Pb206','Pb206U238'] <- out['Pb206U238','U238Pb206']
        out['Pb207Pb206','Pb206U238'] <- out['Pb206U238','Pb207Pb206']
        out['Pb207Pb206','U238Pb206'] <- out['U238Pb206','Pb207Pb206']
    }
    out
}

# as.UPb: returns an 'UPb' data object
get.ratios.UPb <- function(tt,st=0,exterr=TRUE,as.UPb=FALSE){
    if (tt == 0){ tt <- 1e-10 }
    out <- list()
    l8 <- lambda('U238')[1]
    l5 <- lambda('U235')[1]
    R.x <- 1/iratio('U238U235')[1]
    R.e <- R.x*iratio('U238U235')[2]/iratio('U238U235')[1]

    Pb207U235 <- (exp(l5*tt)-1)
    Pb206U238 <- (exp(l8*tt)-1)
    U238Pb206 <- 1/(exp(l8*tt)-1)
    Pb207Pb206 <- R.x*Pb207U235/Pb206U238
    
    out$x <- c(Pb207U235,Pb206U238,U238Pb206,Pb207Pb206)

    E <- matrix(0,4,4)
    E[1,1] <- R.e^2
    E[2,2] <- lambda('U235')[2]^2
    E[3,3] <- lambda('U238')[2]^2
    E[4,4] <- st^2

    J <- matrix(0,nrow=4,ncol=4)
    J[1,2] <- tt*exp(l5*tt)
    J[1,4] <- l5*exp(l5*tt)
    J[2,3] <- tt*exp(l8*tt)
    J[2,4] <- l8*exp(l8*tt)
    J[3,3] <- -tt*exp(l8*tt)/(exp(l8*tt)-1)^2
    J[3,4] <- -l8*exp(l8*tt)/(exp(l8*tt)-1)^2
    J[4,1] <- (exp(l5*tt)-1)/(exp(l8*tt)-1)
    J[4,2] <- R.x*tt*exp(l5*tt)/(exp(l8*tt)-1)
    J[4,3] <- -R.x*(exp(l5*tt)-1)*tt*exp(l8*tt)/(exp(l8*tt)-1)^2
    J[4,4] <- R.x*(l5*exp(l5*tt)*(exp(l8*tt)-1) -
                   l8*exp(l8*tt)*(exp(l5*tt)-1))/(exp(l8*tt)-1)^2
    
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- c('Pb207U235','Pb206U238','U238Pb206','Pb207Pb206')
    rownames(out$cov) <- names(out$x)
    colnames(out$cov) <- names(out$x)
    if (as.UPb){
        result <- list()
        class(result) <- "UPb"
        result$format <- 1
        result$x <- matrix(0,1,8)
        colnames(result$x) <- c('Pb207U235','errPb207U235',
                                'Pb206U238','errPb206U238',
                                'U238Pb206','errU238Pb206',
                                'Pb207Pb206','errPb207Pb206')
        
        result$x[1,c(1,3,5,7)] <- out$x
        result$x[1,c(2,4,6,8)] <- sqrt(diag(out$cov))
        return(result)
    } else {
        return(out)
    }
}

get.Pb207U235age <- function(r75,sr75=0,exterr=TRUE){
    l5 <- lambda('U235')[1]
    sl5 <- lambda('U235')[2]
    t.75 <- log(1+r75)/l5
    J <- matrix(0,1,2)
    J[1,1] <- 1/(l5*(1+r75))
    if (exterr) J[1,2] <- log(1+r75)/l5^2
    E <- matrix(0,2,2)
    E[1,1] <- sr75^2
    E[2,2] <- sl5^2
    st.75 <- sqrt(J %*% E %*% t(J))
    out <- c(t.75,st.75)
    out
}

get.Pb206U238age <- function(r68,sr68=0,exterr=TRUE){
    l8 <- lambda('U238')[1]
    sl8 <- lambda('U238')[2]
    t.68 <- log(1+r68)/l8
    J <- matrix(0,1,2)
    J[1,1] <- 1/(l8*(1+r68))
    if (exterr) J[1,2] <- log(1+r68)/l8^2
    E <- matrix(0,2,2)
    E[1,1] <- sr68^2
    E[2,2] <- sl8^2
    st.68 <- sqrt(J %*% E %*% t(J))
    out <- c(t.68,st.68)
    out
}

#' @importFrom stats optimize
get.Pb207Pb206age <- function(r76,sr76=0,exterr=TRUE){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    sl5 <- lambda('U235')[2]
    sl8 <- lambda('U238')[2]
    R <- iratio('U238U235')[1]
    sR <- iratio('U238U235')[2]
    Pb207Pb206.misfit <- function(tt,x) {
        (get.ratios.UPb(tt)$x['Pb207Pb206']-x)^2
    }
    fit <- optimize(Pb207Pb206.misfit,c(0,4600),x=r76)
    t.76 <- fit$minimum
    J <- matrix(0,1,4)
    J[1,1] <- -(-1)/dD76dt(t.76,l5,l8,R)                      # d76/dt
    if (exterr) {
        J[1,2] <- -dD76dl5(t.76,l5,l8,R)/dD76dt(t.76,l5,l8,R) # d76/dl5
        J[1,3] <- -dD76dl8(t.76,l5,l8,R)/dD76dt(t.76,l5,l8,R) # d76/dl8
        J[1,4] <- -dD76dR(t.76,l5,l8,R)/dD76dt(t.76,l5,l8,R)  # d76/dR
    }
    E <- matrix(0,4,4)
    E[1,1] <- sr76^2
    E[2,2] <- sl5^2
    E[3,3] <- sl8^2
    E[4,4] <- sR^2
    st.76 <- sqrt( J %*% E %*% t(J) )
    out <- c(t.76,st.76)
    out
}

# x an object of class \code{UPb} returns a matrix of 7/5, 6/8, 7/6
# and concordia ages and their uncertainties.
UPb.age <- function(x,exterr=TRUE,i=NA,sigdig=NA){
    labels <- c('t.75','s[t.75]','t.68','s[t.68]',
                't.76','s[t.76]','t.conc','s[t.conc]')
    if (!is.na(i)){
        t.conc <- concordia.age(x,i=i,exterr=exterr)
        t.75 <- get.Pb207U235age(x$x[i,'Pb207U235'],x$x[i,'errPb207U235'],exterr=exterr)
        t.68 <- get.Pb206U238age(x$x[i,'Pb206U238'],x$x[i,'errPb206U238'],exterr=exterr)
        t.76 <- get.Pb207Pb206age(x$x[i,'Pb207Pb206'],x$x[i,'errPb207Pb206'],exterr=exterr)
        t.75.out <- roundit(t.75[1],t.75[2],sigdig=sigdig)
        t.68.out <- roundit(t.68[1],t.68[2],sigdig=sigdig)
        t.76.out <- roundit(t.76[1],t.76[2],sigdig=sigdig)
        t.conc.out <- roundit(t.conc$age,t.conc$age.err,sigdig=sigdig)
        out <- c(t.75.out$x,t.75.out$err,t.68.out$x,t.68.out$err,
                 t.76.out$x,t.76.out$err,t.conc.out$x,t.conc.out$err)
    } else {
        nn <- nrow(x$x)
        out <- matrix(0,nn,8)
        colnames(out) <- labels
        for (i in 1:nn){
            out[i,] <- UPb.age(x,i=i,exterr=exterr,sigdig=sigdig)
        }
    }
    out
}

dD76dt <- function(t.76,l5,l8,R){
    el8t1 <- exp(l8*t.76)-1
    el5t1 <- exp(l5*t.76)-1
    num <- el8t1*l5*exp(l5*t.76)-el5t1*l8*exp(l8*t.76)
    den <- R*el8t1^2
    num/den
}

dD76dl5 <- function(t.76,l5,l8,R){
    el8t1 <- exp(l8*t.76)-1
    t.76*exp(l5*t.76)/(R*el8t1)
}

dD76dl8 <- function(t.76,l5,l8,R){
    el8t1 <- exp(l8*t.76)-1
    t.76*exp(l5*t.76)/(R*el8t1)
}

dD76dR <- function(t.76,l5,l8,R){
    el5t1 <- exp(l5*t.76)-1
    el8t1 <- exp(l8*t.76)-1
    -el5t1/(el8t1*R^2)
}

filter.UPb.ages <- function(x,type=4,cutoff.76=1100,cutoff.disc=c(-15,5),exterr=TRUE){
    tt <- UPb.age(x,exterr=exterr)
    do.76 <- tt[,'t.68'] > cutoff.76
    if (any(is.na(cutoff.disc))) {
        is.concordant <- rep(TRUE,nrow(x))
    } else {
        disc.75.68 <- 100*(1-tt[,'t.75']/tt[,'t.68'])
        disc.68.76 <- 100*(1-tt[,'t.68']/tt[,'t.76'])
        is.concordant <- (disc.75.68>cutoff.disc[1] & disc.75.68<cutoff.disc[2]) |
                         (disc.68.76>cutoff.disc[1] & disc.68.76<cutoff.disc[2])
    }
    if (type==1){
        out <- tt[,c('t.75','s[t.75]')]
    } else if (type==2){
        out <- tt[,c('t.68','s[t.68]')]
    } else if (type==3){
        out <- tt[,c('t.76','s[t.76]')]
    } else if (type==4){
        i.76 <- as.vector(which(do.76 & is.concordant))
        i.68 <- as.vector(which(!do.76 & is.concordant))
        out <- rbind(tt[i.68,c('t.68','s[t.68]')],tt[i.76,c('t.76','s[t.76]')])
    } else if (type==4){
        out <- tt[,c('t.conc','s[t.conc]')]
    }
    out
}
