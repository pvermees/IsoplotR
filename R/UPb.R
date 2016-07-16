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

get.ratios.UPb <- function(age){
    if (age == 0){ age <- 1e-10 }
    out <- list()
    l8 <- lambda('U238')[1]
    l5 <- lambda('U235')[1]
    R.x <- 1/iratio('U238U235')[1]
    R.e <- R.x*iratio('U238U235')[2]/iratio('U238U235')[1]

    Pb207U235 <- (exp(l5*age)-1)
    Pb206U238 <- (exp(l8*age)-1)
    U238Pb206 <- 1/(exp(l8*age)-1)
    Pb207Pb206 <- R.x*Pb207U235/Pb206U238
    
    out$x <- c(Pb207U235,Pb206U238,U238Pb206,Pb207Pb206)

    E <- matrix(0,3,3)
    E[1,1] <- R.e^2
    E[2,2] <- lambda('U235')[2]^2
    E[3,3] <- lambda('U238')[2]^2

    J <- matrix(0,nrow=4,ncol=3)
    J[1,2] <- age*exp(l5*age)
    J[2,3] <- age*exp(l8*age)
    J[3,3] <- -age*exp(l8*age)/(exp(l8*age)-1)^2
    J[4,1] <- (exp(l5*age)-1)/(exp(l8*age)-1)
    J[4,2] <- R.x*age*exp(l5*age)/(exp(l8*age)-1)
    J[4,3] <- -R.x*(exp(l5*age)-1)*age*exp(l8*age)/(exp(l8*age)-1)^2
    
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- c('Pb207U235','Pb206U238','U238Pb206','Pb207Pb206')
    rownames(out$cov) <- names(out$x)
    colnames(out$cov) <- names(out$x)
    out
}

get.Pb207U235age <- function(r75,sr75=0,dcu=TRUE){
    l5 <- lambda('U235')[1]
    sl5 <- lambda('U235')[2]
    t.75 <- log(1+r75)/l5
    J <- matrix(0,1,2)
    J[1,1] <- 1/(l5*(1+r75))
    if (dcu) J[1,2] <- log(1+r75)/l5^2
    E <- matrix(0,2,2)
    E[1,1] <- sr75^2
    E[2,2] <- sl5^2
    st.75 <- sqrt(J %*% E %*% t(J))
    out <- c(t.75,st.75)
    out
}

get.Pb206U238age <- function(r68,sr68=0,dcu=TRUE){
    l8 <- lambda('U238')[1]
    sl8 <- lambda('U238')[2]
    t.68 <- log(1+r68)/l8
    J <- matrix(0,1,2)
    J[1,1] <- 1/(l8*(1+r68))
    if (dcu) J[1,2] <- log(1+r68)/l8^2
    E <- matrix(0,2,2)
    E[1,1] <- sr68^2
    E[2,2] <- sl8^2
    st.68 <- sqrt(J %*% E %*% t(J))
    out <- c(t.68,st.68)
    out
}

#' @importFrom stats optimize
get.Pb207Pb206age <- function(r76,sr76=0,dcu=TRUE){
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
    if (dcu) {
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
UPb.age <- function(x,dcu=TRUE,i=NA,sigdig=2){
    labels <- c('t.75','s[t.75]','t.68','s[t.68]',
                't.76','s[t.76]','t.conc','s[t.conc]')
    if (!is.na(i)){
        t.conc <- concordia.age.UPb(x,i=i,dcu=dcu)
        t.75 <- get.Pb207U235age(x$x[i,'Pb207U235'],x$x[i,'errPb207U235'],dcu=dcu)
        t.68 <- get.Pb206U238age(x$x[i,'Pb206U238'],x$x[i,'errPb206U238'],dcu=dcu)
        t.76 <- get.Pb207Pb206age(x$x[i,'Pb207Pb206'],x$x[i,'errPb207Pb206'],dcu=dcu)
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
            out[i,] <- UPb.age(x,i=i,dcu=dcu)
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
