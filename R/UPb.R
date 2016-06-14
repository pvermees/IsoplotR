get.covmat.UPb <- function(x,i){
    covmat <- matrix(rep(0,16),nrow=4)
    rownames(covmat) <- c('Pb207U235','Pb206U238','U238Pb206','Pb207Pb206')
    colnames(covmat) <- rownames(covmat)
    if (x$format == 1){
        covmat[1,1] <- x$x[i,'errPb207U235']^2
        covmat[2,2] <- x$x[i,'errPb206U238']^2
        covmat[3,3] <- (x$x[i,'errPb206U238']^2)/(x$x[i,'Pb206U238']^4)
        covmat[4,4] <- x$x[i,'errPb207Pb206']^2

        covmat[1,2] <- 0.5 * x$x[i,'Pb207U235'] * x$x[i,'Pb206U238'] * (
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 +
                            (x$x[i,'errPb206U238']/x$x[i,'Pb206U238'])^2 -
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 )
        covmat[1,3] <- 0.5 * x$x[i,'Pb207U235'] * x$x[i,'Pb206U238'] * (
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 -
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 -
                            covmat[3,3]*x$x[i,'Pb206U238']^2 )
        covmat[1,4] <- 0.5 * x$x[i,'Pb207U235']*x$x[i,'Pb207Pb206'] * (
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 +
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 -
                            (x$x[i,'errPb206U238']/x$x[i,'Pb206U238'])^2 )
        covmat[2,1] <- covmat[1,2]
        covmat[2,3] <- -1
        covmat[2,4] <- 0.5 * x$x[i,'Pb207Pb206']*x$x[i,'Pb206U238'] * (
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 -
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 -
                            (x$x[i,'errPb206U238']/x$x[i,'Pb206U238'])^2 )
        covmat[3,1] <- covmat[1,3]
        covmat[3,2] <- -1
        covmat[3,4] <- 0.5 * x$x[i,'Pb207Pb206']/x$x[i,'Pb206U238'] * (
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 +
                             covmat[3,3]*x$x[i,'Pb206U238']^2 -
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 )
        covmat[4,1] <- covmat[1,4]
        covmat[4,2] <- covmat[2,4]
        covmat[4,3] <- covmat[3,4]
    }
    covmat
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

get.Pb207U235age <- function(Pb207U235){
    log(1+Pb207U235)/lambda('U235')[1]
}

get.Pb206U238age <- function(Pb206U238){
    log(1+Pb206U238)/lambda('U238')[1]
}

#' @importFrom stats optimize
get.Pb207Pb206age <- function(Pb207Pb206){
    Pb207Pb206.misfit <- function(x,y) { (get.ratios.UPb(x)$x['Pb207Pb206'] - y)^2 }
    out <- optimize(Pb207Pb206.misfit,c(0,4600),y=Pb207Pb206)
    out$minimum
}

# x an object of class \code{UPb}
# returns a matrix of 7/5, 6/8, 7/6 and concordia ages
# and their uncertainties. If i!=NA, returns a list with those
# same ages of aliquot i and their covariance matrix
UPb.age <- function(x,i=NA){
    labels <- c('t.75','s[t.75]','t.68','s[t.68]',
                't.76','s[t.76]','t.conc','s[t.conc]')
    if (!is.na(i)){
        tc <- concordia.age.UPb(x,i=i)
        t.75 <- get.Pb207U235age(x$x[i,'Pb207U235'])
        t.68 <- get.Pb206U238age(x$x[i,'Pb206U238'])
        t.76 <- get.Pb207Pb206age(x$x[i,'Pb207Pb206'])
        t.conc <- tc$age
        E <- E.UPb.age(x,i)
        J <- J.UPb.age(x,i,t.76)
        covmat <- J %*% E %*% t(J)
        ages <- c(t.75,t.68,t.76,t.conc)
        stderrs <- c(sqrt(diag(covmat)),tc$age.err)
        t.75.out <- roundit(ages[1],stderrs[1])
        t.68.out <- roundit(ages[2],stderrs[2])
        t.76.out <- roundit(ages[3],stderrs[3])
        t.conc.out <- roundit(ages[4],stderrs[4])
        out <- c(t.75.out$x,t.75.out$err,t.68.out$x,t.68.out$err,
                 t.76.out$x,t.76.out$err,t.conc.out$x,t.conc.out$err)
    } else {
        nn <- nrow(x$x)
        out <- matrix(0,nn,8)
        colnames(out) <- labels
        for (i in 1:nn){
            out[i,] <- UPb.age(x,i)
        }
    }
    out
}

E.UPb.age <- function(x,i){
    out <- matrix(0,6,6)
    out[1,1] <- x$x[i,'errPb207U235']^2
    out[2,2] <- x$x[i,'errPb206U238']^2
    out[3,3] <- x$x[i,'errPb207Pb206']^2
    out[4,4] <- lambda('U235')[2]^2
    out[5,5] <- lambda('U238')[2]^2
    out[6,6] <- iratio('U238U235')[2]^2
    out
}

J.UPb.age <- function(x,i,t.76){
    r75 <- x$x[i,'Pb207U235']
    r68 <- x$x[i,'Pb206U238']
    r76 <- x$x[i,'Pb207Pb206']
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- iratio('U238U235')[1]
    out <- matrix(0,3,6)
    out[1,1] <- 1/(l5*(1+r75)) # d(t75)/d(75)
    out[1,4] <- -log(1+r75)/l5^2 # d(t75)/d(l5)
    out[2,2] <- 1/(l8*(1+r68))  # d(t68)/d(68)
    out[2,5] <- -log(1+r68)/l8^2 # d(t68)/d(l8)
    out[3,3] <- -(-1)/dD76dt(t.76,l5,l8,R)
    out[3,4] <- -dD76dl5(t.76,l5,l8,R)/dD76dt(t.76,l5,l8,R)
    out[3,5] <- -dD76dl8(t.76,l5,l8,R)/dD76dt(t.76,l5,l8,R)
    out[3,6] <- -dD76dR(t.76,l5,l8,R)/dD76dt(t.76,l5,l8,R)
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
