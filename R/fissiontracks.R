fissiontrack.age <- function(x,i=NA,sigdig=2,external=TRUE){
    ns <- nrow(x$x)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    if (external){
        zeta <- x$zeta
        rhoD <- x$rhoD
    } else {
        zeta <- c(x$zeta[1],0)
        rhoD <- c(x$rhoD[1],0)
    }
    for (j in 1:ns){
        tt <- get.fissiontrack.age(x$x[j,'Ns'],x$x[j,'Ni'],
                                   zeta,rhoD)
        t.out <- roundit(tt[1],tt[2],sigdig=sigdig)
        out[j,] <- c(t.out$x,t.out$err)
    }
    if (!is.na(i)) out <- out[i,]
    out
}

# zeta and rhoD are two-element vectors
get.fissiontrack.age <- function(Ns,Ni,zeta,rhoD){
    L8 <- lambda('U238')[1]
    out <- c(0,0)
    out[1] <- log(1+0.5*L8*(zeta[1]/1e6)*rhoD[1]*(Ns/Ni))/L8
    out[2] <- out[1]*sqrt(1/Ns + 1/Ni +
                          (rhoD[2]/rhoD[1])^2 +
                          (zeta[2]/zeta[1])^2)
    out
}
