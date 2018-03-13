UThHe.age <- function(x,i=NA,sigdig=NA){
    ns <- nrow(x)
    doSm <- doSm(x)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    for (j in 1:ns){
        if (doSm) tt <- get.UThHe.age(U=x[j,'U'],sU=x[j,'errU'],
                                      Th=x[j,'Th'],sTh=x[j,'errTh'],
                                      He=x[j,'He'],sHe=x[j,'errHe'],
                                      Sm=x[j,'Sm'],sSm=x[j,'errSm'])
        else tt <- get.UThHe.age(U=x[j,'U'],sU=x[j,'errU'],
                                 Th=x[j,'Th'],sTh=x[j,'errTh'],
                                 He=x[j,'He'],sHe=x[j,'errHe'])
        out[j,] <- roundit(tt[1],tt[2],sigdig=sigdig)
    }
    if (!is.na(i)) out <- out[i,]
    out
}

get.UThHe.age <- function(U,sU,Th,sTh,He,sHe,Sm=0,sSm=0){
    R <- iratio('U238U235')[1]
    L8 <- lambda('U238')[1]
    L5 <- lambda('U235')[1]
    L2 <- lambda('Th232')[1]
    L7 <- lambda('Sm147')[1]
    f147 <- f147Sm()[1]
    P <- 8*L8*U*R/(1+R) + 7*L5*U/(1+R) + 6*L2*Th + f147*L7*Sm
    L <- ( 8*L8*L8*U*R/(1+R) + 7*L5*L5*U/(1+R) + 
           6*Th*L2*L2 + f147*Sm*L7*L7 ) / P
    tt <- log(1 + L*He/P)/L;
    D <- 8*(exp(L8*tt)-1)*U*R/(1+R) + 7*(exp(L5*tt)-1)*U/(1+R) +
         6*(exp(L2*tt)-1)*Th + f147*(exp(L7*tt)-1)*Sm - He
    dDdt <- 8*L8*exp(L8*tt)*U*R/(1+R) + 7*L5*exp(L5*tt)*U/(1+R) +
            6*L2*exp(L2*tt)*Th + f147*L7*exp(L7*tt)*Sm
    dDdU <- 8*(exp(L8*tt)-1)*R/(1+R) + 7*(exp(L5*tt)-1)/(1+R)
    dDdTh <- 6*L2*exp(L2*tt)
    dDdSm <- f147*L7*exp(L7*tt)
    dDdHe <- -1
    J <- matrix(0,1,4)
    E <- matrix(0,4,4)
    J[1,1] <- -dDdU/dDdt
    J[1,2] <- -dDdTh/dDdt
    J[1,3] <- -dDdSm/dDdt
    J[1,4] <- -dDdHe/dDdt
    E[1,1] <- sU^2
    E[2,2] <- sTh^2
    E[3,3] <- sSm^2
    E[4,4] <- sHe^2
    st <- sqrt( J %*% E %*% t(J) )
    c(tt,st)
}

flat.uv.table <- function(x,w=0){
    ns <- length(x)
    out <- matrix(0,ns,5)
    for (i in 1:ns){
        uvc <- UThHe2uv.covmat(x,i,w=w)
        out[i,c(1,3)] <- uvc$uv
        out[i,2] <- sqrt(uvc$covmat[1,1])
        out[i,4] <- sqrt(uvc$covmat[2,2])
        out[i,5] <- uvc$covmat[1,2]/(out[i,2]*out[i,4])
    }
    colnames(out) <- c('u','s[u]','v','s[v]','cor[u,v]')
    out
}
flat.uvw.table <- function(x,w=0){
    ns <- length(x)
    out <- matrix(0,ns,9)
    for (i in 1:ns){
        uvwc <- UThHe2uvw.covmat(x,i,w=w)
        out[i,c(1,3,5)] <- uvwc$uvw
        out[i,2] <- sqrt(uvwc$covmat[1,1])
        out[i,4] <- sqrt(uvwc$covmat[2,2])
        out[i,6] <- sqrt(uvwc$covmat[3,3])
        out[i,7] <- uvwc$covmat[1,2]/(out[i,2]*out[i,4])
        out[i,8] <- uvwc$covmat[1,3]/(out[i,2]*out[i,6])
        out[i,9] <- uvwc$covmat[2,3]/(out[i,4]*out[i,6])
    }
    colnames(out) <- c('u','s[u]','v','s[v]','w','s[w]',
                       'rho[u,v]','rho[u,w]','rho[v,w]')
    out
}

get.He <- function(tt,U,Th,Sm=0){
    R <- iratio('U238U235')[1]
    L8 <- lambda('U238')[1]
    L5 <- lambda('U235')[1]
    L2 <- lambda('Th232')[1]
    L7 <- lambda('Sm147')[1]
    f147 <- f147Sm()[1]
    aa <- 8*R*(exp(L8*tt)-1)/(1+R) +
          7*(exp(L5*tt)-1)/(1+R)
    bb <- 6*(exp(L2*tt)-1)
    cc <- f147*(exp(L7*tt)-1)
    aa*U + bb*Th + cc*Sm
}

# atomic abundance of 147Sm (Chang et al., 2002)
f147Sm <- function(){
    c(0.1502,0.0003)    
}

doSm <- function(x){
    (ncol(x) == 8) && all(is.finite(x[,c(7,8)])) && all(x[,c(7,8)]>0)
}
