UThHe.age <- function(x){
    ns <- nrow(x)
    out <- matrix(0,ns,2)
    for (i in 1:ns){
        out[i,] <- get.UThHe.age(U=x[i,'U'],sU=x[i,'errU'],
                                 Th=x[i,'Th'],sTh=x[i,'errTh'],
                                 He=x[i,'He'],sHe=x[i,'errHe'],
                                 Sm=x[i,'Sm'],sSm=x[i,'errSm'])
    }
    out
}

get.UThHe.age <- function(U,sU,Th,sTh,He,sHe,Sm=0,sSm=0){
    R <- iratio('U238U235')[1]
    L8 <- lambda('U238')[1]
    L5 <- lambda('U235')[1]
    L2 <- lambda('Th232')[1]
    L7 <- lambda('Sm147')[1]
    f147 <- f147Sm()[1]
    P <- 8*L8*U*R/(1+R) + 7*L5*U/(1+R) + 6*L2*Th + f147*L7*Sm;
    L <- ( 8*L8*L8*U*R/(1+R) + 7*L5*L5*U/(1+R) + 
           6*Th*L2*L2 + f147*Sm*L7*L7 ) / P;
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

# atomic abundance of 147Sm (Chang et al., 2002)
f147Sm <- function(){
    c(0.1502,0.0003)    
}
