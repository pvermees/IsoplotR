# currently only works for formats 1-3
mcmc.UPb <- function(x,anchor=0,seed=1,burnin=1000,mcmc=9000,...){
    set.seed(seed)
    ydat <- data2york(x,option=2)
    X <- x
    X$d <- diseq()
    init <- isochron(X,exterr=FALSE,model=1,anchor=anchor,plot=FALSE,...)
    lt <- init$logpar[1]
    la <- init$logpar[2]
    l48i <- 0
    V <- matrix(0,3,3)
    V[1:2,1:2] <- init$logcov
    V[3,3] <- 1e-5
    startvalue <- c(lt,la,l48i)
    chain <- MH(startvalue,iterations=burnin+mcmc,V=V,d=x$d,ydat=ydat)
    return(chain)
}

MH <- function(startvalue,iterations,V,...){
    chain <- array(dim = c(iterations+1,4))
    chain[1,1:3] <- startvalue
    chain[1,4] <- posterior(startvalue,...)
    colnames(chain) <- c('log[t]','log[a]','log[U48i]','posterior')
    i <- j <- 1
    while (TRUE){
        j <- j + 1
        prop <- proposal(chain[i,1:3],V=V)
        poster <- posterior(prop,...)
        prob <- exp(poster - chain[i,4])
        if (stats::runif(1) < prob){
            i <- i + 1
            chain[i,1:3] <- prop
            chain[i,4] <- poster
            if (i==iterations){
                message('acceptance rate = ',i/j)
                break
            }
            if (i %% 500 == 0){
                message('step ',i)
            }
        }
        if (j>20*iterations){
            message('acceptance rate too low (',i/j,'), so exiting early')
            break
        }
    }
    return(chain)
}

proposal <- function(param,V){
    MASS::mvrnorm(n=1,mu=param,Sigma=V)
}

prior <- function(param){
    tt <- exp(param[1])
    a <- exp(param[2])
    U48i <- exp(param[3])
    ttprior <- stats::dunif(tt, min=0, max=5000, log=TRUE)
    aprior <- stats::dunif(a, min=0, max=2, log=TRUE)
    U48iprior <- stats::dunif(U48i, min=0, max=20, log=TRUE)
    return(ttprior + aprior + U48iprior)
}

posterior <- function(param,d,ydat){
    return(likelihood(param,d=d,ydat=ydat) + prior(param))
}

likelihood <- function(param,d,ydat){
    tt <- exp(param[1])
    a <- exp(param[2])
    U48i <- exp(param[3])
    U48m <- d$U48$x
    sU48m <- d$U48$sx
    d$U48 <- list(x=U48i,sx=0,option=1)
    pred <- mclean(tt=tt,d=d)
    Pb76 <- pred$nt['Pb207']/pred$nt['Pb206']
    U8Pb6 <- pred$nt['U238']/pred$nt['Pb206']
    b <- (Pb76-a)/U8Pb6
    xy <- get.york.xy(ydat,a=a,b=b)
    dxy <- c(ydat[,'X']-xy[,1],ydat[,'Y']-xy[,2])
    nxy <- length(dxy)
    E <- matrix(0,nxy,nxy)
    diag(E) <- c(ydat[,'sX'],ydat[,'sY'])^2
    diag(E[1:(nxy/2),(nxy/2+1):nxy]) <-
        diag(E[(nxy/2+1):nxy,1:(nxy/2)]) <-
        ydat[,'rXY']*ydat[,'sX']*ydat[,'sY']
    LLUPb <- -LL.norm(dxy,E)
    LLU48 <- stats::dnorm(pred$U48,mean=U48m,sd=sU48m,log=TRUE)
    #message('U-Pb:',LLUPb,', LLU48:',LLU48)
    return(LLUPb+LLU48)
}
