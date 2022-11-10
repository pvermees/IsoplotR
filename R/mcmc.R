# currently only works for formats 1-3
mcmc.UPb <- function(x,anchor=0,seed=1,burnin=1000,mcmc=9000,...){
    ydat <- data2york(x,option=2) # formats 1-3
    yfit <- york(ydat)
    y0 <- yfit$a[1]
    
    
    init <- isochron(x,exterr=FALSE,model=1,anchor=anchor,plot=FALSE,...)
    lt <- init$logpar[1]
    la <- init$logpar[2]
    pred <- mclean(init$par[1],d=x$d)
    l48i <- log(pred$U48i)
    V <- init$logcov
    l34 <- lambda('U234')[1]
    l38 <- lambda('U238')[1]
    for (i in 1:burnin){
    }
    for (i in 1:n){
    }
}

LL.ludwig <- function(lt,la,l48i){
    
}
