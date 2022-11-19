# Ludwig regression with initial disequilibrium
ludwigd <- function(x,bayes=FALSE){
    if (x$format<4){
        init <- init.ludwigd(x)
        fit <- optifix(parms=init$pars,fn=LL.ludwigd,fixed=c(FALSE,FALSE,FALSE),
                       method="L-BFGS-B",lower=init$lower,upper=init$upper,
                       x=x,hessian=TRUE)
    } else if (x$format<7){
        # TODO
    } else {
        # TODO
    }
    fit
}

init.ludwigd <- function(x){
    out <- list()
    if (x$format<4){
        yd <- data2york(x,option=2)
        yfit <- york(yd)
        PbU0 <- -yfit$b[1]/yfit$a[1]
        lt <- log(get.Pb206U238.age(x=PbU0)[1])
        la0 <- log(yfit$a[1])
        out$pars <- c(lt,la0,0) # log[t], log[a0], log[U48i]
        out$lower <- c(out$pars[1]-2,out$pars[2]-1,log(0.01))
        out$upper <- c(out$pars[1]+2,out$pars[2]+1,log(20))
    }
    out
}

LL.ludwigd <- function(pars,x,exterr=FALSE){
    lta0b0w <- pars[1:2]
    X <- x
    U48i <- exp(pars[3])
    X$d$U48$x <- U48i
    X$d$U48$option <- 1
    LL1 <- data2ludwig(X,lta0b0w,exterr=exterr)$LL
    pred <- mclean(tt=exp(pars[1]),d=X$d)
    U48m <- x$d$U48$x
    sU48m <- x$d$U48$sx
    LL2 <- dnorm(x=pred$U48,mean=U48m,sd=sU48m,log=TRUE)
    -(LL1+LL2)
}
