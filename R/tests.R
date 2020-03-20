# to be used in ludwig.default as:
# wtest.UPb(fit$logpar,x=x)
wtest.UPb <- function(lta0b0w,x){
    nn <- 50
    lw <- seq(from=-10,to=2,length.out=nn)
    LL <- rep(0,nn)
    for (i in 1:nn){
        lta0b0w[4] <- lw[i]
        LL[i] <- LL.lud(lta0b0w=lta0b0w,x=x,LL=TRUE)
    }
    graphics::plot(lw,LL,type='l')
}

# to be used in UThPb.R as:
# wtest(fit$par,x=x)
wtest.UThPb <- function(lta0b0wc0,x){
    XYZW <- get_XYZW(x)
    nn <- 50
    lw <- seq(from=-10,to=2,length.out=nn)
    LL <- rep(0,nn)
    for (i in 1:nn){
        lta0b0wc0[4] <- lw[i]
        LL[i] <- LL.lud.UThPb(lta0b0wc0=lta0b0wc0,x=x,
                              LL=TRUE,XYZW=XYZW)
    }
    graphics::plot(lw,LL,type='l')
}
