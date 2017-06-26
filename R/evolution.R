evolution <- function(x,xlim=c(0,2.4),ylim=c(1,3),alpha=0.05,
                      show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),
                      lwd=1,line.col='darksalmon',exterr=TRUE,
                      show.age=1,sigdig=2,...){
    fit <- isochron.ThU(x,type=3,plot=FALSE,exterr=exterr)
    d <- data2evolution(x,fit)
    x.lab <- expression(paste(""^"232","Th/"^"238","U"))
    y.lab <- expression(paste(""^"230","Th/"^"238","U"))
    evolution.lines(d,xlim=xlim,ylim=ylim,xlab=x.lab,ylab=y.lab,...)
    covmat <- matrix(0,2,2)
    for (i in 1:nrow(d)){
        x0 <- d[i,'Th230U238']
        y0 <- d[i,'U234U238']
        diag(covmat) <- d[i,c('errTh230U238','errU234U238')]^2
        covmat[1,2] <- d[i,'rho']*d[i,'errTh230U238']*d[i,'errU234U238']
        covmat[2,1] <- covmat[1,2]
        ell <- ellipse(x0,y0,covmat,alpha=alpha)
        polygon(ell,col=ellipse.col)
        points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) { text(x0,y0,i) }
    }
    title(evolution.title(fit,sigdig=sigdig),xlab=x.lab,ylab=y.lab)
}

evolution.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age[1],fit$age[2],sigdig=sigdig)
    rounded.a0 <- roundit(fit$age[3],fit$age[4],sigdig=sigdig)
    line1 <- substitute('isochron age ='~a%+-%b~'[Ma] (1'~sigma~')',
                        list(a=rounded.age[1], b=rounded.age[2]))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,2),
                             b=signif(fit$p.value,2)))
    line3 <- substitute('234/238_0 ='~a%+-%b~'(1'~sigma~')',
                        list(a=rounded.a0[1], b=rounded.a0[2]))
    graphics::mtext(line1,line=2)
    graphics::mtext(line2,line=1)
    graphics::mtext(line3,line=0)
}

evolution.lines <- function(d,xlim=c(0,2.4),ylim=c(1,3),
                            xlab='',ylab='',bty='n',...){
    plot(xlim,ylim,type='n',xlab=xlab,ylab=ylab,bty=bty,...)
    nn <- 100
    maxt <- 400
    tt <- seq(from=0,to=maxt,by=50)
    a0 <- pretty(ylim)
    a0 <- a0[a0>1]
    for (i in 1:length(a0)){
        ttt <- seq(0,maxt,length.out=nn)
        x <- get.Th230U238(ttt,a0[i])
        y <- get.U234U238(ttt,a0[i])
        lines(x,y)
    }
    for (i in 1:length(tt)){
        x <- get.Th230U238(tt[i],a0)
        y <- get.U234U238(tt[i],a0)
        x0 <- get.Th230U238(tt[i],1)
        lines(c(x0,x),c(1,y))
    }
}

data2evolution <- function(x,fit){
    labels <- c('Th230U238','errTh230U238',
                'U234U238','errU234U238','rho')

    if (x$format == 1){
        ns <- length(x)
        out <- matrix(0,ns,5)
        out[,'Th230U238'] <- x$x[,'Th230Th232']/x$x[,'U238Th232']
        out[,'U234U238'] <- x$x[,'U234Th232']/x$x[,'U238Th232']
        J <- matrix(0,2,3)
        for (i in 1:ns){
            J[1,1] <- -out[,'Th230U238']/x$x[,'U238Th232'] # d/d82
            J[1,2] <- 0                                    # d/d42
            J[1,3] <- 1/x$x[,'U238Th232']                  # d/d02
            J[2,1] <- -out[,'U234U238']/x$x[,'U238Th232']
            J[2,2] <- 1/x$x[,'U238Th232']
            J[2,3] <- 0
            E <- cor2cov3(x$x[,'errU238Th232'],x$x[,'errU234Th232'],
                          x$x[,'errTh230Th232'],x$x[,'rhoXY'],
                          x$x[,'rhoXZ'],x$x[,'rhoYZ'])
            covmat <- J %*% E %*% t(J)
            out[i,'errTh230U238'] <- sqrt(covmat[1,1])
            out[i,'errU234U238'] <- sqrt(covmat[2,2])
            out[i,'rho'] <- covmat[1,2]
        }
    } else {
        out <- x$x[,c('Th230U238','errTh230U238',
                      'U234U238','errU234U238','rhoYZ')]
        colnames(out) <- labels
    }
    out
}
