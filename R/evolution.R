#' Th-U evolution diagram
#'
#' Plots Th-U data on a
#' \eqn{^{234}}U/\eqn{^{238}}U-\eqn{^{230}}Th/\eqn{^{238}}U evolution
#' diagram or a \eqn{^{234}}U/\eqn{^{238}}U-age diagram, calculates
#' isochron ages.
#'
#' @param x an object of class \code{ThU}
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param alpha confidence cutoff for the error ellipses
#' @param transform if \code{TRUE}, plots \eqn{^{234}}U/\eqn{^{238}}U
#'     vs. Th-U age.
#' @param detrital apply a detrital Th correction by projecting the
#'     compositions along an isochron?
#' @param show.numbers label the error ellipses with the grain
#'     numbers?
#' @param ellipse.col background colour of the error ellipses
#' @param line.col colour of the age grid
#' @param isochron fit a 3D isochron to the data?
#' @param exterr propagate the decay constant uncertainty in the
#'     isochron age?
#' @param sigdig number of significant digits for the isochron age
#' @param ... optional arguments to the generic \code{plot} function
#' @examples
#' data(examples)
#' evolution(examples$ThU)
#' @references Ludwig, K.R. and Titterington, D.M., 1994. Calculation
#'     of \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#'     Cosmochimica Acta, 58(22), pp.5031-5042.
#'
#' Ludwig, K.R., 2003. Mathematical-statistical treatment of data and
#'     errors for 230 Th/U geochronology. Reviews in Mineralogy and
#'     Geochemistry, 52(1), pp.631-656.
#' @importFrom grDevices rgb
#' @export
evolution <- function(x,xlim=NA,ylim=NA,alpha=0.05,transform=FALSE,
                      detrital=FALSE,show.numbers=FALSE,
                      ellipse.col=rgb(0,1,0,0.5),
                      line.col='darksalmon',isochron=FALSE,
                      exterr=TRUE,sigdig=2,...){
    if (transform){
        U4U8vst(x,detrital=detrital,xlim=xlim,ylim=ylim,alpha=alpha,
                show.numbers=show.numbers,ellipse.col=ellipse.col,...)
    } else {
        U4U8vsTh0U8(x,isochron=isochron,detrital=detrital,xlim=xlim,
                    ylim=ylim,alpha=alpha, show.numbers=show.numbers,
                    ellipse.col=ellipse.col, line.col=line.col,...)
    }
    if (isochron){
        fit <- isochron.ThU(x,type=3,plot=FALSE,exterr=exterr)
        graphics::title(evolution.title(fit,sigdig=sigdig))
    }
}

U4U8vst <- function(x,detrital=FALSE,xlim=NA,ylim=NA,
                    alpha=0.05,show.numbers=FALSE,
                    ellipse.col=grDevices::rgb(0,1,0,0.5),...){
    ns <- length(x)
    ta0 <- ThU.age(x,exterr=FALSE,i2i=detrital,cor=FALSE)
    nsd <- 3
    if (any(is.na(xlim))) xlim <- range(c(ta0[,'t']-nsd*ta0[,'s[t]'],
                                          ta0[,'t']+nsd*ta0[,'s[t]']))
    if (any(is.na(ylim))) ylim <- range(c(ta0[,'48_0']-nsd*ta0[,'s[48_0]'],
                                          ta0[,'48_0']+nsd*ta0[,'s[48_0]']))
    x.lab <- 'Age [ka]'
    y.lab <- expression(paste("("^"234","U/"^"238","U)"[o]))
    graphics::plot(xlim,ylim,type='n',bty='n',xlab=x.lab,ylab=y.lab)
    covmat <- matrix(0,2,2)
    for (i in 1:ns){
        x0 <- ta0[i,'t']
        y0 <- ta0[i,'48_0']
        diag(covmat) <- ta0[i,c('s[t]','s[48_0]')]^2
        covmat[1,2] <- ta0[i,'cov[t,48_0]']
        covmat[2,1] <- covmat[1,2]
        ell <- ellipse(x0,y0,covmat,alpha=alpha)
        graphics::polygon(ell,col=ellipse.col)
        graphics::points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) graphics::text(x0,y0,i)
    }
}

U4U8vsTh0U8 <- function(x,isochron=FALSE,detrital=FALSE, xlim=NA,
                        ylim=NA, alpha=0.05,show.numbers=FALSE,
                        ellipse.col=grDevices::rgb(0,1,0,0.5),
                        line.col='darksalmon',...){
    ns <- length(x)
    d <- data2evolution(x,detrital=detrital)
    lim <- evolution.lines(d,xlim=xlim,ylim=ylim,...)
    if (isochron){
        fit <- isochron(x,type=2,plot=FALSE)
        b48 <- fit$par['a']
        b08 <- fit$par['A']
        e48 <- 1
        e08 <- fit$par['A'] + fit$par['B']*(e48-fit$par['a'])/fit$par['b']
        graphics::lines(c(b08,e08),c(b48,e48))
    }
    covmat <- matrix(0,2,2)
    for (i in 1:ns){
        x0 <- d[i,'Th230U238']
        y0 <- d[i,'U234U238']
        diag(covmat) <- d[i,c('errTh230U238','errU234U238')]^2
        covmat[1,2] <- d[i,'cov']
        covmat[2,1] <- covmat[1,2]
        ell <- ellipse(x0,y0,covmat,alpha=alpha)
        graphics::polygon(ell,col=ellipse.col)
        graphics::points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) { graphics::text(x0,y0,i) }
    }
    if (isochron){
        sa <- sqrt(fit$cov['a','a'])
        sA <- sqrt(fit$cov['A','A'])
        ell <- matrix(c(fit$par['A'],sA,fit$par['a'],sa,
                        fit$cov['a','A']/(sa*sA)),1,5)
        scatterplot(ell,alpha=alpha,
                    ellipse.col=grDevices::rgb(1,1,1,0.85),
                    line.col='black',new.plot=FALSE)
    }
}

evolution.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age[1],fit$age[2],sigdig=sigdig)
    rounded.a0 <- roundit(fit$y0[1],fit$y0[2],sigdig=sigdig)
    line1 <- substitute('isochron age ='~a%+-%b~'[ka] (1'~sigma~')',
                        list(a=rounded.age[1], b=rounded.age[2]))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,2),
                             b=signif(fit$p.value,2)))
    line3 <- substitute('(234/238)'[o]*' ='~a%+-%b~'(1'~sigma~')',
                        list(a=rounded.a0[1], b=rounded.a0[2]))
    graphics::mtext(line1,line=2)
    graphics::mtext(line2,line=1)
    graphics::mtext(line3,line=0)
}

evolution.lines <- function(d,xlim=NA,ylim=NA,bty='n',
                            line.col='darksalmon',...){
    nn <- 20
    maxt <- 400
    tt <- seq(from=0,to=maxt,by=50)
    nsd <- 3
    if (any(is.na(xlim))) max.dx <- max(d[,'Th230U238']+nsd*d[,'errTh230U238'])
    else max.dx <- xlim[2] # only used if ylim == NA
    if (any(is.na(ylim))){
        max.dy <- max(d[,'U234U238']+nsd*d[,'errU234U238'])
        a0max <- get.ThU.age(max.dy,0,max.dx,0,0,exterr=FALSE)['48_0']
    } else {
        a0max <- ylim[2]
    }    
    if (any(is.na(xlim)))
        xlim <- range(pretty(c(0,max(get.Th230U238(tt,a0max)))))
    a0 <- pretty(c(1,a0max))
    if (any(is.na(ylim))){
        ylim <- range(a0)
        a0 <- a0[a0>1]
    }
    else {
        a0 <- c(a0[a0>1 & a0<ylim[2]],ylim[2])
    }
    x.lab <- expression(paste(""^"230","Th/"^"238","U"))
    y.lab <- expression(paste(""^"234","U/"^"238","U"))
    graphics::plot(xlim,ylim,type='n',xlab=x.lab,ylab=y.lab,bty=bty,...)
    na0 <- length(a0)
    for (i in 1:na0){
        ttt <- seq(0,maxt,length.out=nn)
        x <- get.Th230U238(ttt,a0[i])
        y <- get.U234U238(ttt,a0[i])
        graphics::lines(x,y,col=line.col)
    }
    for (i in 2:nn){
        x <- get.Th230U238(tt[i],a0)
        y <- get.U234U238(tt[i],a0)
        x0 <- get.Th230U238(tt[i],1)
        graphics::lines(c(x0,x),c(1,y),col=line.col)
        graphics::text(x[na0],y[na0],tt[i],pos=4,col=line.col)
    }
    rbind(xlim,ylim)
}

data2evolution <- function(x,detrital=FALSE){
    out <- list()
    labels <- c('Th230U238','errTh230U238',
                'U234U238','errU234U238','cov')
    if (x$format == 1){
        ns <- length(x)
        out <- matrix(0,ns,5)
        colnames(out) <- labels
        out[,'Th230U238'] <- x$x[,'Th230Th232']/x$x[,'U238Th232']
        out[,'U234U238'] <- x$x[,'U234Th232']/x$x[,'U238Th232']
        J <- matrix(0,2,3)
        for (i in 1:ns){
            J[1,1] <- -out[i,'Th230U238']/x$x[i,'U238Th232'] # d/d82
            J[1,2] <- 0                                      # d/d42
            J[1,3] <- 1/x$x[i,'U238Th232']                   # d/d02
            J[2,1] <- -out[i,'U234U238']/x$x[i,'U238Th232']
            J[2,2] <- 1/x$x[i,'U238Th232']
            J[2,3] <- 0
            E <- cor2cov3(x$x[i,'errU238Th232'],x$x[i,'errU234Th232'],
                          x$x[i,'errTh230Th232'],x$x[i,'rhoXY'],
                          x$x[i,'rhoXZ'],x$x[i,'rhoYZ'])
            covmat <- J %*% E %*% t(J)
            out[i,'errTh230U238'] <- sqrt(covmat[1,1])
            out[i,'errU234U238'] <- sqrt(covmat[2,2])
            out[i,'cov'] <- covmat[1,2]
        }
    } else {
        xy <- x$x[,c('Th230U238','errTh230U238',
                      'U234U238','errU234U238')]
        covariance <- x$x[,'errTh230U238']*
            x$x[,'errU234U238']*x$x[,'rhoYZ']
        out <- cbind(xy,covariance)
        colnames(out) <- labels
    }
    if (detrital){
        osmond <- data2tit.ThU(x,osmond=TRUE)
        fit <- titterington(osmond)
        out[,'U234U238'] <- out[,'U234U238'] - fit$par['b']*osmond[,'X']
        out[,'Th230U238'] <- out[,'Th230U238'] - fit$par['B']*osmond[,'X']
    }
    out
}
