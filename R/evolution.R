#' Th-U evolution diagram
#'
#' Plots Th-U data on a
#' \eqn{^{234}}U/\eqn{^{238}}U-\eqn{^{230}}Th/\eqn{^{238}}U evolution
#' diagram, a \eqn{^{234}}U/\eqn{^{238}}U-age diagram, or
#' (if \eqn{^{234}}U/\eqn{^{238}}U is assumed to be in secular equilibrium),
#' a \eqn{^{230}}Th/\eqn{^{232}}Th-\eqn{^{238}}U/\eqn{^{232}}Th diagram,
#' calculates isochron ages.
#'
#' @details
#'
#' Similar to the \code{\link{concordia}} diagram (for U-Pb data) and
#' the \code{\link{helioplot}} diagram (for U-Th-He data), the
#' evolution diagram simultaneously displays the isotopic composition
#' and age of U-series data. For carbonate data (Th-U formats 1 and
#' 2), the Th-U evolution diagram consists of a scatter plot that sets
#' out the \eqn{^{234}}U/\eqn{^{238}}U-activity ratios against the
#' \eqn{^{230}}Th/\eqn{^{238}}U-activity ratios as error ellipses, and
#' displays the initial \eqn{^{234}}U/\eqn{^{238}}U-activity ratios
#' and ages as a set of intersecting lines.  Alternatively, the
#' \eqn{^{234}}U/\eqn{^{238}}U-ratios can also be set out against the
#' \eqn{^{230}}Th-\eqn{^{234}}U-\eqn{^{238}}U-ages.  In both types of
#' evolution diagrams, \code{IsoplotR} provides the option to project
#' the raw measurements along the best fitting isochron line and
#' thereby remove the detrital \eqn{^{230}}Th-component. This
#' procedure allows a visual assessment of the degree of homogeneity
#' within a dataset, as is quantified by the MSWD.
#'
#' Neither the U-series evolution diagram, nor the
#' \eqn{^{234}}U/\eqn{^{238}}U vs. age plot is applicable to igneous
#' datasets (Th-U formats 3 and 4), in which \eqn{^{234}}U and
#' \eqn{^{238}}U are in secular equilibrium.  For such datasets,
#' \code{IsoplotR} produces an Osmond-style regression plot that is
#' decorated with a fanning set of \code{\link{isochron}} lines.
#'
#' @param x an object of class \code{ThU}
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param alpha probability cutoff for the error ellipses and
#'     confidence intervals
#' @param transform if \code{TRUE}, plots \eqn{^{234}}U/\eqn{^{238}}U
#'     vs. Th-U age.
#' @param detrital apply a detrital Th correction by projecting the
#'     compositions along an isochron?
#' @param show.numbers label the error ellipses with the grain
#'     numbers?
#' @param levels a vector with additional values to be displayed as
#'     different background colours within the error ellipses.
#' @param clabel label of the colour legend.
#' @param ellipse.col a vector of two background colours for the error
#'     ellipses. If \code{levels=NA}, then only the first colour will
#'     be used. If \code{levels} is a vector of numbers, then
#'     \code{ellipse.col} is used to construct a colour ramp.
#' @param line.col colour of the age grid
#' @param isochron fit a 3D isochron to the data?
#' @param exterr propagate the decay constant uncertainty in the
#'     isochron age?
#' @param sigdig number of significant digits for the isochron age
#' @param model if \code{isochron=TRUE}, choose one of three
#'     regression models:
#'
#' \code{1}: maximum likelihood regression, using either the modified
#' error weighted least squares algorithm of York et al. (2004) for
#' 2-dimensional data, or the Maximum Likelihood formulation of Ludwig
#' and Titterington (1994) for 3-dimensional data. These algorithms
#' take into account the analytical uncertainties and error
#' correlations, under the assumption that the scatter between the
#' data points is solely caused by the analytical uncertainty. If this
#' assumption is correct, then the MSWD value should be approximately
#' equal to one. There are three strategies to deal with the case
#' where MSWD>1. The first of these is to assume that the analytical
#' uncertainties have been underestimated by a factor
#' \eqn{\sqrt{MSWD}}.
#'
#' \code{2}: ordinary least squares regression: a second way to deal
#' with over- or underdispersed datasets is to simply ignore the
#' analytical uncertainties.
#'
#' \code{3}: maximum likelihood regression with overdispersion:
#' instead of attributing any overdispersion (MSWD > 1) to
#' underestimated analytical uncertainties (model 1), one can also
#' attribute it to the presence of geological uncertainty, which
#' manifests itself as an added (co)variance term.
#'
#' @param ... optional arguments to the generic \code{plot} function
#' @seealso \code{\link{isochron}}
#'
#' @examples
#' data(examples)
#' evolution(examples$ThU)
#'
#' dev.new()
#' evolution(examples$ThU,transform=TRUE,
#'           isochron=TRUE,model=1)
#'
#' @references Ludwig, K.R. and Titterington, D.M., 1994. Calculation
#'     of \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#'     Cosmochimica Acta, 58(22), pp.5031-5042.
#'
#' Ludwig, K.R., 2003. Mathematical-statistical treatment of data and
#'     errors for \eqn{^{230}}Th/U geochronology. Reviews in Mineralogy and
#'     Geochemistry, 52(1), pp.631-656.
#' @export
evolution <- function(x,xlim=NA,ylim=NA,alpha=0.05,transform=FALSE,
                      detrital=FALSE,show.numbers=FALSE,levels=NA,
                      clabel="",ellipse.col=c("#00FF0080","#FF000080"),
                      line.col='darksalmon',isochron=FALSE, model=1,
                      exterr=TRUE,sigdig=2,...){
    if (x$format %in% c(1,2)){
        if (transform){
            U4U8vst(x,detrital=detrital,xlim=xlim,ylim=ylim,
                    alpha=alpha,show.numbers=show.numbers,
                    levels=levels,clabel=clabel,
                    ellipse.col=ellipse.col,show.ellipses=(model!=2),
                    ...)
        } else {
            U4U8vsTh0U8(x,isochron=isochron,model=model,
                        detrital=detrital,xlim=xlim,ylim=ylim,
                        alpha=alpha,show.numbers=show.numbers,
                        levels=levels,clabel=clabel,
                        ellipse.col=ellipse.col, line.col=line.col,
                        show.ellipses=(model!=2), ...)
        }
        if (isochron){
            fit <- isochron.ThU(x,type=3,plot=FALSE,
                                exterr=exterr,model=model)
            fit$n <- length(x)
            graphics::title(evolution.title(fit,sigdig=sigdig))
        }
    } else {
        Th02vsTh0U8(x,isochron=isochron,model=model,xlim=xlim,
                    ylim=ylim,alpha=alpha,show.numbers=show.numbers,
                    exterr=exterr,sigdig=sigdig,levels=levels,
                    clabel=clabel, ellipse.col=ellipse.col,
                    line.col=line.col,...)
    }
}

U4U8vst <- function(x,detrital=FALSE,xlim=NA,ylim=NA,alpha=0.05,
                    show.numbers=FALSE,levels=NA,clabel="",
                    ellipse.col=c("#00FF0080","#FF000080"),
                    show.ellipses=TRUE,...){
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
    ellipse.cols <- set.ellipse.colours(ns=ns,levels=levels,col=ellipse.col)
    x0 <- ta0[,'t']
    y0 <- ta0[,'48_0']
    if (show.ellipses){
        for (i in 1:ns){
            diag(covmat) <- ta0[i,c('s[t]','s[48_0]')]^2
            covmat[1,2] <- ta0[i,'cov[t,48_0]']
            covmat[2,1] <- covmat[1,2]
            ell <- ellipse(x0[i],y0[i],covmat,alpha=alpha)
            graphics::polygon(ell,col=ellipse.cols[i])
            if (show.numbers) graphics::text(x0[i],y0[i],i)
            else graphics::points(x0[i],y0[i],pch=19,cex=0.25)
        }
    } else {
        plot_points(x0,y0,bg=ellipse.cols,show.numbers=show.numbers,...)
    }
    colourbar(z=levels,col=ellipse.col,clabel=clabel)
}

U4U8vsTh0U8 <- function(x,isochron=FALSE,model=1,detrital=FALSE,
                        xlim=NA,ylim=NA,alpha=0.05,
                        show.numbers=FALSE,levels=NA,clabel="",
                        ellipse.col=c("#00FF0080","#FF000080"),
                        line.col='darksalmon',show.ellipses=TRUE,...){
    ns <- length(x)
    d <- data2evolution(x,detrital=detrital)
    lim <- evolution.lines(d,xlim=xlim,ylim=ylim,...)
    if (isochron){
        fit <- isochron(x,type=2,plot=FALSE,model=model)
        b48 <- fit$par['a']
        b08 <- fit$par['A']
        e48 <- 1
        e08 <- fit$par['A'] + fit$par['B']*(e48-fit$par['a'])/fit$par['b']
        graphics::lines(c(b08,e08),c(b48,e48))
    }
    ellipse.cols <- set.ellipse.colours(ns=ns,levels=levels,col=ellipse.col)
    covmat <- matrix(0,2,2)
    x0 <- d[,'Th230U238']
    y0 <- d[,'U234U238']
    if (show.ellipses){
        for (i in 1:ns){
            diag(covmat) <- d[i,c('errTh230U238','errU234U238')]^2
            covmat[1,2] <- d[i,'cov']
            covmat[2,1] <- covmat[1,2]
            ell <- ellipse(x0[i],y0[i],covmat,alpha=alpha)
            graphics::polygon(ell,col=ellipse.cols[i])
            if (show.numbers) graphics::text(x0[i],y0[i],i)
            else graphics::points(x0[i],y0[i],pch=19,cex=0.25)
        }
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
    if (!show.ellipses)
        plot_points(x0,y0,bg=ellipse.cols,show.numbers=show.numbers,...)
    colourbar(z=levels,col=ellipse.col,clabel=clabel)
}

Th02vsTh0U8 <- function(x,isochron=FALSE,model=1,xlim=NA,ylim=NA,
                        alpha=0.05,show.numbers=FALSE,exterr=TRUE,
                        levels=NA,ellipse.col=c("#00FF0080","#FF000080"),
                        sigdig=2,line.col='darksalmon',...){
    d <- data2evolution(x,isochron=isochron)
    scatterplot(d$x,xlim=xlim,ylim=ylim,empty=TRUE)
    ticks <- c(0,1,10,20,50,100,200,300)
    X <- graphics::par('usr')[1:2]
    Y <- X
    graphics::lines(X,Y,col=line.col,...) # equilibrium line
    minY <- graphics::par('usr')[3]
    maxY <- graphics::par('usr')[4]
    if (maxY<X[2]) # add infinity symbol for equilibrium line
        graphics::text(maxY,maxY,'\U221E',pos=1,col=line.col)
    else
        graphics::text(X[2],X[2],'\U221E',pos=2,col=line.col)
    for (tick in ticks){ # plot isolines
        Y <- get.Th230Th232(tick,d$Th230Th232_0x,X)
        graphics::lines(X,Y,col=line.col,...)
        if (Y[2]<minY){
            # do nothing
        } else if (Y[2]>maxY){ # label below upper margin
            xtext <- get.U238Th232(tick,d$Th230Th232_0x,maxY)
            ytext <- maxY
            graphics::text(xtext,ytext,tick,pos=1,col=line.col)
        } else { # label to the left of the right margin
            xtext <- X[2]
            ytext <- Y[2]
            graphics::text(xtext,ytext,tick,pos=2,col=line.col)
        }
    }
    if (isochron){ # plot the data and isochron line fit
        isochron.ThU(x,type=1,plot=TRUE,show.numbers=show.numbers,
                     levels=levels,ellipse.col=ellipse.col,
                     line.col='black',exterr=exterr,sigdig=sigdig,
                     new.plot=FALSE,model=model)
    } else { # plot just the data
        scatterplot(d$x,alpha=alpha,show.numbers=show.numbers,
                    levels=levels,ellipse.col=ellipse.col,
                    new.plot=FALSE)
        xlab <- expression(paste(""^"238","U/"^"232","Th"))
        ylab <- expression(paste(""^"230","Th/"^"232","Th"))
        graphics::title(xlab=xlab,ylab=ylab)
        tit <- expression(paste("[isochrons assume ("^"230","Th/"^
                                "232","Th)"[o]^x*" = 0]"))
        graphics::mtext(tit,line=0)
    }
}


evolution.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age[1],fit$age[2:4],sigdig=sigdig)
    rounded.a0 <- roundit(fit$y0[1],fit$y0[2:4],sigdig=sigdig)
    expr1 <- quote('isochron age =')
    list1 <- list(a=rounded.age[1],
                  b=rounded.age[2],
                  c=rounded.age[3],
                  n=fit$n)
    expr2 <- quote('('^234*'U/'^238*'U)'[o]*~'=')
    list2 <- list(a=rounded.a0[1],
                  b=rounded.a0[2],
                  c=rounded.a0[3])
    if (fit$model==1 && fit$mswd>1){
        args1 <- quote(~a%+-%b~'|'~c~'|'~d~'ka'~'(n='~n~')')
        args2 <- quote(~a%+-%b~'|'~c~'|'~d)
        list1$d <- rounded.age[4]
        list2$d <- rounded.a0[4]
    } else {
        args1 <- quote(~a%+-%b~'|'~c~'ka')
        args2 <- quote(~a%+-%b~'|'~c)
    }
    call1 <- substitute(e~a,list(e=expr1,a=args1))
    line1 <- do.call(substitute,list(eval(call1),list1))
    call2 <- substitute(e~a,list(e=expr2,a=args2))
    line2 <- do.call(substitute,list(eval(call2),list2))
    if (fit$model==1 && fit$mswd>1){
        line3 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                            list(a=signif(fit$mswd,2),
                                 b=signif(fit$p.value,2)))
        graphics::mtext(line1,line=2)
        graphics::mtext(line2,line=1)
        graphics::mtext(line3,line=0)
    } else if (fit$model==2) {
        graphics::mtext(line1,line=1)
        graphics::mtext(line2,line=0)
    } else if (fit$model==3) {
        rounded.disp <- roundit(fit$w[1],fit$w[2:3],sigdig=sigdig)
        expr3 <- quote('('^232*'Th/'^238*'U)'-dispersion~'=')
        args3 <- quote(a+b-c)
        list3 <- list(a=rounded.disp[1],b=rounded.disp[3],c=rounded.disp[2])
        call3 <- substitute(e~a,list(e=expr3,a=args3))
        line3 <- do.call(substitute,list(eval(call3),list3))
        graphics::mtext(line1,line=2)
        graphics::mtext(line2,line=1)
        graphics::mtext(line3,line=0)
    }
}

evolution.lines <- function(d,xlim=NA,ylim=NA,bty='n',
                            line.col='darksalmon',...){
    nn <- 20
    maxt <- 400
    tt <- seq(from=0,to=maxt,by=50)
    nsd <- 3
    if (any(is.na(xlim)))
        max.dx <- max(d[,'Th230U238']+nsd*d[,'errTh230U238'])
    else max.dx <- xlim[2] # only used if ylim == NA
    if (any(is.na(ylim))){
        max.dy <- max(d[,'U234U238']+nsd*d[,'errU234U238'])
        a0max <- get.ThU.age(max.dx,0,max.dy,0,0,exterr=FALSE)['48_0']
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

data2evolution <- function(x,detrital=FALSE,isochron=FALSE){
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
    } else if (x$format == 2){
        xy <- subset(x$x,select=c('Th230U238','errTh230U238',
                                  'U234U238','errU234U238'))
        covariance <- x$x[,'errTh230U238']*
            x$x[,'errU234U238']*x$x[,'rhoYZ']
        out <- cbind(xy,covariance)
        colnames(out) <- labels
    } else if (x$format %in% c(3,4)){
        if (isochron){ # calculate initio 230Th from isochron intercept
            fit <- isochron.ThU(x,type=1,plot=FALSE,exterr=FALSE)
            out$x <- fit$d
            out$Th230Th232_0x <- fit$y0[1]
            out$Th230Th232_0 <- fit$a[1]
        } else { # assume no initial 230Th
            out$x <- data2york(x,type=1)
            out$Th230Th232_0x <- 0
            out$Th230Th232_0 <- 0
        }
        colnames(out$x) <- c('U238Th232','errU238Th232',
                             'Th230Th232','errTh230Th232',
                             'rho')
    }
    if (detrital){
        osmond <- data2tit.ThU(x,osmond=TRUE)
        fit <- titterington(osmond)
        out[,'U234U238'] <- out[,'U234U238'] - fit$par['b']*osmond[,'X']
        out[,'Th230U238'] <- out[,'Th230U238'] - fit$par['B']*osmond[,'X']
    }
    out
}
