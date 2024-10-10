#' @title Th-U evolution diagram
#'
#' @description Plots Th-U data on a
#' \eqn{^{234}}U/\eqn{^{238}}U-\eqn{^{230}}Th/\eqn{^{238}}U evolution
#' diagram, a \eqn{^{234}}U/\eqn{^{238}}U-age diagram, or (if
#' \eqn{^{234}}U/\eqn{^{238}}U is assumed to be in secular
#' equilibrium), a
#' \eqn{^{230}}Th/\eqn{^{232}}Th-\eqn{^{238}}U/\eqn{^{232}}Th diagram;
#' calculates isochron ages.
#'
#' @details
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
#' @param tticks time intervals of the evolution grid
#' @param aticks initial activity ratio ticks of the evolution grid
#' @param oerr indicates whether the analytical uncertainties of the
#'     output are reported in the plot title as:
#' 
#' \code{1}: 1\eqn{\sigma} absolute uncertainties.
#' 
#' \code{2}: 2\eqn{\sigma} absolute uncertainties.
#' 
#' \code{3}: absolute (1-\eqn{\alpha})\% confidence intervals, where
#' \eqn{\alpha} equales the value that is stored in
#' \code{settings('alpha')}.
#'
#' \code{4}: 1\eqn{\sigma} relative uncertainties (\eqn{\%}).
#' 
#' \code{5}: 2\eqn{\sigma} relative uncertainties (\eqn{\%}).
#'
#' \code{6}: relative (1-\eqn{\alpha})\% confidence intervals, where
#' \eqn{\alpha} equales the value that is stored in
#' \code{settings('alpha')}.
#' @param transform if \code{TRUE}, plots \eqn{^{234}}U/\eqn{^{238}}U
#'     vs. Th-U age.
#' @param Th0i initial \eqn{^{230}}Th correction.
#'
#' \code{0}: no correction
#'
#' \code{1}: if \code{x$format} is \code{1} or \code{2}, project the
#' data along an isochron fit. If \code{x$format} is \code{3} or
#' \code{4}, infer the initial \eqn{^{230}}Th/\eqn{^{238}}U activity
#' ratio from the isochron.
#'
#' \code{2}: if \code{x$format} is \code{1} or \code{2}, correct the
#' data using the measured present day \eqn{^{230}}Th/\eqn{^{238}}U,
#' \eqn{^{232}}Th/\eqn{^{238}}U and \eqn{^{234}}U/\eqn{^{238}}U
#' activity ratios in the detritus. If \code{x$format} is \code{3} or
#' \code{4}, anchor the isochrons to the equiline, based on the
#' measured \eqn{^{238}}U/\eqn{^{232}}Th activity ratio of the whole
#' rock, as stored in \code{x} by the \code{read.data()} function.
#'
#' \code{3}: correct the data using an assumed initial
#' \eqn{^{230}}Th/\eqn{^{232}}Th-ratio for the detritus (only relevant
#' if \code{x$format} is \code{1} or \code{2}).
#'
#' @param show.numbers label the error ellipses with the grain
#'     numbers?
#' @param levels a vector with additional values to be displayed as
#'     different background colours within the error ellipses.
#' @param clabel label of the colour legend.
#' @param ellipse.fill fill colour for the error ellipses. This can
#'     either be a single colour or multiple colours to form a colour
#'     ramp. Examples:
#'
#' a single colour: \code{rgb(0,1,0,0.5)}, \code{'#FF000080'},
#' \code{'white'}, etc.;
#'
#' multiple colours: \code{c(rbg(1,0,0,0.5)},
#' \code{rgb(0,1,0,0.5))}, \code{c('#FF000080','#00FF0080')},
#' \code{c('blue','red')}, \code{c('blue','yellow','red')}, etc.;
#'
#' a colour palette: \code{rainbow(n=100)},
#' \code{topo.colors(n=100,alpha=0.5)}, etc.; or
#'
#' a reversed palette: \code{rev(topo.colors(n=100,alpha=0.5))},
#' etc.
#'
#' For empty ellipses, set \code{ellipse.fill=NA}
#' 
#' @param ellipse.stroke the stroke colour for the error
#'     ellipses. Follows the same formatting guidelines as
#'     \code{ellipse.fill}
#' @param line.col colour of the age grid
#' @param isochron fit an isochron to the data?
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
#' uncertainties have been underestipmated by a factor
#' \eqn{\sqrt{MSWD}}.
#'
#' \code{2}: total least squares regression: a second way to deal
#' with over- or underdispersed datasets is to simply ignore the
#' analytical uncertainties.
#'
#' \code{3}: maximum likelihood regression with overdispersion:
#' instead of attributing any overdispersion (MSWD > 1) to
#' underestimated analytical uncertainties (model 1), one can also
#' attribute it to the presence of geological uncertainty, which
#' manifests itself as an added (co)variance term.
#'
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the isochron age calculation.
#' @param omit.fill fill colour that should be used for the omitted
#'     aliquots.
#' @param omit.stroke stroke colour that should be used for the omitted
#'     aliquots.
#' @param ... optional arguments to the generic \code{plot} function
#' @seealso \code{\link{isochron}}
#'
#' @examples
#' attach(examples)
#' evolution(ThU)
#'
#' dev.new()
#' evolution(ThU,transform=TRUE,isochron=TRUE,model=1)
#'
#' @references Ludwig, K.R. and Titterington, D.M., 1994. Calculation
#'     of \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#'     Cosmochimica Acta, 58(22), pp.5031-5042.
#'
#' Ludwig, K.R., 2003. Mathematical-statistical treatment of data and
#'     errors for \eqn{^{230}}Th/U geochronology. Reviews in Mineralogy and
#'     Geochemistry, 52(1), pp.631-656.
#' @export
evolution <- function(x,xlim=NULL,ylim=NULL,tticks=NULL,aticks=NULL,oerr=3,
                      transform=FALSE,Th0i=0,show.numbers=FALSE,levels=NA,
                      clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                      ellipse.stroke='black',line.col='darksalmon',
                      isochron=FALSE,model=1,exterr=FALSE,sigdig=2,
                      hide=NULL,omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    if (x$format %in% c(1,2)){
        if (transform){
            U4U8vst(x,Th0i=Th0i,xlim=xlim,ylim=ylim,oerr=oerr,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,
                    show.ellipses=(model!=2),hide=hide,omit=omit,
                    omit.fill=omit.fill,omit.stroke=omit.stroke,...)
        } else {
            U4U8vsTh0U8(x,isochron=isochron,model=model,xlim=xlim,
                        ylim=ylim,tticks=tticks,aticks=aticks,
                        oerr=oerr,Th0i=Th0i,show.numbers=show.numbers,
                        levels=levels,clabel=clabel,
                        ellipse.fill=ellipse.fill,
                        ellipse.stroke=ellipse.stroke,
                        line.col=line.col,show.ellipses=(model!=2),
                        hide=hide,omit=omit,omit.fill=omit.fill,
                        omit.stroke=omit.stroke,...)
        }
        if (isochron){
            fit <- isochron.ThU(x,type=3,plot=FALSE,exterr=exterr,
                                model=model,hide=hide,omit=omit,oerr=oerr)
            fit$n <- length(x)-length(hide)-length(omit)
            graphics::title(evolution_title(fit,sigdig=sigdig,oerr=oerr))
        }
    } else {
        Th02vsU8Th2(x,isochron=isochron,model=model,Th0i=Th0i,
                    xlim=xlim,ylim=ylim,tticks=tticks,aticks=aticks,
                    oerr=oerr,show.numbers=show.numbers,exterr=exterr,
                    sigdig=sigdig,levels=levels,clabel=clabel,
                    ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,line.col=line.col,
                    hide=hide,omit=omit,omit.fill=omit.fill,
                    omit.stroke=omit.stroke,...)
    }
}

U4U8vst <- function(x,Th0i=0,xlim=NULL,ylim=NULL,oerr=3,
                    show.numbers=FALSE,levels=NA,clabel="",
                    ellipse.fill=c("#00FF0080","#FF000080"),
                    ellipse.stroke='black',show.ellipses=TRUE,
                    hide=NULL,omit=NULL,omit.fill=NA,
                    omit.stroke='grey',...){
    ns <- length(x)
    plotit <- (1:ns)%ni%hide
    ta0 <- get_ThU_age_corals(x,exterr=FALSE,cor=FALSE,Th0i=Th0i)
    nsd <- 3
    if (is.null(xlim))
        xlim <- range(c(ta0[plotit,'t']-nsd*ta0[plotit,'s[t]'],
                        ta0[plotit,'t']+nsd*ta0[plotit,'s[t]']))
    if (is.null(ylim))
        ylim <- range(c(ta0[plotit,'48_0']-nsd*ta0[plotit,'s[48_0]'],
                        ta0[plotit,'48_0']+nsd*ta0[plotit,'s[48_0]']))
    x.lab <- 'Age [ka]'
    y.lab <- expression(paste("("^"234","U/"^"238","U)"[0]))
    graphics::plot(xlim,ylim,type='n',bty='n',xlab=x.lab,ylab=y.lab)
    d <- ta0
    colnames(d) <- c('X','sX','Y','sY','rXY')
    d[,'rXY'] <- ta0[,'cov[t,48_0]']/(ta0[,'s[t]']*ta0[,'s[48_0]'])
    scatterplot(d,oerr=oerr,show.numbers=show.numbers,
                show.ellipses=show.ellipses,levels=levels,
                clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,add=TRUE, hide=hide,
                omit=omit,omit.fill=omit.fill,omit.stroke=omit.stroke,...)
}

U4U8vsTh0U8 <- function(x,isochron=FALSE,model=1,Th0i=0,
                        xlim=NULL,ylim=NULL,oerr=3,
                        show.numbers=FALSE,levels=NA,clabel="",
                        ellipse.fill=c("#00FF0080","#FF000080"),
                        ellipse.stroke='black',line.col='darksalmon',
                        show.ellipses=TRUE,hide=NULL,omit=NULL,
                        omit.fill=NA,omit.stroke='grey',...){
    ns <- length(x)
    plotit <- (1:ns)%ni%hide
    calcit <- (1:ns)%ni%c(hide,omit)
    y <- data2evolution(x,Th0i=Th0i,omit4c=unique(c(hide,omit)))
    d2plot <- subset(y,subset=plotit)
    lim <- evolution_lines(d2plot,xlim=xlim,ylim=ylim,...)
    if (isochron){
        fit <- isochron(x,type=3,plot=FALSE,
                        model=model,hide=hide,omit=omit)
        b48 <- fit$par['b']
        b08 <- fit$par['B']
        initial <- matrix(0,1,5)
        initial[1] <- b08
        initial[2] <- sqrt(fit$cov['b','b'])
        initial[3] <- b48
        initial[4] <- sqrt(fit$cov['B','B'])
        initial[5] <- fit$cov['b','B']/(initial[2]*initial[4])
        scatterplot(initial,oerr=oerr,line.col='black',add=TRUE,
                    ellipse.fill=grDevices::rgb(1,1,1,0.85))
        e48 <- 1
        e08 <- b08 + fit$par['A']*(e48-b48)/fit$par['a']
        graphics::lines(c(b08,e08),c(b48,e48))
    }
    pdat <- y[,c('Th230U238','sTh230U238',
                 'U234U238','sU234U238','rYZ'),drop=FALSE]
    scatterplot(pdat,oerr=oerr,show.numbers=show.numbers,
                show.ellipses=show.ellipses,levels=levels,
                clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,add=TRUE,hide=hide,
                omit=omit,omit.fill=omit.fill,omit.stroke=omit.stroke,...)
    colourbar(z=levels[calcit],fill=ellipse.fill,
              stroke=ellipse.stroke,clabel=clabel)
}

Th02vsU8Th2 <- function(x,isochron=FALSE,model=1,Th0i=0,xlim=NULL,
                        ylim=NULL,tticks=NULL,oerr=3,show.numbers=FALSE,
                        exterr=FALSE,clabel="",levels=NA,
                        ellipse.fill=c("#00FF0080","#FF000080"),
                        ellipse.stroke='black',sigdig=2,
                        line.col='darksalmon',hide=NULL,omit=NULL,
                        omit.fill=NA,omit.stroke='grey',...){
    ns <- length(x)
    plotit <- (1:ns)%ni%hide
    calcit <- (1:ns)%ni%c(hide,omit)
    d <- data2evolution(x,omit4c=unique(c(hide,omit)))
    d2plot <- subset(d,subset=plotit)
    scatterplot(d2plot,xlim=xlim,ylim=ylim,empty=TRUE)
    if (is.null(tticks)) tticks <- c(0,5,10,20,50,100,200,Inf)
    nt <- length(tticks)
    X <- graphics::par('usr')[1:2]
    Y <- X
    minY <- graphics::par('usr')[3]
    maxY <- graphics::par('usr')[4]
    l0 <- lambda('Th230')[1]
    if (isochron|Th0i==1){
        fit <- isochron.ThU(x,type=1,plot=FALSE,exterr=FALSE,
                            hide=hide,omit=omit,omit.fill=omit.fill,
                            omit.stroke=omit.stroke)
        anchor <- cbind(0,fit$y0[1]*exp(-l0*tticks))
    } else if (Th0i==2){
        anchor <- matrix(x$U8Th2,nt-1,2)
        tticks <- tticks[is.finite(tticks)]
    } else {
        anchor <- matrix(0,nt,2)
    }
    for (i in seq_along(tticks)){ # plot isolines
        if (is.finite(tticks[i])) ticktext <- tticks[i]
        else ticktext <- expression(infinity)
        slope <- 1-exp(-l0*tticks[i])
        Y <- anchor[i,2] + slope*(X-anchor[i,1])
        graphics::lines(X,Y,col=line.col,...)
        if (Y[2]<minY){
            # do nothing
        } else if (Y[2]>maxY){ # label below upper margin
            xtext <- anchor[i,1] + (maxY-anchor[i,2])/slope
            ytext <- maxY
            graphics::text(xtext,ytext,ticktext,pos=1,col=line.col)
        } else { # label to the left of the right margin
            xtext <- X[2]
            ytext <- Y[2]
            graphics::text(xtext,ytext,ticktext,pos=2,col=line.col)
        }
    }
    if (isochron){ # plot the data and isochron line fit
        isochron.ThU(x,type=1,oerr=oerr,plot=TRUE,show.numbers=show.numbers,
                     levels=levels,ellipse.fill=ellipse.fill,
                     ellipse.stroke=ellipse.stroke,
                     line.col='black',exterr=exterr,sigdig=sigdig,
                     add=TRUE,model=model,hide=hide,
                     omit=omit,omit.fill=omit.fill,omit.stroke=omit.stroke)
    } else { # plot just the data
        scatterplot(d,oerr=oerr,show.numbers=show.numbers,
                    levels=levels,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,
                    add=TRUE,hide=hide,omit=omit,
                    omit.fill=omit.fill,omit.stroke=omit.stroke)
        xlab <- expression(paste(""^"238","U/"^"232","Th"))
        ylab <- expression(paste(""^"230","Th/"^"232","Th"))
        graphics::title(xlab=xlab,ylab=ylab)
        if (Th0i==0){
            tit <- expression(paste("[isochrons assume ("^"230","Th/"^
                                    "232","Th)"[i]*" = 0]"))
            mymtext(tit,line=0,...)
        }
        if (Th0i==2){ # add equiline
            middle <- max(X[1],minY)/2 + min(X[2],maxY)/2
            graphics::text(middle,middle,'1:1',pos=3)
            graphics::lines(X,X)
            graphics::points(x$U8Th2,x$U8Th2,pch=16)
        }
    }
    invisible(colourbar(z=levels[calcit],fill=ellipse.fill,
                        stroke=ellipse.stroke,clabel=clabel))
}

evolution_title <- function(fit,sigdig=2,oerr=3,...){
    content <- list()
    content[[1]] <- maintit(x=fit$age[1],sx=fit$age[-1],
                            sigdig=sigdig,n=fit$n,
                            oerr=oerr,prefix='isochron age =',
                            units=' ka',df=fit$df)
    content[[2]] <- maintit(x=fit$y0[1],sx=fit$y0[-1],sigdig=sigdig,ntit='',
                            oerr=oerr,prefix=fit$y0label,units='',df=fit$df)
    if (fit$model==1){
        content[[3]] <- mswdtit(mswd=fit$mswd,p=fit$p.value,sigdig=sigdig)
    }
    if (fit$model==3) {
        prefix <- quote('('^234*'U/'^238*'U)'[a]*'-dispersion = ')
        content[[3]] <- disptit(w=fit$disp[1],sw=fit$disp[-1],sigdig=sigdig,
                                    oerr=oerr,units='',prefix=prefix)
    }
    nl <- length(content)
    for (i in 1:nl){
        mymtext(content[[i]],line=nl-i)
    }
}

evolution_lines <- function(d,xlim=NULL,ylim=NULL,
                            tticks=NULL,aticks=NULL,
                            bty='n',line.col='darksalmon',...){
    if (is.null(tticks)) tticks <- seq(0,500,by=50)
    nsd <- 3
    if (is.null(xlim)){
        min.dx <- 0
        max.dx <- max(d[,'Th230U238']+nsd*d[,'sTh230U238'])
    } else {
        min.dx <- xlim[1]
        max.dx <- xlim[2]
    }
    if (is.null(ylim)){
        min.dy <- min(d[,'U234U238']-nsd*d[,'sU234U238'])
        max.dy <- max(d[,'U234U238']+nsd*d[,'sU234U238'])
        a01 <- get_ThU_age(min.dx,0,min.dy,0,0,exterr=FALSE)['48_0']
        a02 <- get_ThU_age(min.dx,0,max.dy,0,0,exterr=FALSE)['48_0']
        a03 <- get_ThU_age(max.dx,0,min.dy,0,0,exterr=FALSE)['48_0']
        a04 <- get_ThU_age(max.dx,0,max.dy,0,0,exterr=FALSE)['48_0']
        a0min <- min(a01,a02,a03,a04)
        a0max <- max(a01,a02,a03,a04)
    } else {
        a0min <- ylim[1]
        a0max <- ylim[2]
    }
    if ((a0max-a0min)<0.01){
        a0min <- 0
        a0max <- 1.5
    }
    if (is.null(xlim)){
        xlim <- range(pretty(c(min.dx,
                               min(get_Th230U238_ratio(tticks,a0min)),
                               max(get_Th230U238_ratio(tticks,a0max)),
                               max.dx)))
    }
    if (is.null(aticks)) aticks <- pretty(c(a0min,a0max))
    if (is.null(ylim)){
        ylim <- range(aticks)
        aticks <- aticks[aticks>0]
    } else {
        aticks <- c(aticks[aticks>0 & aticks<ylim[2]],ylim[2])
    }
    x.lab <- expression(paste(""^"230","Th/"^"238","U"))
    y.lab <- expression(paste(""^"234","U/"^"238","U"))
    graphics::plot(xlim,ylim,type='n',xlab=x.lab,ylab=y.lab,bty=bty,...)
    na0 <- length(aticks)
    nn <- 100
    for (i in 1:na0){
        tt <- seq(0,2000,length.out=nn)
        x <- get_Th230U238_ratio(tt,aticks[i])
        y <- get_U234U238_ratio(tt,aticks[i])
        graphics::lines(x,y,col=line.col)
    }
    for (i in 1:nn){
        x <- get_Th230U238_ratio(tticks[i],aticks)
        y <- get_U234U238_ratio(tticks[i],aticks)
        x0 <- get_Th230U238_ratio(tticks[i],1)
        graphics::lines(c(x0,x),c(1,y),col=line.col)
        if (is.finite(tticks[i])) ticktext <- tticks[i]
        else ticktext <- expression(infinity)
        graphics::text(x[na0],y[na0],ticktext,pos=4,col=line.col)
    }
    rbind(xlim,ylim)
}

data2evolution <- function(x,Th0i=0,omit4c=NULL){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format %in% c(1,2)){
        td <- data2tit(x,osmond=TRUE,generic=FALSE) # 2/8 - 4/8 - 0/8
        if (Th0i==1){
            out <- Th230correction_isochron(td,omit4c=omit4c)
        } else if (Th0i==2){
            out <- Th230correction_measured_detritus(td,Th02U48=x$Th02U48)
        } else if (Th0i==3){
            tt <- get_ThU_age_corals(x,Th0i=Th0i)[,'t']
            out <- Th230correction_assumed_detritus(td,age=tt,Th02i=x$Th02i)
        } else {
            out <- td
        }
    } else if (x$format %in% c(3,4)){
        out <- data2york(x,type=1,generic=FALSE) # 8/2 - 0/2
    }
    out
}

Th230correction_isochron <- function(td,omit4c=NULL){
    fit <- titterington(clear(td,omit4c))
    out <- td
    out[,'U234U238'] <- td[,'U234U238'] - fit$par['b']*td[,'Th232U238']
    out[,'Th230U238'] <- td[,'Th230U238'] - fit$par['B']*td[,'Th232U238']
    out
}
Th230correction_assumed_detritus <- function(td,age=Inf,Th02i=c(0,0)){
    out <- td
    l0 <- lambda('Th230')[1]
    A <- Th02i[1]*exp(-l0*age)*td[,'Th232U238']
    out[,'Th230U238'] <- td[,'Th230U238'] - A
    dA.dTh02i <- -exp(-l0*age)*td[,'Th232U238']
    dA.Th2U8 <- -Th02i[1]*exp(-l0*age)
    sA <- errorprop1x2(dA.dTh02i,dA.Th2U8,
                       Th02i[2]^2,td[,'sTh232U238']^2,0)
    out[,'sTh230U238'] <- sqrt(td[,'sTh230U238']^2 + sA^2)         
    out
}
Th230correction_measured_detritus <- function(td,Th02U48=c(0,0,1e6,0,0,0,0,0,0)){
    X1 <- td[,'Th232U238']
    sX1 <- td[,'sTh232U238']
    Y1 <- td[,'U234U238']
    sY1 <- td[,'sU234U238']
    Z1 <- td[,'Th230U238']
    sZ1 <- td[,'sTh230U238']
    rX1Y1 <- td[,'rXY']
    rX1Z1 <- td[,'rXZ']
    rY1Z1 <- td[,'rYZ']
    covX1Y1 <- rX1Y1*sX1*sY1
    covX1Z1 <- rX1Z1*sX1*sZ1
    covY1Z1 <- rY1Z1*sY1*sZ1
    X2 <- Th02U48[3]
    Y2 <- Th02U48[5]
    Z2 <- Th02U48[1]
    sX2 <- Th02U48[4]^2
    sY2 <- Th02U48[6]^2
    sZ2 <- Th02U48[2]^2
    covX2Y2 <- Th02U48[9]*Th02U48[4]*Th02U48[6]
    covX2Z2 <- Th02U48[7]*Th02U48[2]*Th02U48[4]
    covY2Z2 <- Th02U48[8]*Th02U48[2]*Th02U48[6]
    b <- (Y2-Y1)/(X2-X1)
    B <- (Z2-Z1)/(X2-X1)
    r1 <- X1/(X2-X1)
    r2 <- X2/(X2-X1)
    a <- Y1 - b * X1
    A <- Z1 - B * X1
    sa <- sqrt( (b^2)*((r2*sX1)^2+(r1*sX2)^2) +
                (r2*sY1)^2 + (r1*sY2)^2 -
                2*r2*b*covX1Y1 - 2*(r1^2)*covX2Y2
               )
    sA <- sqrt( (B^2)*((r2*sX1)^2+(r1*sX2)^2) +
                (r2*sZ1)^2 + (r1*sZ2)^2 -
                2*r2*B*covX1Z1 - 2*(r1^2)*covX2Z2
               )
    covaA <- b*B*((r2*sX1)^2+(r1*sX2)^2) +
        (r2^2)*(covY1Z1-B*covX1Y1-b*covX1Z1) +
        (r1^2)*(covY2Z2-B*covX2Y2-b*covX2Z2)
    out <- td
    out[,'Th230U238'] <- A
    out[,'sTh230U238'] <- sA
    out[,'U234U238'] <- a
    out[,'sU234U238'] <- sa
    out[,'rYZ'] <- covaA/(sA*sa)
    out
}
