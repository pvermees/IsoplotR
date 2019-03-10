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
#' @param detritus detrital \eqn{^{230}}Th correction (only applicable
#'     when \code{x$format} is \code{2} or \code{3}.
#'
#' \code{0}: no correction
#'
#' \code{1}: project the data along an isochron fit
#'
#' \code{2}: correct the data using an assumed initial
#' \eqn{^{230}}Th/\eqn{^{232}}Th-ratio for the detritus.
#'
#' \code{3}: correct the data using the measured present day
#' \eqn{^{230}}Th/\eqn{^{238}}U, \eqn{^{232}}Th/\eqn{^{238}}U and
#' \eqn{^{234}}U/\eqn{^{238}}U-ratios in the detritus.
#' 
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
#' uncertainties have been underestipmated by a factor
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
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the isochron age calculation.
#' @param omit.col colour that should be used for the omitted
#'     aliquots.
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
                      detritus=0,show.numbers=FALSE,levels=NA,
                      clabel="",ellipse.col=c("#00FF0080","#FF000080"),
                      line.col='darksalmon',isochron=FALSE,model=1,
                      exterr=TRUE,sigdig=2,hide=NULL,omit=NULL,
                      omit.col=NA,...){
    if (x$format %in% c(1,2)){
        if (transform){
            U4U8vst(x,detritus=detritus,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.col=ellipse.col,
                    show.ellipses=(model!=2),hide=hide,omit=omit,
                    omit.col=omit.col,...)
        } else {
            U4U8vsTh0U8(x,isochron=isochron,model=model,xlim=xlim,
                        ylim=ylim,alpha=alpha,detritus=detritus,
                        show.numbers=show.numbers,levels=levels,
                        clabel=clabel,ellipse.col=ellipse.col,
                        line.col=line.col,show.ellipses=(model!=2),
                        hide=hide,omit=omit,omit.col=omit.col,...)
        }
        if (isochron){
            fit <- isochron.ThU(x,type=3,plot=FALSE,exterr=exterr,
                                model=model,hide=hide,omit=omit)
            fit$n <- length(x)-length(hide)-length(omit)
            graphics::title(evolution.title(fit,sigdig=sigdig))
        }
    } else {
        Th02vsU8Th2(x,isochron=isochron,model=model,xlim=xlim,
                    ylim=ylim,alpha=alpha,show.numbers=show.numbers,
                    exterr=exterr,sigdig=sigdig,levels=levels,
                    clabel=clabel,ellipse.col=ellipse.col,
                    line.col=line.col,hide=hide,omit=omit,
                    omit.col=omit.col,...)
    }
}

U4U8vst <- function(x,detritus=0,xlim=NA,ylim=NA,alpha=0.05,
                    show.numbers=FALSE,levels=NA,clabel="",
                    ellipse.col=c("#00FF0080","#FF000080"),
                    show.ellipses=TRUE,hide=NULL,omit=NULL,omit.col=NA,...){
    ns <- length(x)
    plotit <- (1:ns)%ni%hide
    calcit <- (1:ns)%ni%c(hide,omit)
    ta0 <- get.ThU.age.corals(x,exterr=FALSE,cor=FALSE,detritus=detritus)
    nsd <- 3
    if (any(is.na(xlim)))
        xlim <- range(c(ta0[plotit,'t']-nsd*ta0[plotit,'s[t]'],
                        ta0[plotit,'t']+nsd*ta0[plotit,'s[t]']))
    if (any(is.na(ylim)))
        ylim <- range(c(ta0[plotit,'48_0']-nsd*ta0[plotit,'s[48_0]'],
                        ta0[plotit,'48_0']+nsd*ta0[plotit,'s[48_0]']))
    x.lab <- 'Age [ka]'
    y.lab <- expression(paste("("^"234","U/"^"238","U)"[o]))
    graphics::plot(xlim,ylim,type='n',bty='n',xlab=x.lab,ylab=y.lab)
    d <- ta0
    colnames(d) <- c('X','sX','Y','sY','rXY')
    d[,'rXY'] <- ta0[,'cov[t,48_0]']/(ta0[,'s[t]']*ta0[,'s[48_0]'])
    scatterplot(d,alpha=alpha,show.numbers=show.numbers,
                show.ellipses=show.ellipses,levels=levels,
                clabel=clabel,ellipse.col=ellipse.col,add=TRUE,
                hide=hide,omit=omit,omit.col=omit.col,...)
}

U4U8vsTh0U8 <- function(x,isochron=FALSE,model=1,detritus=0,
                        xlim=NA,ylim=NA,alpha=0.05,
                        show.numbers=FALSE,levels=NA,clabel="",
                        ellipse.col=c("#00FF0080","#FF000080"),
                        line.col='darksalmon',show.ellipses=TRUE,
                        hide=NULL,omit=NULL,omit.col=NA,...){
    ns <- length(x)
    plotit <- (1:ns)%ni%hide
    calcit <- (1:ns)%ni%c(hide,omit)
    y <- data2evolution(x,detritus=detritus)
    d2plot <- subset(y,subset=plotit)
    lim <- evolution.lines(d2plot,xlim=xlim,ylim=ylim,...)
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
        scatterplot(initial,alpha=alpha,
                    ellipse.col=grDevices::rgb(1,1,1,0.85),
                    line.col='black',add=TRUE)
        e48 <- 1
        e08 <- b08 + fit$par['A']*(e48-b48)/fit$par['a']
        graphics::lines(c(b08,e08),c(b48,e48))
    }
    pdat <- y[,c('Th230U238','sTh230U238','U234U238','sU234U238','rYZ')]
    scatterplot(pdat,alpha=alpha,show.numbers=show.numbers,
                show.ellipses=show.ellipses,levels=levels,
                clabel=clabel,ellipse.col=ellipse.col,add=TRUE,
                hide=hide,omit=omit,omit.col=omit.col,...)
    colourbar(z=levels[calcit],col=ellipse.col,clabel=clabel)
}

Th02vsU8Th2 <- function(x,isochron=FALSE,model=1,xlim=NA,ylim=NA,
                        alpha=0.05,show.numbers=FALSE,exterr=TRUE,
                        clabel="",levels=NA,
                        ellipse.col=c("#00FF0080","#FF000080"),
                        sigdig=2,line.col='darksalmon',
                        hide=NULL,omit=NULL,omit.col=NA,...){
    ns <- length(x)
    plotit <- (1:ns)%ni%hide
    calcit <- (1:ns)%ni%c(hide,omit)
    d <- data2evolution(x)
    d2plot <- subset(d,subset=plotit)
    scatterplot(d2plot,xlim=xlim,ylim=ylim,empty=TRUE)
    ticks <- c(0,1,10,20,50,100,200,300)
    X <- graphics::par('usr')[1:2]
    Y <- X
    graphics::lines(X,Y,col=line.col,...) # equilibrium line
    minY <- graphics::par('usr')[3]
    maxY <- graphics::par('usr')[4]
    if (isochron){
        fit <- isochron.ThU(x,type=1,plot=FALSE,exterr=FALSE,
                            hide=hide,omit=omit,omit.col=omit.col)
        Th230Th232_0x <- fit$y0[1]
    } else {
        Th230Th232_0x <- 0
    }
    if (maxY<X[2]) # add infinity symbol for equilibrium line
        graphics::text(maxY,maxY,'\U221E',pos=1,col=line.col)
    else
        graphics::text(X[2],X[2],'\U221E',pos=2,col=line.col)
    for (tick in ticks){ # plot isolines
        Y <- get.Th230Th232(tick,Th230Th232_0x,X)
        graphics::lines(X,Y,col=line.col,...)
        if (Y[2]<minY){
            # do nothing
        } else if (Y[2]>maxY){ # label below upper margin
            xtext <- get.U238Th232(tick,Th230Th232_0x,maxY)
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
                     add=TRUE,model=model,hide=hide,
                     omit=omit,omit.col=omit.col)
    } else { # plot just the data
        scatterplot(d,alpha=alpha,show.numbers=show.numbers,
                    levels=levels,ellipse.col=ellipse.col,
                    add=TRUE,hide=hide,omit=omit,
                    omit.col=omit.col)
        xlab <- expression(paste(""^"238","U/"^"232","Th"))
        ylab <- expression(paste(""^"230","Th/"^"232","Th"))
        graphics::title(xlab=xlab,ylab=ylab)
        tit <- expression(paste("[isochrons assume ("^"230","Th/"^
                                "232","Th)"[o]^x*" = 0]"))
        mymtext(tit,line=0,...)
    }
    colourbar(z=levels[calcit],col=ellipse.col,clabel=clabel)
}

evolution.title <- function(fit,sigdig=2,...){
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
        args1 <- quote(~a%+-%b~'|'~c~'|'~d~'ka'~'(n='*n*')')
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
        line3 <- substitute('MSWD ='~a~', p('*chi^2*')='~b,
                            list(a=signif(fit$mswd,2),
                                 b=signif(fit$p.value,2)))
        mymtext(line1,line=2,...)
        mymtext(line2,line=1,...)
        mymtext(line3,line=0,...)
    } else if (fit$model==2) {
        mymtext(line1,line=1,...)
        mymtext(line2,line=0,...)
    } else if (fit$model==3) {
        rounded.disp <- roundit(fit$w[1],fit$w[2:3],sigdig=sigdig)
        expr3 <- quote('('^232*'Th/'^238*'U)'-dispersion~'=')
        args3 <- quote(a+b-c)
        list3 <- list(a=rounded.disp[1],b=rounded.disp[3],c=rounded.disp[2])
        call3 <- substitute(e~a,list(e=expr3,a=args3))
        line3 <- do.call(substitute,list(eval(call3),list3))
        mymtext(line1,line=2,...)
        mymtext(line2,line=1,...)
        mymtext(line3,line=0,...)
    }
}

evolution.lines <- function(d,xlim=NA,ylim=NA,bty='n',
                            line.col='darksalmon',...){
    nn <- 20
    maxt <- 400
    tt <- seq(from=0,to=maxt,by=50)
    nsd <- 3
    if (any(is.na(xlim))){
        min.dx <- 0
        max.dx <- max(d[,'Th230U238']+nsd*d[,'sTh230U238'])
    } else {
        min.dx <- xlim[1]
        max.dx <- xlim[2]
    }
    if (any(is.na(ylim))){
        min.dy <- min(d[,'U234U238']-nsd*d[,'sU234U238'])
        max.dy <- max(d[,'U234U238']+nsd*d[,'sU234U238'])
        a01 <- get.ThU.age(min.dx,0,min.dy,0,0,exterr=FALSE)['48_0']
        a02 <- get.ThU.age(min.dx,0,max.dy,0,0,exterr=FALSE)['48_0']
        a03 <- get.ThU.age(max.dx,0,min.dy,0,0,exterr=FALSE)['48_0']
        a04 <- get.ThU.age(max.dx,0,max.dy,0,0,exterr=FALSE)['48_0']
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
    if (any(is.na(xlim))){
        xlim <- range(pretty(c(min(get.Th230U238(tt,a0min)),
                               max(get.Th230U238(tt,a0max)))))
    }
    a0 <- pretty(c(a0min,a0max))
    if (any(is.na(ylim))){
        ylim <- range(a0)
        a0 <- a0[a0>0]
    } else {
        a0 <- c(a0[a0>0 & a0<ylim[2]],ylim[2])
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

data2evolution <- function(x,detritus=0){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format %in% c(1,2)){
        td <- data2tit(x,osmond=TRUE,generic=FALSE) # 2/8 - 4/8 - 0/8
        out <- Th230correction(td,option=detritus,dat=x)
    } else if (x$format %in% c(3,4)){
        out <- data2york(x,type=1) # 8/2 - 0/2
        covariance <- out[,'sX']*out[,'sY']*out[,'rXY']
        out[,5] <- covariance
        colnames(out) <- c('U238Th232','errU238Th232',
                           'Th230Th232','errTh230Th232',
                           'cov')
    }
    out
}

# x = table with 'Th230U238','errTh230U238', 'U234U238','errU234U238'
#                (and 'Th232U238','errTh232U238' if option==2)
Th230correction <- function(x,option=0,dat=NA){
    out <- x
    if (option==1){
        out <- Th230correction.isochron(x,dat=dat)
    } else if (option==2){
        tt <- get.ThU.age.corals(dat,detritus=2)[,'t']
        out <- Th230correction.assumed.detritus(x,age=tt,Th02=dat$Th02)
    } else if (option==3){
        out <- Th230correction.measured.detritus(dat)
    }
    out
}
Th230correction.isochron <- function(x,dat){
    osmond <- data2tit.ThU(dat,osmond=TRUE)
    fit <- titterington(osmond)
    out <- x
    out[,'U234U238'] <- x[,'U234U238'] - fit$par['b']*osmond[,'X']
    out[,'Th230U238'] <- x[,'Th230U238'] - fit$par['B']*osmond[,'X']
    out
}
Th230correction.assumed.detritus <- function(x,age=Inf,Th02=c(0,0)){
    out <- x
    l0 <- lambda('Th230')[1]
    A <- Th02[1]*exp(-l0*age)*x[,'Th232U238']
    out[,'Th230U238'] <- x[,'Th230U238'] - A
    dA.dTh02 <- -exp(-l0*age)*x[,'Th232U238']
    dA.Th2U8 <- -Th02[1]*exp(-l0*age)
    sA <- errorprop1x2(dA.dTh02,dA.Th2U8,
                       Th02[2]^2,x[,'sTh232U238']^2,0)
    out[,'sTh230U238'] <- sqrt(x[,'sTh230U238']^2 + sA^2)         
    out
}
Th230correction.measured.detritus <- function(x){
    Th02U48 <- x$Th02U48
    osmond <- data2tit.ThU(x,osmond=TRUE,generic=FALSE) # 2/8 - 4/8 - 0/8
    X1 <- osmond[,'Th232U238']
    sX1 <- osmond[,'sTh232U238']
    Y1 <- osmond[,'U234U238']
    sY1 <- osmond[,'sU234U238']
    Z1 <- osmond[,'Th230U238']
    sZ1 <- osmond[,'sTh230U238']
    rX1Y1 <- osmond[,'rXY']
    rX1Z1 <- osmond[,'rXZ']
    rY1Z1 <- osmond[,'rYZ']
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
    out <- osmond
    out[,'Th230U238'] <- A
    out[,'sTh230U238'] <- sA
    out[,'U234U238'] <- a
    out[,'sU234U238'] <- sa
    out[,'rYZ'] <- covaA/(sA*sa)
    out    
}
