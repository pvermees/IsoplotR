#' @title
#' Visualise heteroscedastic data on a radial plot
#' 
#' @description
#' Implementation of a graphical device developed by Rex Galbraith to
#' display several estimates of the same quantity that have different
#' standard errors. Serves as a vehicle to display finite and continuous
#' mixture models.
#' 
#' @details
#' The radial plot (Galbraith, 1988, 1990) is a graphical device that
#' was specifically designed to display heteroscedastic data, and is
#' constructed as follows.  Consider a set of dates
#' \eqn{\{t_1,...,t_i,...,t_n\}} and uncertainties
#' \eqn{\{s[t_1],...,s[t_i],...,s[t_n]\}}. Define \eqn{z_i = z[t_i]}
#' to be a transformation of \eqn{t_i} (e.g., \eqn{z_i = log[t_i]}),
#' and let \eqn{s[z_i]} be its propagated analytical uncertainty
#' (i.e., \eqn{s[z_i] = s[t_i]/t_i} in the case of a logarithmic
#' transformation). Create a scatter plot of \eqn{(x_i,y_i)} values,
#' where \eqn{x_i = 1/s[z_i]} and \eqn{y_i = (z_i-z_\circ)/s[z_i]},
#' where \eqn{z_\circ} is some reference value such as the mean. The
#' slope of a line connecting the origin of this scatter plot with any
#' of the \eqn{(x_i,y_i)}s is proportional to \eqn{z_i} and, hence,
#' the date \eqn{t_i}.
#'
#' These dates can be more easily visualised by drawing a radial scale
#' at some convenient distance from the origin and annotating it with
#' labelled ticks at the appropriate angles. While the angular
#' position of each data point represents the date, its horizontal
#' distance from the origin is proportional to the
#' precision. Imprecise measurements plot on the left hand side of the
#' radial plot, whereas precise age determinations are found further
#' towards the right. Thus, radial plots allow the observer to assess
#' both the magnitude and the precision of quantitative data in one
#' glance.
#' 
#' @param x Either an \code{[nx2]} matix of (transformed) values z
#'     and their standard errors s
#'
#' OR
#'
#' and object of class \code{fissiontracks}, \code{UThHe},
#' \code{ArAr}, \code{KCa}, \code{ReOs}, \code{SmNd}, \code{RbSr},
#' \code{LuHf}, \code{ThU}, \code{PbPb}, \code{ThPb} or \code{UPb}
#' @param from minimum age limit of the radial scale
#' @param to maximum age limit of the radial scale
#' @param z0 central value
#' @param transformation one of either \code{log}, \code{linear},
#'     \code{sqrt} or \code{arcsin} (if \code{x} has class
#'     \code{fissiontracks} and \code{fissiontracks$format} \eqn{\neq
#'     1}).
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output.
#' @param show.numbers boolean flag (\code{TRUE} to show grain
#'     numbers)
#' @param pch plot character (default is a filled circle)
#' @param levels a vector with additional values to be displayed as
#'     different background colours of the plot symbols.
#' @param clabel label of the colour legend
#' @param bg Fill colour for the plot symbols. This can either be a
#'     single colour or multiple colours to form a colour ramp (to be
#'     used if \code{levels!=NA}):
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
#' a reversed palette: \code{rev(topo.colors(n=100,alpha=0.5))}, etc.
#'
#' for plot symbols, set \code{bg=NA}
#'
#' @param col text colour to be used if \code{show.numbers=TRUE}
#' @param k number of peaks to fit using the finite mixture models of
#'     Galbraith and Laslett (1993). Setting \code{k='auto'}
#'     automatically selects an optimal number of components based on
#'     the Bayes Information Criterion (BIC). Setting \code{k='min'}
#'     estimates the minimum value using a three parameter model
#'     consisting of a Normal distribution truncated by a discrete
#'     component.
#' @param np number of parameters for the minimum age
#'     model. Must be either 3 or 4.
#' @param markers vector of ages of radial marker lines to add to the
#'     plot.
#' 
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
#'
#' @param units measurement units to be displayed in the legend.
#' @param hide vector with indices of aliquots that should be removed
#'     from the radial plot.
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the central age calculation or mixture models.
#' @param omit.col colour that should be used for the omitted
#'     aliquots.
#' @param ... additional arguments to the generic \code{points}
#'     function
#' @return does not produce any numerical output, but does report the
#'     central age and the results of any mixture modelling in the
#'     title. An asterisk is added to the plot title if the initial
#'     daughter correction is based on an isochron regression, to
#'     highlight the circularity of using an isochron to compute a
#'     central age, and to indicate that the reported uncertainties do
#'     not include the uncertainty of the initial daughter correction.
#'     This is because this uncertainty is neither purely random nor
#'     purely systematic.
#' @seealso \code{\link{peakfit}}, \code{\link{central}}
#' @references Galbraith, R.F., 1988. Graphical display of estimates
#'     having differing standard errors. Technometrics, 30(3),
#'     pp.271-281.
#'
#' Galbraith, R.F., 1990. The radial plot: graphical assessment of
#' spread in ages. International Journal of Radiation Applications and
#' Instrumentation. Part D. Nuclear Tracks and Radiation Measurements,
#' 17(3), pp.207-214.
#'
#' Galbraith, R.F. and Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear Tracks and Radiation
#' Measurements, 21(4), pp.459-470.
#' @examples
#' attach(examples)
#' radialplot(FT1)
#'
#' dev.new()
#' radialplot(LudwigMixture,k='min')
#' @rdname radialplot
#' @export
radialplot <- function(x,...){ UseMethod("radialplot",x) }
#' @rdname radialplot
#' @export
radialplot.default <- function(x,from=NA,to=NA,z0=NA,transformation='log',
                               sigdig=2,show.numbers=FALSE,pch=21,levels=NA,
                               clabel="",bg=c("yellow","red"),col='black',
                               k=0,np=4,markers=NULL,oerr=3,units='',
                               hide=NA,omit=NULL,omit.col=NA,...){
    x <- x[,c(1,2),drop=FALSE]
    ns <- nrow(x)
    calcit <- (1:ns)%ni%c(hide,omit)
    x2calc <- clear(x,hide,omit)
    peaks <- peakfit(x2calc,k=k,np=np,sigdig=sigdig,oerr=oerr)
    markers <- c(markers,peaks$peaks['t',])
    radial_helper(x,from=from,to=to,z0=z0,transformation=transformation,
                  show.numbers=show.numbers,pch=pch,levels=levels,
                  clabel=clabel,bg=bg,col=col,markers=markers,
                  hide=hide,omit=omit,omit.col=omit.col,...)
    fit <- central(x2calc)
    graphics::title(radial.title(fit,sigdig=sigdig,oerr=oerr,
                                 units=units,ntit=get.ntit(x2calc[,1])))
    if (!is.null(peaks$legend))
        graphics::legend('bottomleft',legend=peaks$legend,bty='n')
}
#' @param exterr include the external sources of uncertainty into
#'     the error propagation for the central age and mixture models?
#' @rdname radialplot
#' @export
radialplot.fissiontracks <- function(x,from=NA,to=NA,z0=NA,
                                     transformation='arcsin',sigdig=2,
                                     show.numbers=FALSE,pch=21,
                                     levels=NA,clabel="",bg=c("yellow","red"),
                                     col='black',markers=NULL,k=0,np=4,
                                     exterr=TRUE,oerr=3,hide=NULL,
                                     omit=NULL,omit.col=NA,...){
    ns <- length(x)
    calcit <- (1:ns)%ni%c(hide,omit)
    plotit <- (1:ns)%ni%hide
    x2calc <- clear(x,hide,omit)
    x2plot <- clear(x,hide)
    peaks <- peakfit(x2calc,k=k,np=np,exterr=exterr,sigdig=sigdig,oerr=oerr)
    markers <- c(markers,peaks$peaks['t',])
    X <- x2zs(x2plot,z0=z0,from=from,to=to,transformation=transformation)
    pcol <- set.ellipse.colours(ns=ns,levels=levels[plotit],
                                col=bg,hide=hide,omit=omit,
                                omit.col=omit.col)
    tcol <- rep(col,ns)
    tcol[omit] <- 'grey'
    radial.plot(X,zeta=x$zeta[1],rhoD=x$rhoD[1],
                show.numbers=show.numbers,pch=pch,
                levels=levels[plotit],clabel=clabel,
                markers=markers,bg=pcol[plotit],
                col=tcol[plotit],sn=(1:ns)[plotit],...)
    colourbar(z=levels[calcit],fill=bg,stroke=col,clabel=clabel)
    fit <- central(x2calc,exterr=exterr)
    graphics::title(radial.title(fit,sigdig=sigdig,oerr=oerr,
                                 units=' Ma',ntit=get.ntit(x2calc)))
    if (!is.null(peaks$legend))
        graphics::legend('bottomleft',legend=peaks$legend,bty='n')
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (\code{type}=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (\code{type}=4), the concordia age (\code{type}=5), or the
#'     \eqn{^{208}}Pb/\eqn{^{232}}Th age (\code{type}=6).
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age is used. This parameter is
#'     only used if \code{type=4}.
#' @param cutoff.disc discordance cutoff filter. This is an object of
#'     class \code{\link{discfilter}}.
#' 
#' @param common.Pb common lead correction:
#'
#' \code{0}: none
#'
#' \code{1}: use the Pb-composition stored in
#' 
#' \code{settings('iratio','Pb207Pb206')} (if \code{x} has class
#' \code{UPb} and \code{x$format<4});
#' 
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')} (if \code{x} has class
#' \code{PbPb} or \code{x} has class \code{UPb} and
#' \code{3<x$format<7}); or
#'
#' \code{settings('iratio','Pb208Pb206')} and
#' \code{settings('iratio','Pb208Pb207')} (if \code{x} has class
#' \code{UPb} and \code{x$format=7} or \code{8}).
#'
#' \code{2}: remove the common Pb by projecting the data along an
#' inverse isochron. Note: choosing this option introduces a degree of
#' circularity in the central age calculation. In this case the radial
#' plot just serves as a way to visualise the residuals of the data
#' around the isochron, and one should be careful not to
#' over-interpret the numerical output.
#'
#' \code{3}: use the Stacey-Kramers two-stage model to infer the
#' initial Pb-composition
#' @rdname radialplot
#' @export
radialplot.UPb <- function(x,from=NA,to=NA,z0=NA,
                           transformation='log',type=4,
                           cutoff.76=1100,cutoff.disc=discfilter(),
                           show.numbers=FALSE,pch=21,sigdig=2,
                           levels=NA,clabel="",bg=c("yellow","red"),
                           col='black',markers=NULL,k=0,np=4,
                           exterr=TRUE,common.Pb=0,oerr=3,hide=NULL,
                           omit=NULL,omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,type=type,
               cutoff.76=cutoff.76,cutoff.disc=cutoff.disc,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,
               k=k,np=np,exterr=exterr,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,common.Pb=common.Pb,...)
}
#' @param i2i `isochron to intercept': calculates the initial
#' (aka `inherited', `excess', or `common') \eqn{^{40}}Ar/\eqn{^{36}}Ar,
#' \eqn{^{40}}Ca/\eqn{^{44}}Ca, \eqn{^{207}}Pb/\eqn{^{204}}Pb,
#' \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd,
#' \eqn{^{187}}Os/\eqn{^{188}}Os, \eqn{^{230}}Th/\eqn{^{232}}Th,
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf or \eqn{^{204}}Pb/\eqn{^{208}}Pb
#' ratio from an isochron fit. Setting \code{i2i} to \code{FALSE} uses
#' the default values stored in \code{settings('iratio',...)}.
#' 
#' Note that choosing this option introduces a degree of circularity
#' in the central age calculation. In this case the radial plot plot
#' just serves as a way to visualise the residuals of the data around
#' the isochron, and one should be careful not to over-interpret the
#' numerical output.
#' @rdname radialplot
#' @export
radialplot.PbPb <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",bg=c("yellow","red"),
                            col='black',markers=NULL,k=0,np=4,
                            exterr=TRUE,common.Pb=2,oerr=3,
                            hide=NULL,omit=NULL,omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,k=k,
               np=np,exterr=exterr,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,common.Pb=common.Pb,...)
}
#' @rdname radialplot
#' @export
radialplot.ArAr <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("yellow","red"),col='black',
                            markers=NULL,k=0,np=4,exterr=TRUE,i2i=FALSE,
                            oerr=3,hide=NULL,omit=NULL,omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,k=k,
               np=np,exterr=exterr,i2i=i2i,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,...)
}
#' @rdname radialplot
#' @export
radialplot.KCa <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                           transformation='log',show.numbers=FALSE,
                           pch=21,levels=NA,clabel="",
                           bg=c("yellow","red"),col='black',
                           markers=NULL,k=0,np=4,exterr=TRUE,
                           i2i=FALSE,oerr=3,hide=NULL,omit=NULL,
                           omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,k=k,
               np=np,exterr=exterr,i2i=i2i,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,...)
}
#' @rdname radialplot
#' @export
radialplot.ThPb <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("yellow","red"),col='black',
                            markers=NULL,k=0,np=4,exterr=TRUE,i2i=TRUE,
                            oerr=3,hide=NULL,omit=NULL,omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,k=k,
               np=np,exterr=exterr,i2i=i2i,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,...)
}
#' @rdname radialplot
#' @export
radialplot.UThHe <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                             transformation='log',show.numbers=FALSE,
                             pch=21,levels=NA,clabel="",
                             bg=c("yellow","red"),col='black',
                             markers=NULL,k=0,np=4,oerr=3,
                             hide=NULL,omit=NULL,omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,k=k,
               np=np,exterr=FALSE,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,...)
}
#' @rdname radialplot
#' @export
radialplot.ReOs <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("yellow","red"),col='black',
                            markers=NULL,k=0,np=4,exterr=TRUE,i2i=TRUE,
                            oerr=3,hide=NULL,omit=NULL,omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,k=k,
               np=np,exterr=exterr,i2i=i2i,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,...)
}
#' @rdname radialplot
#' @export
radialplot.SmNd <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("yellow","red"),col='black',
                            markers=NULL,k=0,np=4,exterr=TRUE,i2i=TRUE,
                            oerr=3,hide=NULL,omit=NULL,omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,k=k,
               np=np,exterr=exterr,i2i=i2i,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,...)
}
#' @rdname radialplot
#' @export
radialplot.RbSr <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("yellow","red"),col='black',
                            markers=NULL,k=0,np=4,exterr=TRUE,
                            i2i=TRUE,oerr=3,hide=NULL,omit=NULL,
                            omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,k=k,
               np=np,exterr=exterr,i2i=i2i,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,...)
}
#' @rdname radialplot
#' @export
radialplot.LuHf <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("yellow","red"),col='black',
                            markers=NULL,k=0,np=4,exterr=TRUE,i2i=TRUE,
                            oerr=3,hide=NULL,omit=NULL,omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,
               transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,col=col,markers=markers,k=k,
               np=np,exterr=exterr,i2i=i2i,oerr=oerr,hide=hide,
               omit=omit,omit.col=omit.col,...)
}
#' @param Th0i initial \eqn{^{230}}Th correction.
#'
#' \code{0}: no correction
#'
#' \code{1}: project the data along an isochron fit
#'
#' \code{2}: if \code{x$format} is \code{1} or \code{2}, correct the
#' data using the measured present day \eqn{^{230}}Th/\eqn{^{238}}U,
#' \eqn{^{232}}Th/\eqn{^{238}}U and \eqn{^{234}}U/\eqn{^{238}}U
#' activity ratios in the detritus. If \code{x$format} is \code{3} or
#' \code{4}, correct the data using the measured
#' \eqn{^{238}}U/\eqn{^{232}}Th activity ratio of the whole rock, as
#' stored in \code{x} by the \code{read.data()} function.
#'
#' \code{3}: correct the data using an assumed initial
#' \eqn{^{230}}Th/\eqn{^{232}}Th-ratio for the detritus (only relevant
#' if \code{x$format} is \code{1} or \code{2}).
#' 
#' @rdname radialplot
#' @export
radialplot.ThU <- function(x,from=NA,to=NA,z0=NA,sigdig=2,
                           transformation='log',show.numbers=FALSE,
                           pch=21,levels=NA,clabel="",bg=c("yellow","red"),
                           col='black',markers=NULL,k=0,np=4,Th0i=0,oerr=3,
                           hide=NULL,omit=NULL,omit.col=NA,...){
    age2radial(x,from=from,to=to,z0=z0,sigdig=sigdig,transformation=transformation,
               show.numbers=show.numbers,pch=pch,levels=levels,clabel=clabel,bg=bg,
               col=col,markers=markers,k=k,np=np,exterr=FALSE,Th0i=Th0i,oerr=oerr,
               units=' ka',hide=hide,omit=omit,omit.col=omit.col,...)
}
age2radial <- function(x,from=NA,to=NA,z0=NA,transformation='log',type=4,
                       cutoff.76=1100,cutoff.disc=discfilter(),
                       show.numbers=FALSE,pch=21,levels=NA,sigdig=2,
                       clabel="",bg=c("yellow","red"),col='black',markers=NULL,
                       k=0,np=4,exterr=TRUE,i2i=FALSE,Th0i=0,common.Pb=0,
                       oerr=3,units=' Ma',hide=NULL,omit=NULL,omit.col=NA,...){
    x2calc <- clear(x,hide,omit)
    peaks <- peakfit(x2calc,k=k,np=np,exterr=exterr,sigdig=sigdig,
                     oerr=oerr,i2i=i2i,type=type,cutoff.76=cutoff.76,
                     cutoff.disc=cutoff.disc,common.Pb=common.Pb,Th0i=Th0i)
    markers <- c(markers,peaks$peaks['t',])
    tt <- get.ages(x,type=type,cutoff.76=cutoff.76,cutoff.disc=cutoff.disc,i2i=i2i,
                   omit4c=unique(c(hide,omit)),Th0i=Th0i,common.Pb=common.Pb)
    radial_helper(tt,from=from,to=to,z0=z0,transformation=transformation,
                  show.numbers=show.numbers,pch=pch,levels=levels,
                  clabel=clabel,bg=bg,col=col,markers=markers,
                  hide=hide,omit=omit,omit.col=omit.col,...)
    selection <- clear(1:nrow(tt),unique(c(hide,omit)))
    fit <- central(tt[selection,])
    if (exterr){
        fit$age[1:2] <- add.exterr(x2calc,tt=fit$age[1],st=fit$age[2],
                                   cutoff.76=cutoff.76,type=type)
    }
    graphics::title(radial.title(fit,sigdig=sigdig,oerr=oerr,
                                 units=units,ntit=get.ntit(tt[,1]),
                                 caveat=(i2i|common.Pb==2|Th0i==1)))
    if (!is.null(peaks$legend))
        graphics::legend('bottomleft',legend=peaks$legend,bty='n')
}

radial_helper <- function(x,from=NA,to=NA,z0=NA,transformation='log',
                          show.numbers=FALSE,pch=21,levels=NA,clabel="",
                          bg=c("yellow","red"),col='black',markers=NULL,
                          hide=NULL,omit=NULL,omit.col=NA,...){
    x <- x[,c(1,2),drop=FALSE]
    ns <- nrow(x)
    calcit <- (1:ns)%ni%c(hide,omit)
    plotit <- (1:ns)%ni%hide
    x2calc <- clear(x,hide,omit)
    x2plot <- clear(x,hide)
    X <- x2zs(x2plot,z0=z0,from=from,to=to,transformation=transformation)
    pcol <- set.ellipse.colours(ns=ns,levels=levels,col=bg,hide=hide,
                                omit=omit,omit.col=omit.col)
    tcol <- rep(col,ns)
    tcol[omit] <- 'grey'
    radial.plot(X,show.numbers=show.numbers,pch=pch,
                levels=levels[plotit],clabel=clabel,markers=markers,
                bg=pcol[plotit],col=tcol[plotit],sn=(1:ns)[plotit],...)
    colourbar(z=levels[calcit],fill=bg,stroke=col,clabel=clabel)
}

radial.plot <- function(x,zeta=0,rhoD=0,asprat=3/4,
                        show.numbers=FALSE,levels=NA,clabel="",
                        markers=NULL,pch=21,bg='yellow',col='black',
                        sn=1:length(x$z),...){
    if (show.numbers & all(is.na(levels))) show.points <- FALSE
    else show.points <- TRUE
    exM <- radial.scale(x,zeta,rhoD)
    tticks <- get.radial.tticks(x)
    labelpos <- 4
    if (validLevels(levels)) labelpos <- 2
    plot_radial_lines(tticks,l=0.025,x,exM[1],exM[2],
                      zeta,rhoD,label=TRUE,pos=labelpos)
    if (!is.null(markers)){
        plot_radial_lines(markers,x,exM[1],exM[2],
                          zeta,rhoD,label=FALSE)
    }
    plot_radial_axes(x)
    plot_radial_points(x,show.points=show.points,
                       show.numbers=show.numbers,
                       pch=pch,bg=bg,col=col,sn=sn,...)
}

plot_radial_points <- function(x,show.points=TRUE,show.numbers=FALSE,
                               bg='yellow',pch=21,col='black',
                               sn=1:length(x$z),...){
    rxy <- data2rxry(x)
    rx <- rxy[,1]
    ry <- rxy[,2]
    if (show.numbers & show.points){
        graphics::text(rx,ry,labels=sn,pos=1)
    } else if (show.numbers){
        graphics::text(rx,ry,labels=sn,col=col)
    }
    if (show.points){
        graphics::points(rx,ry,bg=bg,pch=pch,...)
    }

}

plot_radial_axes <- function(x,...){
    xs <- stats::na.omit(x$s)
    graphics::Axis(side=2,at=c(-2,0,2),labels=c(-2,0,2),...)
    if (x$transformation %in% c('arcsin','arctan')){
        plabels <- pretty(c(0,range(1/(2*xs)^2 - 1/2)))
        pticks <- (2*sqrt(plabels+1/2))
    } else {
        plabels <- pretty(c(0,1/xs))
        pticks <- plabels
    }
    graphics::Axis(side=1,at=pticks,labels=plabels,...)
}

data2rxry <- function(x){
    rx <- 1/x$s
    ry <- (x$z-x$z0)/x$s
    cbind(rx,ry)
}

radial.scale <- function(x,zeta=0,rhoD=0,...){
    zm <- t2z(x$from,x,zeta,rhoD)
    zM <- t2z(x$to,x,zeta,rhoD)
    padding <- 1.1
    N <- 50
    a <- grDevices::dev.size()[2]/grDevices::dev.size()[1] # aspect ratio
    e <- a/(zM-zm) # ellipticity of the arc
    # get rxM
    theta <- atan(e*(x$z-x$z0))
    rxy <- data2rxry(x)
    xM <- padding * sqrt(max(
          (rxy[,1]^2+rxy[,2]^2)/(cos(theta)^2+(sin(theta)/e)^2), na.rm=TRUE
          ))
    # plot arc
    Z <- seq(zm-x$z0,zM-x$z0,length.out=N)
    rxy <- z2rxy(Z,e,xM)
    graphics::plot(rxy[,1],rxy[,2],type='l',xlim=c(0,xM),
                   axes=FALSE,bty='n',xlab=x$xlab,
                   ylab='standardised estimate',...)
    c(e,xM)
}

plot_radial_lines <- function(tt,x,e,xM,zeta=0,rhoD=0,l=1,
                              label=FALSE,pos=4,...){
    z <- t2z(tt,x,zeta,rhoD)
    rxyb <- z2rxy(z-x$z0,e,xM)
    rxye <- z2rxy(z-x$z0,e,(1-l)*xM)
    for (i in 1:length(tt)){
        graphics::lines(c(rxyb[i,1],rxye[i,1]),
                        c(rxyb[i,2],rxye[i,2]),...)
        if (label) {
            if (pos==2)
                graphics::text(rxye[i,1],rxye[i,2],
                               labels=tt[i],pos=pos,xpd=NA)
            else
                graphics::text(rxyb[i,1],rxyb[i,2],
                               labels=tt[i],pos=4,xpd=NA)
        }
    }
}

z2rxy <- function(Z,e,xM){
    theta <- atan(e*Z)
    rx <- xM*cos(theta)
    ry <- (xM/e)*sin(theta)
    cbind(rx,ry)
}

t2z <- function(tt,x,zeta,rhoD){
    if (identical(x$transformation,'log')){
        out <- log(tt+x$offset)
    } else if (identical(x$transformation,'arcsin')){
        out <- att(tt,zeta,rhoD)
    } else if (identical(x$transformation,'arctan')){
        out <- atan(sqrt(tt))
    } else if (identical(x$transformation,'linear')){
        out <- tt
    } else if (identical(x$transformation,'sqrt')){
        out <- sqrt(tt)
    }
    out
}

get.radial.tticks <- function(x){
    if (identical(x$transformation,'linear')){
        out <- pretty(c(x$from,x$to))
    } else if (identical(x$transformation,'log')){
        logrange <- log10(c(x$from,x$to)+x$offset)
        out <- grDevices::axisTicks(usr=logrange,log=TRUE)-x$offset
    } else if (x$transformation %in% c('arcsin','arctan')){
        logrange <- log10(c(x$from,x$to))
        out <- grDevices::axisTicks(usr=logrange,log=TRUE)
    }  else if (identical(x$transformation,'sqrt')){
        out <- pretty(sqrt(c(x$from,x$to)))^2
    }
    nt <- length(out)
    reldiff <- (out[1]-x$from)/(out[nt]-out[1])
    sigdig <- ceiling(1-log10(abs(1-x$from/out[1])))
    firsttick <- signif(x$from,sigdig)
    if (out[1] < x$from) {
        out[1] <- firsttick
    } else if (reldiff > 0.2) {
        out <- c(firsttick,out)
        nt <- nt+1
    }
    reldiff <- (x$to-out[nt])/(out[nt]-out[1])
    sigdig <- ceiling(1-log10(abs(1-out[nt]/x$to)))
    lasttick <- signif(x$to,sigdig)
    if (out[nt] > x$to) {
        out[nt] <- lasttick
    } else if (reldiff > 0.2) {
        out <- c(out,lasttick)
    }
    out
}

x2zs <- function(x,...){ UseMethod("x2zs",x) }
x2zs.default <- function(x,z0=NA,from=NA,to=NA,transformation=NA,...){
    out <- list()
    if (is.na(transformation)) out$transformation <- 'log'
    else out$transformation <- transformation
    if (identical(transformation,'log')){
        out$offset <- get.offset(x[,1],from)
        out$z <- log(x[,1]+out$offset)
        out$s <- x[,2]/(x[,1]+out$offset)
        if (out$offset>0){
            out$xlab <- substitute('t/('~sigma~'+'~a~')',list(a=out$offset))
        } else {
            out$xlab <- expression(t/sigma)
        }
    } else if (identical(transformation,'sqrt')){
        out$z <- sqrt(x[,1])
        out$s <- 0.5*x[,2]/out$z
        out$xlab <- 'precision'
    } else  if (identical(transformation,'linear')){
        out$z <- x[,1]
        out$s <- x[,2]
        out$xlab <- expression(1/sigma)
    }
    out$z0 <- get.z0(out,z0,from,to)
    # reset limits if necessary
    if (is.na(from)){
        min.z <- get.min.z(out)
        if (identical(transformation,'log'))
            out$from <- exp(min.z)-out$offset
        else if (identical(transformation,'sqrt'))
            out$from <- min.z^2
        else if (identical(transformation,'linear'))
            out$from <- min.z
    } else {
        out$from <- from
    }
    if (is.na(to)){
        max.z <- get.max.z(out)
        if (identical(transformation,'log'))
            out$to <- exp(max.z)-out$offset
        else if (identical(transformation,'linear'))
            out$to <- max.z
        else if (identical(transformation,'sqrt'))
            out$to <- max.z^2
    } else {
        out$to <- to
    }
    out
}
x2zs.fissiontracks <- function(x,z0=NA,from=NA,to=NA,transformation=NA,...){
    out <- list()
    if (transformation %in% c('linear','log','arcsin','sqrt')){
        out$transformation <- transformation
    } else {
        if (x$format==1) out$transformation <- 'arcsin'
        else out$transformation <- 'log'
    }
    if (x$format==1){
        Ns <- x$x[,'Ns']
        Ni <- x$x[,'Ni']
        if (identical(transformation,'linear')){
            tt <- fissiontrack.age(x,exterr=FALSE)
            out$z <- tt[,'t']
            out$s <- tt[,'s[t]']
            out$z0 <- get.z0(out,z0,from,to)
            out$xlab <- expression(1/sigma)
        }
        if (identical(transformation,'sqrt')){
            tt <- fissiontrack.age(x,exterr=FALSE)
            out$z <- sqrt(tt[,'t'])
            out$s <- 0.5*tt[,'s[t]']/out$z
            out$z0 <- get.z0(out,z0,from,to)
            out$xlab <- 'precision'
        }
        if (identical(transformation,'log')){
            tt <- fissiontrack.age(x,exterr=FALSE)
            if (any(tt[,'t']<=0)) {
                out$transformation <- 'arcsin'
            } else {
                out$offset <- 0
                out$z <- log(tt[,'t'])
                out$s <- tt[,'s[t]']/tt[,'t']
                out$z0 <- get.z0(out,z0,from,to)
                out$xlab <- expression(t/sigma)
            }
        }
        if (identical(out$transformation,'arcsin')){
            out$z <- atan(sqrt((Ns+3/8)/(Ni+3/8)))
            out$s <- sqrt(1/(Ns+Ni+1/2))/2
            if (is.na(z0))
                if (is.na(from) | is.na(to)) {
                    out$z0 <- atan(sqrt(sum(Ns,na.rm=TRUE)/
                                        sum(Ni,na.rm=TRUE)))
                } else {
                    zmin <- att(from,x$zeta[1],x$rhoD[1])
                    zmax <- att(to,x$zeta[1],x$rhoD[1])
                    out$z0 <- mean(c(zmin,zmax),na.rm=TRUE)
                }
            else {
                out$z0 <- att(z0,x$zeta[1],x$rhoD[1])
            }
            out$xlab <- 'Ns+Ni'
        }
        # reset limits if necessary
        if (is.na(from)){
            if (identical(transformation,'log'))
                out$from <- exp(min(out$z,na.rm=TRUE))
            else if (identical(transformation,'arcsin'))
                out$from <- iatt(min(out$z,na.rm=TRUE),
                                 x$zeta[1],x$rhoD[1])
            else if (identical(transformation,'linear'))
                out$from <- min(out$z,na.rm=TRUE)
            else if (identical(transformation,'sqrt'))
                out$from <- min(out$z,na.rm=TRUE)^2
        } else {
            out$from <- from
        }
        if (is.na(to)){
            if (identical(transformation,'log'))
                out$to <- exp(max(out$z,na.rm=TRUE))
            else if (identical(transformation,'arcsin'))
                out$to <- iatt(max(out$z,na.rm=TRUE),
                               x$zeta[1],x$rhoD[1])
            else if (identical(transformation,'linear'))
                out$to <- max(out$z,na.rm=TRUE)
            else if (identical(transformation,'sqrt'))
                out$to <- max(out$z,na.rm=TRUE)^2
        } else {
            out$to <- to
        }
    } else {
        tt <- fissiontrack.age(x,exterr=FALSE)
        if (identical(transformation,'arcsin')) transformation <- 'log'
        out <- x2zs.default(tt,z0=z0,from=from,to=to,
                            transformation=transformation)
    }
    out
}

# only for log and lin transformation
# x is a list containing the items z, s, transformation
# and (if transformation=='log) offset
get.z0 <- function(x,z0=NA,from=NA,to=NA){
    if (is.na(z0)){
        if (is.na(from) | is.na(to)) {
            z0 <- mean(x$z,na.rm=TRUE)
        } else if (identical(x$transformation,'log')) {
            z0 <- mean(log(c(from,to)+x$offset),na.rm=TRUE)
        } else if (identical(x$transformation,'linear')) {
            z0 <- mean(c(from,to),na.rm=TRUE)
        } else if (identical(x$transformation,'sqrt')) {
            z0 <- mean(sqrt(c(from,to)),na.rm=TRUE)
        } else {
            stop('illegal input')
        }
    } else if (identical(x$transformation,'log')){
        z0 <- log(z0+x$offset)
    } else if (identical(x$transformation,'linear')){
        z0 <- z0
    } else if (identical(x$transformation,'sqrt')){
        z0 <- sqrt(z0)
    } else {
        stop('illegal input')
    }
    z0
}

# arctan transformation
att <- function(tt,zeta,rhoD){
    L8 <- lambda('U238')[1]
    atan(sqrt((exp(L8*tt)-1)/(L8*(zeta/1e6)*rhoD/2)))
}
# inverse arctan transformation
iatt <- function(z,zeta,rhoD){
    L8 <- lambda('U238')[1]
    log(1+L8*(zeta/2e6)*rhoD*tan(z)^2)/L8
}

radial.title <- function(fit,sigdig=2,oerr=3,units='',ntit,caveat=FALSE,...){
    ast <- ifelse(caveat,'*','')
    line1 <- maintit(fit$age[1],fit$age[2],ntit=ntit,
                     units=units,sigdig=sigdig,oerr=oerr,
                     prefix=paste0('central age',ast,' ='))
    line2 <- mswdtit(mswd=fit$mswd,p=fit$p.value,sigdig=sigdig)
    if (fit$p.value<alpha()){
        line3 <- disptit(w=fit$disp[1],sw=fit$disp[2],sigdig=sigdig,oerr=oerr)
        mymtext(line3,line=0,...)
        line1line <- 2
        line2line <- 1
    } else {
        line1line <- 1
        line2line <- 0
    }
    mymtext(line1,line=line1line,...)
    mymtext(line2,line=line2line,...)    
}

get.offset <- function(x,from=NA){
    m <- min(c(x,from),na.rm=TRUE)
    if (m>0){
        offset <- 0;
    } else if (m==0){
        offset <- 1;
    } else {
        offset = 10^(floor(log10(-m))+1);
    }
    offset
}

get.min.z <- function(x){
    rxry <- data2rxry(x)
    min.ry <- min(rxry[,2],na.rm=TRUE)
    if (min.ry > -2){ # if the data are underdispersed
        max.rx <- max(rxry[,1],na.rm=TRUE)
        min.z <- x$z0-2/max.rx
    } else {
        min.z <- min(x$z,na.rm=TRUE)
    }
    min.z
}

get.max.z <- function(x){
    rxry <- data2rxry(x)
    max.ry <- max(rxry[,2],na.rm=TRUE)
    if (max.ry < 2){ # if the data are underdispersed
        max.rx <- max(rxry[,1],na.rm=TRUE)
        max.z <- x$z0+2/max.rx
    } else {
        max.z <- max(x$z,na.rm=TRUE)
    }
    max.z
}
