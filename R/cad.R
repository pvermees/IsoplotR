#' Plot continuous data as cumulative age distributions
#'
#' Plot a dataset as a Cumulative Age Distribution (CAD), also known
#' as a `empirical cumulative distribution function'.
#'
#' @details
#' Empirical cumulative distribution functions or cumulative age
#' distributions CADs are the most straightforward way to visualise
#' the probability distribution of multiple dates.  Suppose that we
#' have a set of \eqn{n} dates \eqn{t_i}. The the CAD is a step
#' function that sets out the rank order of the dates against their
#' numerical value:
#'
#' \eqn{CAD(t) = \sum_i 1(t<t_i)/n}
#'
#' where 1(\eqn{\ast}) = 1 if \eqn{\ast} is true and 1(\eqn{\ast}) = 0
#' if \eqn{\ast} is false. CADs have two desirable properties
#' (Vermeesch, 2007). First, they do not require any pre-treatment or
#' smoothing of the data. This is not the case for histograms or
#' kernel density estimates. Second, it is easy to superimpose several
#' CADs on the same plot. This facilitates the intercomparison of
#' multiple samples. The interpretation of CADs is straightforward but
#' not very intuitive. The prominence of individual age components is
#' proportional to the steepness of the CAD. This is different from
#' probability density estimates such as histograms, in which such
#' components stand out as peaks.
#'
#' @param x a numerical vector OR an object of class \code{UPb},
#'     \code{PbPb}, \code{ArAr}, \code{KCa}, \code{UThHe},
#'     \code{fissiontracks}, \code{ReOs}, \code{RbSr}, \code{SmNd},
#'     \code{LuHf}, \code{ThU} or \code{detritals}
#' @seealso \code{\link{kde}}, \code{\link{radialplot}}
#' @rdname cad
#' @export
cad <- function(x,...){ UseMethod("cad",x) }
#' @param pch plot character to mark the beginning of each CAD step
#' @param verticals logical flag indicating if the horizontal lines of
#'     the CAD should be connected by vertical lines
#' @param xlab x-axis label
#' @param colmap an optional string with the name of one of \code{R}'s
#'     built-in colour palettes (e.g., \code{heat.colors},
#'     \code{terrain.colors}, \code{topo.colors}, \code{cm.colors}),
#'     which are to be used for plotting data of class
#'     \code{detritals}.
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param ... optional arguments to the generic \code{plot} function
#' @examples
#' data(examples)
#' cad(examples$DZ,verticals=FALSE,pch=20)
#' @rdname cad
#' @export
cad.default <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                        colmap='heat.colors',col='black',
                        hide=NULL,...){
    calcit <- (1:length(x))%ni%hide
    d <- x[calcit]
    graphics::plot(range(d,na.rm=TRUE),c(0,1),type='n',
                   ylab='cumulative probability',xlab=xlab,...)
    graphics::lines(stats::ecdf(d),pch=pch,verticals=verticals,col=col)
    mymtext(paste0('n=',length(d)),line=0,adj=1)
}
#' @param col colour to give to single sample datasets (not applicable
#'     if \code{x} has class \code{detritals})
#' @rdname cad
#' @export
cad.detritals <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                          colmap='heat.colors',hide=NULL,...){
    if (is.character(hide)) hide <- which(names(x)%in%hide)
    x2plot <- clear(x,hide)
    ns <- length(x2plot)
    col <- do.call(colmap,list(ns))
    graphics::plot(range(x2plot,na.rm=TRUE),c(0,1),type='n',xlab=xlab,
                   ylab='cumulative probability',...)
    for (i in 1:ns){
        graphics::lines(stats::ecdf(x2plot[[i]]),pch=pch,
                        verticals=verticals,col=col[i])
    }
    graphics::legend("bottomright",legend=names(x2plot),lwd=1,col=col)
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (\code{type}=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (\code{type}=4), or the (Wetherill) concordia age
#'     (\code{type}=5)
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U-age and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-age is used. This parameter is
#'     only used if \code{type=4}.
#' @param cutoff.disc two element vector with the maximum and minimum
#'     percentage discordance allowed between the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U and \eqn{^{206}}Pb/\eqn{^{238}}U
#'     age (if \eqn{^{206}}Pb/\eqn{^{238}}U < cutoff.76) or between
#'     the \eqn{^{206}}Pb/\eqn{^{238}}U and
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (if
#'     \eqn{^{206}}Pb/\eqn{^{238}}U > cutoff.76).  Set
#'     \code{cutoff.disc=NA} if you do not want to use this filter.
#' @param common.Pb apply a common lead correction using one of three
#'     methods:
#'
#' \code{1}: use the isochron intercept as the initial Pb-composition
#'
#' \code{2}: use the Stacey-Kramer two-stage model to infer the initial
#' Pb-composition
#'
#' \code{3}: use the Pb-composition stored in
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')}
#'
#' @rdname cad
#' @export
cad.UPb <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                    col='black',type=4,cutoff.76=1100,
                    cutoff.disc=c(-15,5),common.Pb=0,
                    hide=NULL,...){
    tt <- filter.UPb.ages(x,type,cutoff.76,cutoff.disc,common.Pb=common.Pb)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @param i2i `isochron to intercept': calculates the initial (aka
#'     `inherited', `excess', or `common')
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar, \eqn{^{40}}Ca/\eqn{^{44}}Ca,
#'     \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#'     \eqn{^{143}}Nd/\eqn{^{144}}Nd, \eqn{^{187}}Os/\eqn{^{188}}Os,
#'     \eqn{^{230}}Th/\eqn{^{232}}Th or \eqn{^{176}}Hf/\eqn{^{177}}Hf
#'     ratio from an isochron fit. Setting \code{i2i} to \code{FALSE}
#'     uses the default values stored in \code{settings('iratio',...)}
#'     or zero (for the Pb-Pb method). When applied to data of class
#'     \code{ThU}, setting \code{i2i} to \code{TRUE} applies a
#'     detrital Th-correction.
#' 
#' @references Vermeesch, P., 2007. Quantitative geomorphology of the
#'     White Mountains (California) using detrital apatite fission
#'     track thermochronology. Journal of Geophysical Research: Earth
#'     Surface, 112(F3).
#' @rdname cad
#' @export
cad.PbPb <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                     col='black',common.Pb=1,hide=NULL,...){
    tt <- PbPb.age(x,common.Pb=common.Pb)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @rdname cad
#' @export
cad.ArAr <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                     col='black',i2i=FALSE,hide=NULL,...){
    tt <- ArAr.age(x,i2i=i2i)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @rdname cad
#' @export
cad.KCa <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                    col='black',i2i=FALSE,hide=NULL,...){
    tt <- KCa.age(x,i2i=i2i)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @param detritus detrital \eqn{^{230}}Th correction (only applicable
#'     when \code{x$format == 1} or \code{2}).
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
#' @rdname cad
#' @export
cad.ThU <- function(x,pch=NA,verticals=TRUE, xlab='age [ka]',
                    col='black',i2i=FALSE,detritus=0,hide=NULL,...){
    tt <- ThU.age(x,i2i=i2i,detritus=detritus)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @rdname cad
#' @export
cad.ReOs <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                     col='black',i2i=TRUE,hide=NULL,...){
    tt <- ReOs.age(x,i2i=i2i)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @rdname cad
#' @export
cad.SmNd <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                     col='black',i2i=TRUE,hide=NULL,...){
    tt <- SmNd.age(x,i2i=i2i)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @rdname cad
#' @export
cad.RbSr <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                     col='black',i2i=TRUE,hide=NULL,...){
    tt <- RbSr.age(x,i2i=i2i)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @rdname cad
#' @export
cad.LuHf <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                     col='black',i2i=TRUE,hide=NULL,...){
    tt <- LuHf.age(x,i2i=i2i)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @rdname cad
#' @export
cad.UThHe <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                      col='black',hide=NULL,...){
    tt <- UThHe.age(x)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
#' @rdname cad
#' @export
cad.fissiontracks <- function(x,pch=NA,verticals=TRUE,xlab='age [Ma]',
                              col='black',hide=NULL,...){
    tt <- fissiontrack.age(x)[,1]
    cad.default(tt,pch=pch,verticals=verticals,xlab=xlab,
                col=col,hide=hide,...)
}
