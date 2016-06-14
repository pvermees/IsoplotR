#' Plot continuous data as cumulative age distributions
#'
#' Plot a dataset as Cumulative Age Distributions (CAD), also known as
#' `empirical cumulative distribution function'.
#' 
#' @param x an object of class \code{UPb} or \code{detritals}
#' @param method a string indicating what kind of age should be plotted.
#'
#' If \code{x} has class \code{UPb}, \code{type} could be one of
#' either \code{t.75}, \code{t.68} (default), \code{t.76} or
#' \code{t.conc}
#'
#' @param pch (optional) plot character
#' @param verticals boolean flag indicating if the horizontal lines of
#'     the CAD should be connected by vertical lines
#' @param xlab x-axis label
#' @param colmap an optional string with the name of one of \code{R}'s
#'     built-in colour palettes (e.g., heat.colors, terrain.colors,
#'     topo.colors, cm.colors), which are to be used for plotting data
#'     of class \code{detritals}.
#' @param col colour to give to single sample datasets (i.e. not of
#'     class \code{detritals})
#' @param ... optional arguments to the generic \code{plot} function
#' @examples
#' data(examples)
#' cadplot(examples$DZ)
#' @export
cadplot <- function(x,method=NA,pch=NA,verticals=TRUE,xlab='age [Ma]',
                    colmap='heat.colors',col='black',...){
    X <- age(x)
    ns <- 1
    if (is(x,'detritals')) {        
        ns <- length(x)
        snames <- names(x)
        col <- do.call(colmap,list(ns))
    } else if (is.na(method)){
        if (is(x,'UPb')) method <- 't.68'
        X <- list(X[,method])
    }
    graphics::plot(range(X),c(0,1),type='n',xlab=xlab,ylab='cumulative probability',...)
    for (i in 1:ns){
        graphics::lines(stats::ecdf(X[[i]]),verticals=verticals,
                        pch=pch,col=col[i],...)
    }
    if (is(x,'detritals'))
        graphics::legend("bottomright",legend=snames,lwd=1,col=col)
}
