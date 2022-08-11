#' @title
#' Calculate the weighted mean age
#' 
#' @description
#' Averages heteroscedastic data either using the ordinary weighted
#' mean, or using a random effects model with two sources of variance.
#' Computes the MSWD of a normal fit without
#' overdispersion. Implements a modified Chauvenet criterion to detect
#' and reject outliers. Only propagates the systematic uncertainty
#' associated with decay constants and calibration factors after
#' computing the weighted mean isotopic composition. Does not propagate
#' the uncertainty of any initial daughter correction, because this is
#' neither a purely random or purely systematic uncertainty.
#' 
#' @details
#' Let \eqn{\{t_1, ..., t_n\}} be a set of n age estimates
#' determined on different aliquots of the same sample, and let
#' \eqn{\{s[t_1], ..., s[t_n]\}} be their analytical
#' uncertainties. \code{IsoplotR} then calculates the weighted mean of
#' these data using one of two methods:
#'
#' \enumerate{
#'
#' \item The ordinary error-weighted mean:
#'
#' \eqn{\mu = \sum(t_i/s[t_i]^2)/\sum(1/s[t_i]^2)}
#'
#' \item A random effects model with two sources of variance:
#'
#' \eqn{\log[t_i] \sim N(\log[\mu], \sigma^2 = (s[t_i]/t_i)^2 + \omega^2 )}
#'
#' where \eqn{\mu} is the mean, \eqn{\sigma^2} is the total variance
#' and \eqn{\omega} is the 'overdispersion'. This equation can be
#' solved for \eqn{\mu} and \eqn{\omega} by the method of maximum
#' likelihood.
#' 
#' }
#'
#' IsoplotR uses a modified version of Chauvenet's criterion for
#' outlier detection:
#'
#' \enumerate{
#'
#' \item Compute the error-weighted mean (\eqn{\mu}) of the \eqn{n}
#' age determinations \eqn{t_i} using their analytical uncertainties
#' \eqn{s[t_i]}
#'
#' \item For each \eqn{t_i}, compute the probability \eqn{p_i} that
#' that \eqn{|t-\mu|>|t_i-\mu|} for \eqn{t \sim N(\mu, s[t_i]^2 MSWD)}
#' (ordinary weighted mean) or \eqn{\log[t] \sim
#' N(\log[\mu],s[t_i]^2+\omega^2)} (random effects model)
#'
#' \item Let \eqn{p_j \equiv \min(p_1, ..., p_n)}. If
#' \eqn{p_j<0.05/n}, then reject the j\eqn{^{th}} date, reduce \eqn{n}
#' by one (i.e., \eqn{n \rightarrow n-1}) and repeat steps 1 through 3
#' until the surviving dates pass the third step.  }
#'
#' If the analytical uncertainties are small compared to the scatter
#' between the dates (i.e. if \eqn{\omega \gg s[t]} for all \eqn{i}),
#' then this generalised algorithm reduces to the conventional
#' Chauvenet criterion. If the analytical uncertainties are large and
#' the data do not exhibit any overdispersion, then the heuristic
#' outlier detection method is equivalent to Ludwig (2003)'s `2-sigma'
#' method.
#'
#' The uncertainty budget of the weighted mean does not include the
#' uncertainty of the initial daughter correction (if any). This
#' uncertainty is neither a purely systematic nor a purely random
#' uncertainty and cannot easily be propagated with conventional
#' geochronological data processing algorithms. This caveat is
#' especially pertinent to chronometers whose initial daughter
#' composition is determined by isochron regression. You may note that
#' the uncertainties of the weighted mean are usually much smaller
#' than those of the isochron. In this case the isochron errors are
#' more meaningful, and the weighted mean plot should just be used to
#' inspect the residuals of the data around the isochron.
#'
#' @param x a two column matrix of values (first column) and their
#'     standard errors (second column) OR an object of class
#'     \code{UPb}, \code{PbPb}, \code{ThPb}, \code{ArAr}, \code{KCa},
#'     \code{ReOs}, \code{SmNd}, \code{RbSr}, \code{LuHf}, \code{ThU},
#'     \code{fissiontracks} or \code{UThHe}
#' @param random.effects if \code{TRUE}, computes the weighted mean
#'     using a random effects model with two parameters: the mean and
#'     the dispersion. This is akin to a `model-3' isochron
#'     regression.
#' 
#'     if \code{FALSE}, attributes any excess dispersion to an
#'     underestimation of the analytical uncertainties. This akin to a
#'     `model-1' isochron regression.
#' @param ... optional arguments
#' @seealso \code{\link{central}}
#' 
#' @return Returns a list with the following items:
#'
#' \describe{
#'
#' \item{mean}{a two or three element vector with:
#'
#' \code{t}: the weighted mean. An asterisk is added to the plot title
#' if the initial daughter correction is based on an isochron
#' regression, to mark the circularity of using an isochron to compute
#' a weighted mean.
#'
#' \code{s[t]}: the standard error of the weighted mean, excluding the
#' uncertainty of the initial daughter correction.  This is because
#' this uncertainty is neither purely random nor purely systematic.
#'
#' }
#'
#' \item{disp}{a two-element vector with the (over)dispersion and its
#' standard error.}
#'
#' \item{mswd}{the Mean Square of the Weighted Deviates
#' (a.k.a. `reduced Chi-square' statistic)}
#'
#' \item{df}{the number of degrees of freedom of the Chi-square test
#' for homogeneity (\eqn{df=n-1}, where \eqn{n} is the number of
#' samples).}
#'
#' \item{p.value}{the p-value of a Chi-square test with \eqn{df}
#' degrees of freedom, testing the null hypothesis that the underlying
#' population is not overdispersed.}
#'
#' \item{valid}{vector of logical flags indicating which steps are
#' included into the weighted mean calculation}
#'
#' \item{plotpar}{list of plot parameters for the weighted mean
#' diagram, including \code{mean} (the mean value), \code{ci} (a grey
#' rectangle with the (1 s.e., 2 s.e. or 100[1-\eqn{\alpha}]\%,
#' depending on the value of \code{oerr}) confidence interval ignoring
#' systematic errors), \code{ci.exterr} (a grey rectangle with the
#' confidence interval including systematic errors), \code{dash1} and
#' \code{dash2} (lines marking the confidence interval augmented by
#' \eqn{\sqrt{mswd}} overdispersion if \code{random.effects=FALSE}),
#' and marking the confidence limits of a normal distribution whose
#' standard deviation equals the overdispersion parameter if
#' \code{random.effects=TRUE}). }
#'
#' }
#'
#' @rdname weightedmean
#' @export
weightedmean <- function(x,...){
    UseMethod("weightedmean",x)
}
#' @param detect.outliers logical flag indicating whether outliers
#'     should be detected and rejected using Chauvenet's Criterion.
#' @param plot logical flag indicating whether the function should
#'     produce graphical output or return numerical values to the
#'     user.
#' @param from minimum y-axis limit. Setting \code{from=NA} scales the
#'     plot automatically.
#' @param to maximum y-axis limit. Setting \code{to=NA} scales the
#'     plot automatically.
#' @param levels a vector with additional values to be displayed as
#'     different background colours of the plot symbols.
#' @param clabel label of the colour legend
#' @param rect.col Fill colour for the measurements or age estimates. This can
#'     either be a single colour or multiple colours to form a colour
#'     ramp (to be used if \code{levels!=NA}):
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
#' For empty boxes, set \code{rect.col=NA}
#' 
#' @param outlier.col if \code{detect.outliers=TRUE}, the outliers are
#'     given a different colour.
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output.
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
#' @param ranked plot the aliquots in order of increasing age?
#' @param hide vector with indices of aliquots that should be removed
#'     from the weighted mean plot.
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the weighted mean calculation.
#' @param omit.col colour that should be used for the omitted
#'     aliquots.
#' @importFrom grDevices rgb
#' @rdname weightedmean
#' @export
weightedmean.default <- function(x,from=NA,to=NA,random.effects=FALSE,
                                 detect.outliers=TRUE,plot=TRUE,
                                 levels=NA,clabel="",
                                 rect.col=c("#00FF0080","#FF000080"),
                                 outlier.col="#00FFFF80",sigdig=2,
                                 oerr=3,ranked=FALSE,hide=NULL,
                                 omit=NULL,omit.col=NA,...){
    ns <- nrow(x)
    calcit <- (1:ns)%ni%c(hide,omit)
    X <- x[,1]
    sX <- x[,2]
    valid <- !is.na(X) & !is.na(sX) & calcit
    nvalid <- count(valid)
    if (detect.outliers){
        while (TRUE & nvalid>2){
            valid <- chauvenet(X,sX,valid=valid,
                               random.effects=random.effects)
            if (count(valid) < nvalid) { nvalid <- count(valid) }
            else { break }
        }
    }
    out <- get.weightedmean(X,sX,random.effects=random.effects,
                            valid=valid,oerr=oerr)
    if (plot){
        plot_weightedmean(X,sX,fit=out,from=from,to=to,levels=levels,
                          clabel=clabel,rect.col=rect.col,
                          outlier.col=outlier.col,sigdig=sigdig,
                          oerr=oerr,ranked=ranked,hide=hide,
                          omit=omit,omit.col=omit.col,...)
    }
    invisible(out)
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (\code{type}=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (\code{type}=4), the concordia age (\code{type}=5), or the
#'     \eqn{^{208}}Pb/\eqn{^{232}}Th age (\code{type}=6).
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age is used. This parameter is
#'     only used if \code{type=4}.
#' @param cutoff.disc discordance cutoff filter. This is an object of
#'     class \code{discfilter}
#' @param exterr propagate decay constant uncertainties?
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
#' circularity in the weighted age calculation. In this case the
#' weighted mean plot just serves as a way to visualise the residuals
#' of the data around the isochron, and one should be careful not to
#' over-interpret the numerical output.
#'
#' \code{3}: use the Stacey-Kramers two-stage model to infer the
#' initial Pb-composition (only applicable if \code{x} has class
#' \code{UPb})
#'
#' @examples
#' ages <- c(251.9,251.59,251.47,251.35,251.1,251.04,250.79,250.73,251.22,228.43)
#' errs <- c(0.28,0.28,0.63,0.34,0.28,0.63,0.28,0.4,0.28,0.33)
#' weightedmean(cbind(ages,errs))
#'
#' attach(examples)
#' weightedmean(LudwigMean)
#' @rdname weightedmean
#' @export
weightedmean.UPb <- function(x,random.effects=FALSE,
                             detect.outliers=TRUE,plot=TRUE,from=NA,
                             to=NA,levels=NA,clabel="",
                             rect.col=c("#00FF0080","#FF000080"),
                             outlier.col="#00FFFF80",sigdig=2,type=4,
                             cutoff.76=1100,oerr=3,cutoff.disc=discfilter(),
                             exterr=TRUE,ranked=FALSE,common.Pb=0,
                             hide=NULL,omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        type=type,cutoff.76=cutoff.76,
                        cutoff.disc=cutoff.disc,sigdig=sigdig,
                        oerr=oerr,exterr=exterr,units=' Ma',
                        ranked=ranked,hide=hide,omit=omit,
                        omit.col=omit.col,common.Pb=common.Pb,...)
}
#' @rdname weightedmean
#' @export
weightedmean.PbPb <- function(x,random.effects=FALSE,
                              detect.outliers=TRUE,plot=TRUE, from=NA,
                              to=NA,levels=NA,clabel="",
                              rect.col=c("#00FF0080","#FF000080"),
                              outlier.col="#00FFFF80",sigdig=2,
                              oerr=3,exterr=TRUE,common.Pb=2,
                              ranked=FALSE,hide=NULL,omit=NULL,
                              omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,exterr=exterr,
                        units=' Ma',ranked=ranked,hide=hide,omit=omit,
                        omit.col=omit.col,common.Pb=common.Pb,...)
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
#' in the weighted age calculation. In this case the weighted mean
#' plot just serves as a way to visualise the residuals of the data
#' around the isochron, and one should be careful not to
#' over-interpret the numerical output.
#' 
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
#' @rdname weightedmean
#' @export
weightedmean.ThU <- function(x,random.effects=FALSE,
                             detect.outliers=TRUE,plot=TRUE, from=NA,
                             to=NA,levels=NA,clabel="",
                             rect.col=c("#00FF0080","#FF000080"),
                             outlier.col="#00FFFF80",sigdig=2,
                             oerr=3,ranked=FALSE,Th0i=0,hide=NULL,
                             omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,ranked=ranked,
                        Th0i=Th0i,units=' ka',hide=hide,
                        omit=omit,omit.col=omit.col,...)
}
#' @rdname weightedmean
#' @export
weightedmean.ArAr <- function(x,random.effects=FALSE,
                              detect.outliers=TRUE,plot=TRUE, from=NA,
                              to=NA,levels=NA,clabel="",
                              rect.col=c("#00FF0080","#FF000080"),
                              outlier.col="#00FFFF80",sigdig=2,
                              oerr=3,exterr=TRUE,ranked=FALSE,i2i=FALSE,
                              hide=NULL,omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,exterr=exterr,i2i=i2i,
                        units=' Ma',ranked=ranked,hide=hide,omit=omit,
                        omit.col=omit.col,...)
}
#' @rdname weightedmean
#' @export
weightedmean.KCa <- function(x,random.effects=FALSE,
                             detect.outliers=TRUE,plot=TRUE, from=NA,
                             to=NA,levels=NA,clabel="",
                             rect.col=c("#00FF0080","#FF000080"),
                             outlier.col="#00FFFF80",sigdig=2,
                             oerr=3,exterr=TRUE,ranked=FALSE,i2i=FALSE,
                             hide=NULL,omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,exterr=exterr,
                        i2i=i2i,units=' Ma',ranked=ranked,hide=NULL,
                        omit=omit,omit.col=omit.col,...)
}
#' @rdname weightedmean
#' @export
weightedmean.ThPb <- function(x,random.effects=FALSE,
                              detect.outliers=TRUE,plot=TRUE, from=NA,
                              to=NA,levels=NA,clabel="",
                              rect.col=c("#00FF0080","#FF000080"),
                              outlier.col="#00FFFF80",sigdig=2,
                              oerr=3,exterr=TRUE,ranked=FALSE,i2i=TRUE,
                              hide=NULL,omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig, oerr=oerr,exterr=exterr,
                        i2i=i2i,units=' Ma',ranked=ranked,hide=hide,
                        omit=omit,omit.col=omit.col,...)
}
#' @rdname weightedmean
#' @export
weightedmean.ReOs <- function(x,random.effects=FALSE,
                              detect.outliers=TRUE,plot=TRUE, from=NA,
                              to=NA,levels=NA,clabel="",
                              rect.col=c("#00FF0080","#FF000080"),
                              outlier.col="#00FFFF80",sigdig=2,
                              oerr=3,exterr=TRUE,ranked=FALSE,i2i=TRUE,
                              hide=NULL,omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,exterr=exterr,
                        i2i=i2i,units=' Ma',ranked=ranked,hide=hide,
                        omit=omit,omit.col=omit.col,...)
}
#' @rdname weightedmean
#' @export
weightedmean.SmNd <- function(x,random.effects=FALSE,
                              detect.outliers=TRUE,plot=TRUE,from=NA,
                              to=NA,levels=NA,clabel="",
                              rect.col=c("#00FF0080","#FF000080"),
                              outlier.col="#00FFFF80",sigdig=2,
                              oerr=3,exterr=TRUE,ranked=FALSE,
                              i2i=TRUE,hide=NULL,omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,exterr=exterr,
                        i2i=i2i,units=' Ma',ranked=ranked,hide=hide,
                        omit=omit,omit.col=omit.col,...)
}
#' @rdname weightedmean
#' @export
weightedmean.RbSr <- function(x,random.effects=FALSE,
                              detect.outliers=TRUE,plot=TRUE,from=NA,
                              to=NA,levels=NA,clabel="",
                              rect.col=c("#00FF0080","#FF000080"),
                              outlier.col="#00FFFF80",sigdig=2,
                              oerr=3,exterr=TRUE,i2i=TRUE,ranked=FALSE,
                              hide=NULL,omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,exterr=exterr,
                        i2i=i2i,units=' Ma',ranked=ranked,hide=hide,
                        omit=omit,omit.col=omit.col,...)
}
#' @rdname weightedmean
#' @export
weightedmean.LuHf <- function(x,random.effects=FALSE,
                              detect.outliers=TRUE,plot=TRUE,
                              from=NA,to=NA,levels=NA,clabel="",
                              rect.col=c("#00FF0080","#FF000080"),
                              outlier.col="#00FFFF80",sigdig=2,
                              oerr=3,exterr=TRUE,i2i=TRUE,ranked=FALSE,
                              hide=NULL,omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,exterr=exterr,
                        i2i=i2i,units=' Ma',ranked=ranked,hide=hide,
                        omit=omit,omit.col=omit.col,...)
}
#' @rdname weightedmean
#' @export
weightedmean.UThHe <- function(x,random.effects=FALSE,
                               detect.outliers=TRUE,plot=TRUE,
                               from=NA,to=NA,levels=NA,clabel="",
                               rect.col=c("#00FF0080","#FF000080"),
                               outlier.col="#00FFFF80",sigdig=2,
                               oerr=3,ranked=FALSE,hide=NULL,
                               omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,exterr=FALSE,
                        units=' Ma',ranked=ranked,hide=hide,
                        omit=omit,omit.col=omit.col,...)
}
#' @rdname weightedmean
#' @export
weightedmean.fissiontracks <- function(x,random.effects=FALSE,
                                       detect.outliers=TRUE,plot=TRUE,
                                       from=NA,to=NA,levels=NA,clabel="",
                                       rect.col=c("#00FF0080","#FF000080"),
                                       outlier.col="#00FFFF80",
                                       sigdig=2,oerr=3,exterr=TRUE,ranked=FALSE,
                                       hide=NULL,omit=NULL,omit.col=NA,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,levels=levels,clabel=clabel,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,oerr=oerr,exterr=exterr,
                        units=' Ma',ranked=ranked,hide=hide,
                        omit=omit,omit.col=omit.col,...)
}
weightedmean_helper <- function(x,random.effects=FALSE,
                                detect.outliers=TRUE,plot=TRUE,
                                from=NA,to=NA,levels=NA,clabel="",
                                rect.col=c("#00FF0080","#FF000080"),
                                outlier.col="#00FFFF80",type=4,
                                cutoff.76=1100,cutoff.disc=discfilter(),
                                sigdig=2,oerr=3,exterr=TRUE,ranked=FALSE,
                                i2i=FALSE,common.Pb=1,units='',Th0i=0,
                                hide=NULL,omit=NULL,omit.col=NA,...){
    tt <- get.ages(x,type=type,cutoff.76=cutoff.76,cutoff.disc=cutoff.disc,
                   i2i=i2i,omit4c=unique(c(hide,omit)),
                   common.Pb=common.Pb,Th0i=Th0i)
    fit <- weightedmean.default(tt,random.effects=random.effects,
                                detect.outliers=detect.outliers,
                                oerr=oerr,plot=FALSE,hide=hide,omit=omit)
    if (exterr)
        out <- add.exterr.to.wtdmean(x,fit,oerr=oerr,
                                     cutoff.76=cutoff.76,type=type)
    else out <- fit
    if (plot){
        plot_weightedmean(tt[,1],tt[,2],from=from,to=to,fit=out,
                          levels=levels,clabel=clabel,
                          rect.col=rect.col,outlier.col=outlier.col,
                          sigdig=sigdig,oerr=oerr,units=units,
                          ranked=ranked,hide=hide,omit=omit,omit.col=omit.col,
                          caveat=(i2i|common.Pb==2|Th0i==1),...)
    }
    invisible(out)
}

get.weightedmean <- function(X,sX,random.effects=FALSE,valid=TRUE,oerr=1){
    ns <- length(X)
    x <- X[valid]
    sx <- sX[valid]
    out <- list()
    out$valid <- valid
    if (length(x)<=1){
        out$mean <- c(x,sx)
        out$mswd <- 0
        out$p.value <- 1
        return(out)
    }
    out$random.effects <- random.effects
    out$df <- length(x)-1 # degrees of freedom for the homogeneity test
    if (random.effects){ # random effects model:
        if (all(x>0)){
            fit <- central(cbind(x,sx))
            out$mean <- fit$age
            out$disp <- fit$disp*out$mean['t']
            out$mswd <- fit$mswd
            out$p.value <- fit$p.value
        } else {
            out$mean <- rep(NA,2)
            out$disp <- rep(NA,2)
            names(out$mean) <- c('t','s[t]')
            names(out$disp) <- c('w','s[w]')
            fit <- continuous_mixture(x,sx)
            out$mean['t'] <- fit$mu[1]
            out$mean['s[t]'] <- fit$mu[2]
            out$disp['w'] <- fit$sigma[1]
            out$disp['s[w]'] <- fit$sigma[2]
            SS <- sum(((x-out$mean['t'])/sx)^2)
            out$mswd <- SS/out$df
            out$p.value <- 1-stats::pchisq(SS,out$df)
        }
    } else { # Ludwig's Isoplot approach:
        out$mean <- NULL
        w <- 1/sx^2
        out$mean['t'] <- sum(w*x)/sum(w)
        out$mean['s[t]'] <- 1/sqrt(sum(w))
        SS <- sum(((x-out$mean['t'])/sx)^2)
        out$mswd <- SS/out$df
        out$p.value <- 1-stats::pchisq(SS,out$df)
        if (inflate(out))
            out$mean['disp[t]'] <- sqrt(out$mswd)*out$mean['s[t]']
    }
    plotpar <- list()
    plotpar$mean <- list(x=c(0,ns+1),y=rep(out$mean['t'],2))
    cit <- ci(x=out$mean['t'],sx=out$mean['s[t]'],oerr=oerr,absolute=TRUE)
    plotpar$ci <- list(x=c(0,ns+1,ns+1,0),
                       y=c(rep(out$mean['t']+cit,2),rep(out$mean['t']-cit,2)))
    plotpar$ci.exterr <- NULL # to be defined later
    if (out$random.effects){
        cid <- ci(x=out$mean['t'],sx=out$disp['w'],
                  oerr=oerr,absolute=TRUE)
    } else if (length(out$mean)>2){
        cid <- ci(x=out$mean['t'],sx=out$mean['disp[t]'],
                  oerr=oerr,absolute=TRUE)
    } else {
        cid <- cit
    }
    plotpar$dash1 <- list(x=c(0,ns+1),y=rep(out$mean['t']-cid,2))
    plotpar$dash2 <- list(x=c(0,ns+1),y=rep(out$mean['t']+cid,2))
    out$plotpar <- plotpar
    out
}

plot_weightedmean <- function(X,sX,fit,from=NA,to=NA,levels=NA,clabel="",
                              rect.col=c("#00FF0080","#FF000080"),
                              outlier.col="#00FFFF80",sigdig=2,
                              oerr=3,units='',ranked=FALSE,hide=NULL,
                              omit=NULL,omit.col=NA,caveat=FALSE,...){
    NS <- length(X)
    plotit <- (1:NS)%ni%hide
    calcit <- (1:NS)%ni%c(hide,omit)
    colour <- set.ellipse.colours(ns=NS,levels=levels,col=rect.col,
                                  hide=hide,omit=which(!fit$valid),
                                  omit.col=omit.col)
    Xerr <- ci(X,sX,oerr=oerr,absolute=TRUE)
    x <- X[plotit]
    xerr <- Xerr[plotit]
    valid <- fit$valid[plotit]
    calcit <- calcit[plotit]
    colour <- colour[plotit]
    ns <- length(x)
    if (ranked){
        i <- order(x)
        x <- x[i]
        xerr <- xerr[i]
        valid <- valid[i]
        levels <- levels[i]
        calcit <- calcit[i]
        colour <- colour[i]
    }
    if (is.na(from)) minx <- min(x-xerr,fit$plotpar$dash1$y,na.rm=TRUE)
    else minx <- from
    if (is.na(to)) maxx <- max(x+xerr,fit$plotpar$dash2$y,na.rm=TRUE)
    else maxx <- to
    graphics::plot(c(0,ns+1),c(minx,maxx),type='n',
                   axes=FALSE,xlab='N',ylab='',...)
    if (!is.null(fit$plotpar$ci.exterr))
        graphics::polygon(fit$plotpar$ci.exterr,col='gray90',border=NA)
    graphics::polygon(fit$plotpar$ci,col='gray75',border=NA)
    graphics::lines(fit$plotpar$mean)
    if (fit$random.effects | (fit$p.value<alpha())){
        graphics::lines(fit$plotpar$dash1,lty=3)
        graphics::lines(fit$plotpar$dash2,lty=3)
    }
    graphics::axis(side=1,at=1:ns)
    graphics::axis(side=2)
    for (i in 1:ns){
        if (!calcit[i]){
            col <- omit.col
        } else if (valid[i]){
            col <- colour[i]
        } else {
            col <- outlier.col
        }
        graphics::rect(xleft=i-0.4,ybottom=x[i]-xerr[i],
                       xright=i+0.4,ytop=x[i]+xerr[i],col=col)
    }
    colourbar(z=levels[valid],fill=rect.col,clabel=clabel)
    graphics::title(wtdmean.title(fit,oerr=oerr,sigdig=sigdig,
                                  units=units,caveat=caveat))
}

wtdmean.title <- function(fit,oerr=3,sigdig=2,units='',caveat=FALSE,...){
    ast <- ifelse(caveat,'*','')
    line1 <- maintit(x=fit$mean[1],sx=fit$mean[-1],
                     ntit=ntit.valid(fit$valid),sigdig=sigdig,df=fit$df,
                     oerr=oerr,units=units,prefix=paste0('mean',ast,' ='))
    if (fit$random.effects){
        line2 <- disptit(fit$disp[1],fit$disp[-1],sigdig=sigdig,
                         oerr=oerr,units=units)
    } else {
        line2 <- mswdtit(mswd=fit$mswd,p=fit$p.value,sigdig=sigdig)
    }
    mymtext(line1,line=1,...)
    mymtext(line2,line=0,...)
}

# prune the data if necessary
# X and sX are some measurements and their standard errors
# valid is a vector of logical flags indicating whether the corresponding
# measurements have already been rejected or not
chauvenet <- function(X,sX,valid,random.effects=FALSE){
    if (sum(valid)<2) return(valid)
    fit <- get.weightedmean(X,sX,random.effects=random.effects,valid=valid)
    if (random.effects){
        if (all(X>0,na.rm=TRUE)){
            x <- log(X)
            mu <- log(fit$mean[1])
            sigma <- sqrt((fit$disp[1]/fit$mean[1])^2 + (sX/X)^2)
        } else {
            x <- X
            mu <- fit$mean[1]
            sigma <- sqrt(fit$disp[1]^2 + sX^2)
        }
    } else {
        x <- X
        mu <- fit$mean[1]
        sigma <- sqrt(fit$mean[2]^2 + max(1,fit$mswd)*sX^2)
    }
    misfit <- abs(x-mu)/sigma
    prob <- 2*(1-stats::pnorm(misfit[valid]))
    iworst <- which.max(misfit[valid])
    ivalid <- which(valid)
    minp <- prob[iworst]
    ns <- length(which(valid))
    if (ns*minp < 0.5) {
        valid[ivalid[iworst]] <- FALSE # remove outlier
    } 
    valid
}

add.exterr.to.wtdmean <- function(x,fit,oerr=3,cutoff.76=1100,type=4){
    out <- fit
    out$mean[c('t','s[t]')] <-
        add.exterr(x,
                   tt=fit$mean['t'],
                   st=fit$mean['s[t]'],
                   cutoff.76=cutoff.76,type=type)
    if (inflate(c(fit,model=1+2*fit$random.effects))){
        out$mean['disp[t]'] <-
            add.exterr(x,
                       tt=fit$mean['t'],
                       st=sqrt(fit$mswd)*fit$mean['s[t]'],
                       cutoff.76=cutoff.76,type=type)[2]
    }
    ns <- length(x)
    cit <- ci(x=out$mean['t'],
              sx=out$mean['s[t]'],
              oerr=oerr,absolute=TRUE)
    ci.exterr <- list(x=c(0,ns+1,ns+1,0),
                      y=c(rep(out$mean['t']+cit,2),
                          rep(out$mean['t']-cit,2)))
    out$plotpar$ci.exterr <- ci.exterr
    out
}
