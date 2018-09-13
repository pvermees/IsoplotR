#' Calculate the weighted mean age
#'
#' Models the data as a Normal distribution with two sources of
#' variance.  Estimates the mean and `overdispersion' using the method
#' of Maximum Likelihood. Computes the MSWD of a Normal fit without
#' overdispersion. Implements a modified Chauvenet Criterion to detect
#' and reject outliers. Only propagates the analytical uncertainty
#' associated with decay constants and \eqn{\zeta} and J-factors after
#' computing the weighted mean isotopic composition.
#'
#' @details
#' Let \eqn{\{t_1, ..., t_n\}} be a set of n age estimates
#' determined on different aliquots of the same sample, and let
#' \eqn{\{s[t_1], ..., s[t_n]\}} be their analytical
#' uncertainties. \code{IsoplotR} then calculates the weighted mean of
#' these data assuming a Normal distribution with two sources of
#' variance:
#'
#' \eqn{t_i \sim N(\mu, \sigma^2 = s[t_i]^2 + \omega^2 )}
#'
#' where \eqn{\mu} is the mean, \eqn{\sigma^2} is the total variance
#' and \eqn{\omega} is the 'overdispersion'. This equation can be
#' solved for \eqn{\mu} and \eqn{\omega} by the method of maximum
#' likelihood. IsoplotR uses a modified version of Chauvenet's
#' criterion for outlier detection:
#'
#' \enumerate{
#'
#' \item Compute the error-weighted mean (\eqn{\mu}) of the \eqn{n}
#' age determinations \eqn{t_i} using their analytical uncertainties
#' \eqn{s[t_i]}
#'
#' \item For each \eqn{t_i}, compute the probability \eqn{p_i} that
#' that \eqn{|t-\mu|>|t_i-\mu|} for \eqn{t \sim
#' N(0,\sqrt{s[t_i]^2+\omega^2) }}
#'
#' \item Let \eqn{p_j \equiv \min(p_1, ..., p_n)}. If
#' \eqn{p_j<0.05/n}, then reject the j\eqn{^{th}} date, reduce \eqn{n}
#' by one (i.e., \eqn{n \rightarrow n-1}) and repeat steps 1 through 3
#' until the surviving dates pass the third step.  }
#'
#' If the analtyical uncertainties are small compared to the scatter
#' between the dates (i.e. if \eqn{\omega \gg s[t]} for all \eqn{i}),
#' then this generalised algorithm reduces to the conventional
#' Chauvenet criterion. If the analytical uncertainties are large and
#' the data do not exhibit any overdispersion, then the heuristic
#' outlier detection method is equivalent to Ludwig (2003)'s `2-sigma'
#' method.
#'
#' @param x a two column matrix of values (first column) and their
#'     standard errors (second column) OR an object of class
#'     \code{UPb}, \code{PbPb}, \code{ArAr}, \code{KCa}, \code{ReOs},
#'     \code{SmNd}, \code{RbSr}, \code{LuHf}, \code{ThU},
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
#' @return Returns a list with the following items:
#'
#' \describe{
#'
#' \item{mean}{a three element vector with:
#'
#' \code{x}: the weighted mean
#'
#' \code{s[x]}: the standard error of the weighted mean
#'
#' \code{ci[x]}: the \eqn{100(1-\alpha)\%} confidence interval for
#' \code{x}
#'
#' }
#'
#' \item{disp}{a three-element vector with the (over)dispersion and
#' the lower and upper half-widths of its \eqn{100(1-\alpha)\%}
#' confidence interval.}
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
#' diagram}
#' }
#' @rdname weightedmean
#' @export
weightedmean <- function(x,...){ UseMethod("weightedmean",x) }
#' @param detect.outliers logical flag indicating whether outliers
#'     should be detected and rejected using Chauvenet's Criterion.
#' @param plot logical flag indicating whether the function should
#'     produce graphical output or return numerical values to the
#'     user.
#' @param from minimum y-axis limit. Setting \code{from=NA} scales the
#'     plot automatically.
#' @param to maximum y-axis limit. Setting \code{to=NA} scales the
#'     plot automatically.
#' @param rect.col the fill colour of the rectangles used to show the
#'     measurements or age estimates.
#' @param outlier.col if \code{detect.outliers=TRUE}, the outliers are
#'     given a different colour.
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output.
#' @param alpha the confidence limits of the error bars/rectangles.
#' @param ranked plot the aliquots in order of increasing age?
#' @importFrom grDevices rgb
#' @rdname weightedmean
#' @export
weightedmean.default <- function(x,from=NA,to=NA,random.effects=TRUE,
                                 detect.outliers=TRUE,plot=TRUE,
                                 rect.col=rgb(0,1,0,0.5),
                                 outlier.col=rgb(0,1,1,0.5),sigdig=2,
                                 alpha=0.05,ranked=FALSE,...){
    X <- x[,1]
    sX <- x[,2]
    valid <- !is.na(X) & !is.na(sX)
    nvalid <- count(valid)
    if (detect.outliers){
        while (TRUE){
            valid <- chauvenet(X,sX,valid=valid,
                               random.effects=random.effects)
            if (count(valid) < nvalid)
                nvalid <- count(valid)
            else
                break
        }
    }
    out <- get.weightedmean(X,sX,random.effects=random.effects,
                            valid=valid,alpha=alpha)
    if (plot){
        plot_weightedmean(X,sX,fit=out,from=from,to=to,
                          rect.col=rect.col,outlier.col=outlier.col,
                          sigdig=sigdig,alpha=alpha,ranked=ranked,...)
    }
    invisible(out)
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (\code{type}=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (\code{type}=4), or the (Wetherill) concordia age (\code{type}=5)
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age is used. This parameter is
#'     only used if \code{type=4}.
#' @param cutoff.disc two element vector with the maximum and minimum
#'     percentage discordance allowed between the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U and \eqn{^{206}}Pb/\eqn{^{238}}U
#'     age (if \eqn{^{206}}Pb/\eqn{^{238}}U < \code{cutoff.76}) or
#'     between the \eqn{^{206}}Pb/\eqn{^{238}}U and
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (if
#'     \eqn{^{206}}Pb/\eqn{^{238}}U > \code{cutoff.76}).  Set
#'     \code{cutoff.disc=NA} if you do not want to use this filter.
#' @param exterr propagate decay constant uncertainties?
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
#' @examples
#' ages <- c(251.9,251.59,251.47,251.35,251.1,251.04,250.79,250.73,251.22,228.43)
#' errs <- c(0.28,0.28,0.63,0.34,0.28,0.63,0.28,0.4,0.28,0.33)
#' weightedmean(cbind(ages,errs))
#'
#' data(examples)
#' weightedmean(examples$LudwigMean)
#' @rdname weightedmean
#' @export
weightedmean.UPb <- function(x,random.effects=TRUE,
                             detect.outliers=TRUE,plot=TRUE,
                             from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                             outlier.col=rgb(0,1,1,0.5),sigdig=2,
                             type=4,cutoff.76=1100,
                             cutoff.disc=c(-15,5),alpha=0.05,
                             exterr=TRUE,ranked=FALSE,common.Pb=0,...){
    if (common.Pb %in% c(1,2,3))
        X <- common.Pb.correction(x,option=common.Pb)
    else
        X <- x
    weightedmean_helper(X,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,rect.col=rect.col,
                        outlier.col=outlier.col,type=type,
                        cutoff.76=cutoff.76,cutoff.disc=cutoff.disc,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        units='Ma',ranked=ranked,...)
}
#' @rdname weightedmean
#' @export
weightedmean.PbPb <- function(x,random.effects=TRUE,
                              detect.outliers=TRUE,plot=TRUE,
                              from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5),sigdig=2,
                              alpha=0.05,exterr=TRUE,common.Pb=1,
                              ranked=FALSE,...){
    if (common.Pb %in% c(1,2,3))
        X <- common.Pb.correction(x,option=common.Pb)
    else
        X <- x
    weightedmean_helper(X,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,rect.col=rect.col,
                        outlier.col=outlier.col,sigdig=sigdig,
                        alpha=alpha,exterr=exterr,
                        units='Ma',ranked=ranked,...)
}
#' @param i2i `isochron to intercept': calculates the initial (aka
#'     `inherited', `excess', or `common')
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar, \eqn{^{40}}Ca/\eqn{^{44}}Ca,
#'     \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#'     \eqn{^{143}}Nd/\eqn{^{144}}Nd, \eqn{^{187}}Os/\eqn{^{188}}Os,
#'     \eqn{^{230}}Th/\eqn{^{232}}Th or \eqn{^{176}}Hf/\eqn{^{177}}Hf
#'     ratio from an isochron fit. Setting \code{i2i} to \code{FALSE}
#'     uses the default values stored in
#'     \code{settings('iratio',...)}.
#' @param detritus detrital \eqn{^{230}}Th correction (only applicable
#'     when \code{x$format == 1} or \code{2}.
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
#' @param Th02 2-element vector with the assumed initial
#'     \eqn{^{230}}Th/\eqn{^{232}}Th-ratio of the detritus and its
#'     standard error. Only used if \code{detritus==2}
#' @param Th02U48 9-element vector with the measured composition of
#'     the detritus, containing \code{X=0/8}, \code{sX}, \code{Y=2/8},
#'     \code{sY}, \code{Z=4/8}, \code{sZ}, \code{rXY}, \code{rXZ},
#'     \code{rYZ}. Only used if \code{isochron==FALSE} and
#'     \code{detritus==3}
#' @rdname weightedmean
#' @export
weightedmean.ThU <- function(x,random.effects=TRUE,
                             detect.outliers=TRUE,plot=TRUE,
                             from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                             outlier.col=rgb(0,1,1,0.5),sigdig=2,
                             alpha=0.05,ranked=FALSE,i2i=TRUE,
                             detritus=0,Th02=c(0,0),
                             Th02U48=c(0,0,1e6,0,0,0,0,0,0),...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,rect.col=rect.col,
                        outlier.col=outlier.col,sigdig=sigdig,alpha=alpha,
                        ranked=ranked,i2i=i2i,units='ka',detritus=detritus,
                        Th02=Th02,Th02U48=Th02U48,...)
}
#' @rdname weightedmean
#' @export
weightedmean.ArAr <- function(x,random.effects=TRUE,
                              detect.outliers=TRUE,plot=TRUE,
                              from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5),sigdig=2,
                              alpha=0.05,exterr=TRUE,ranked=FALSE,
                              i2i=FALSE,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,rect.col=rect.col,
                        outlier.col=outlier.col,sigdig=sigdig,
                        alpha=alpha,exterr=exterr,i2i=i2i,
                        units='Ma',ranked=ranked,...)
}
#' @rdname weightedmean
#' @export
weightedmean.KCa <- function(x,random.effects=TRUE,
                             detect.outliers=TRUE,plot=TRUE,
                             from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                             outlier.col=rgb(0,1,1,0.5),sigdig=2,
                             alpha=0.05,exterr=TRUE,ranked=FALSE,
                             i2i=FALSE,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,rect.col=rect.col,
                        outlier.col=outlier.col,sigdig=sigdig,
                        alpha=alpha,exterr=exterr,i2i=i2i,
                        units='Ma',ranked=ranked,...)
}
#' @rdname weightedmean
#' @export
weightedmean.ReOs <- function(x,random.effects=TRUE,
                              detect.outliers=TRUE,plot=TRUE,
                              from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5),sigdig=2,
                              alpha=0.05,exterr=TRUE,ranked=FALSE,
                              i2i=TRUE,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,rect.col=rect.col,
                        outlier.col=outlier.col,sigdig=sigdig,
                        alpha=alpha,exterr=exterr,i2i=i2i,
                        units='Ma',ranked=ranked,...)
}
#' @rdname weightedmean
#' @export
weightedmean.SmNd <- function(x,random.effects=TRUE,
                              detect.outliers=TRUE,plot=TRUE,
                              from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5),sigdig=2,
                              alpha=0.05,exterr=TRUE,ranked=FALSE,
                              i2i=TRUE,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,rect.col=rect.col,
                        outlier.col=outlier.col,sigdig=sigdig,
                        alpha=alpha,exterr=exterr,i2i=i2i,
                        units='Ma',ranked=ranked,...)
}
#' @rdname weightedmean
#' @export
weightedmean.RbSr <- function(x,random.effects=TRUE,
                              detect.outliers=TRUE,plot=TRUE,
                              from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5),sigdig=2,
                              alpha=0.05,exterr=TRUE,
                              i2i=TRUE,ranked=FALSE,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,rect.col=rect.col,
                        outlier.col=outlier.col,sigdig=sigdig,
                        alpha=alpha,exterr=exterr,i2i=i2i,
                        units='Ma',ranked=ranked,...)
}
#' @rdname weightedmean
#' @export
weightedmean.LuHf <- function(x,random.effects=TRUE,
                              detect.outliers=TRUE,plot=TRUE,
                              from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5),sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=TRUE,
                              ranked=FALSE,...){
    weightedmean_helper(x,random.effects=random.effects,
                        detect.outliers=detect.outliers,plot=plot,
                        from=from,to=to,rect.col=rect.col,
                        outlier.col=outlier.col,sigdig=sigdig,
                        alpha=alpha,exterr=exterr,i2i=i2i,
                        units='Ma',ranked=ranked,...)
}
#' @rdname weightedmean
#' @export
weightedmean.UThHe <- function(x,random.effects=TRUE,
                               detect.outliers=TRUE,plot=TRUE,
                               from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                               outlier.col=rgb(0,1,1,0.5),sigdig=2,
                               alpha=0.05,ranked=FALSE,...){
    tt <- UThHe.age(x)
    fit <- weightedmean.default(tt,detect.outliers=detect.outliers,
                                alpha=alpha,plot=FALSE)
    if (plot){
        plot_weightedmean(tt[,1],tt[,2],fit=fit,from=from,to=to,
                          rect.col=rect.col,outlier.col=outlier.col,
                          sigdig=sigdig,alpha=alpha,units='Ma',
                          ranked=ranked,...)
    }
    invisible(fit)
}
#' @rdname weightedmean
#' @export
weightedmean.fissiontracks <- function(x,random.effects=TRUE,
                                       detect.outliers=TRUE,plot=TRUE,
                                       from=NA,to=NA,rect.col=rgb(0,1,0,0.5),
                                       outlier.col=rgb(0,1,1,0.5),
                                       sigdig=2,alpha=0.05,
                                       exterr=TRUE,ranked=FALSE,...){
    tt <- fissiontrack.age(x,exterr=FALSE)
    # calculated weighted mean age ignoring zeta and rhoD uncertainties
    fit <- weightedmean.default(tt,detect.outliers=detect.outliers,
                                alpha=alpha,plot=FALSE,...)
    out <- fit
    if (exterr){
        if (x$format==1) {
            rhoD <- x$rhoD
            zeta <- x$zeta
        } else if (x$format==2) {
            rhoD <- c(1,0)
            zeta <- x$zeta
        } else {
            rhoD <- c(1,0)
            zeta <- c(1,0)
        }
        out$mean['s[x]'] <- fit$mean['x']*
            sqrt( (fit$mean['s[x]']/fit$mean['x'])^2 +
                  (rhoD[2]/rhoD[1])^2 +
                  (zeta[2]/zeta[1])^2
                )
        out$mean['ci[x]'] <- nfact(alpha)*out$mean['s[x]']
    }
    if (plot){
        plot_weightedmean(tt[,1],tt[,2],fit=out,from=from,to=to,
                          rect.col=rect.col,outlier.col=outlier.col,
                          sigdig=sigdig,alpha=alpha,units='Ma',
                          ranked=ranked,...)
    }
    invisible(out)
}
weightedmean_helper <- function(x,random.effects=TRUE,detect.outliers=TRUE,
                                plot=TRUE,from=NA,to=NA,
                                rect.col=grDevices::rgb(0,1,0,0.5),type=4,
                                cutoff.76=1100,cutoff.disc=c(-15,5),
                                outlier.col=grDevices::rgb(0,1,1,0.5),
                                sigdig=2,alpha=0.05,exterr=TRUE,
                                ranked=FALSE,i2i=FALSE,common.Pb=1,
                                units='',detritus=0,Th02=c(0,0),
                                Th02U48=c(0,0,1e6,0,0,0,0,0,0),...){
    tt <- get.ages(x,type=type,cutoff.76=cutoff.76,
                   cutoff.disc=cutoff.disc,i2i=i2i,
                   detritus=detritus,Th02=Th02,Th02U48=Th02U48)
    fit <- weightedmean.default(tt,random.effects=random.effects,
                                detect.outliers=detect.outliers,
                                alpha=alpha,plot=FALSE)
    out <- fit
    if (exterr){
        out$mean[c('x','s[x]')] <-
            add.exterr(x,tt=fit$mean['x'],st=fit$mean['s[x]'],
                       cutoff.76=cutoff.76,type=type)
        out$mean['ci[x]'] <- nfact(alpha)*out$mean['s[x]']
    }
    if (plot){
        plot_weightedmean(tt[,1],tt[,2],from=from,to=to,fit=out,
                          rect.col=rect.col,outlier.col=outlier.col,
                          sigdig=sigdig,alpha=alpha,units=units,
                          ranked=ranked,...)
    }
    invisible(out)
}

get.weightedmean <- function(X,sX,random.effects=TRUE,
                             valid=TRUE,alpha=0.05){
    ns <- length(X)
    x <- X[valid]
    sx <- sX[valid]
    if (length(x)<=1){
        out$mean <- c(x,sx,nfact(alpha)*sx)
        out$mswd <- 0
        out$p.value <- 1
        return(out)
    }
    out <- list()
    out$random.effects <- random.effects
    out$df <- length(x)-1 # degrees of freedom for the homogeneity test
    if (random.effects){ # random effects model:
        out$mean <- rep(NA,3)
        out$disp <- rep(NA,3)
        names(out$mean) <- c('x','s[x]','ci[x]')
        names(out$disp) <- c('w','ll','ul')    
        fit <- continuous_mixture(x,sx)
        out$mean['x'] <- fit$mu[1]
        out$mean['s[x]'] <- fit$mu[2]
        out$mean['ci[x]'] <- nfact(alpha)*out$mean['s[x]']
        out$disp['w'] <- fit$sigma
        out$disp[c('ll','ul')] <-
            profile_LL_weightedmean_disp(fit,x,sx,alpha)
    } else { # Ludwig's Isoplot approach:
        out$mean <- rep(NA,4)
        names(out$mean) <- c('x','s[x]','ci[x]','disp[x]')
        w <- 1/sx^2
        out$mean['x'] <- sum(w*x)/sum(w)
        out$mean['s[x]'] <- 1/sqrt(sum(w))
        SS <- sum(((x-out$mean['x'])/sx)^2)
        out$mswd <- SS/out$df
        out$p.value <- 1-stats::pchisq(SS,out$df)
        out$mean['ci[x]'] <- tfact(alpha,out$df)*out$mean['s[x]']
        out$mean['disp[x]'] <- sqrt(out$mswd)*out$mean['ci[x]']
        out$valid <- valid
    }
    out$plotpar <-
        list(mean=list(x=c(0,ns+1),
                       y=rep(out$mean['x'],2)),
             rect=list(x=c(0,ns+1,ns+1,0),
                       y=c(rep(out$mean['x']+out$mean['ci[x]'],2),
                           rep(out$mean['x']-out$mean['ci[x]'],2)))
             )
    if (random.effects){
        out$plotpar$dash1=list(x=c(0,ns+1),y=rep(out$mean['x']+
                               nfact(alpha)*out$disp['w'],2))
        out$plotpar$dash2=list(x=c(0,ns+1),y=rep(out$mean['x']-
                               nfact(alpha)*out$disp['w'],2))
    }
    SS <- sum(((x-out$mean['x'])/sx)^2)
    out$mswd <- SS/out$df
    out$p.value <- 1-stats::pchisq(SS,out$df)
    out$valid <- valid
    out
}

wtdmean.title <- function(fit,sigdig=2,units='',...){
    if (fit$random.effects){
        rounded.mean <- roundit(fit$mean['x'],
                                fit$mean[c('s[x]','ci[x]')],
                                sigdig=sigdig)
        line1 <- substitute('mean ='~a%+-%b~'|'~c~u~'(n='~n/N~')',
                            list(a=rounded.mean['x'],
                                 b=rounded.mean['s[x]'],
                                 c=rounded.mean['ci[x]'],
                                 u=units,
                                 n=sum(fit$valid),
                                 N=length(fit$valid)))
        rounded.disp <- roundit(fit$disp['w'],
                                fit$disp[c('ll','ul')],
                                sigdig=sigdig)
        line3 <- substitute('dispersion ='~a+b/-c~u,
                            list(a=rounded.disp['w'],
                                 b=rounded.disp['ul'],
                                 c=rounded.disp['ll'],
                                 u=units))
        mymtext(line3,line=0,...)
        line1line <- 2
        line2line <- 1
    } else {
        rounded.mean <- roundit(fit$mean['x'],
                                fit$mean[c('s[x]','ci[x]','disp[x]')],
                                sigdig=sigdig)
        line1 <- substitute('mean ='~a%+-%b~'|'~c~'|'~d~u~'(n='~n/N~')',
                            list(a=rounded.mean['x'],
                                 b=rounded.mean['s[x]'],
                                 c=rounded.mean['ci[x]'],
                                 d=rounded.mean['disp[x]'],
                                 u=units,
                                 n=sum(fit$valid),
                                 N=length(fit$valid)))
        line1line <- 1
        line2line <- 0
    }
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=roundit(fit$mswd,fit$mswd,sigdig=sigdig),
                             b=roundit(fit$p.value,fit$p.value,
                                       sigdig=sigdig)))
    mymtext(line1,line=line1line,...)
    mymtext(line2,line=line2line,...)
}

plot_weightedmean <- function(X,sX,fit,from=NA,to=NA,
                              rect.col=grDevices::rgb(0,1,0,0.5),
                              outlier.col=grDevices::rgb(0,1,1,0.5),
                              sigdig=2,alpha=0.05,units='',
                              ranked=FALSE,...){
    if (ranked){
        i <- order(X)
        x <- X[i]
        sx <- sX[i]
        valid <- fit$valid[i]
    } else {
        x <- X
        sx <- sX
        valid <- fit$valid
    }
    ns <- length(x)
    fact <- nfact(alpha)
    if (is.na(from)) minx <- min(c(x-fact*sx,x-fact*fit$disp['w']),na.rm=TRUE)
    else minx <- from
    if (is.na(to)) maxx <- max(c(x+fact*sx,x+fact*fit$disp['w']),na.rm=TRUE)
    else maxx <- to
    graphics::plot(c(0,ns+1),c(minx,maxx),type='n',
                   axes=FALSE,xlab='N',ylab='',...)
    graphics::polygon(fit$plotpar$rect,col='gray80',border=NA)
    graphics::lines(fit$plotpar$mean)
    if (fit$random.effects){
        graphics::lines(fit$plotpar$dash1,lty=3)
        graphics::lines(fit$plotpar$dash2,lty=3)
    }
    graphics::axis(side=1,at=1:ns)
    graphics::axis(side=2)
    for (i in 1:ns){
        if (valid[i]) col <- rect.col
        else col <- outlier.col
        graphics::rect(xleft=i-0.4,ybottom=x[i]-fact*sx[i],
                       xright=i+0.4,ytop=x[i]+fact*sx[i],col=col)
    }
    graphics::title(wtdmean.title(fit,sigdig=sigdig,units=units))
}

# prune the data if necessary
# X and sX are some measurements and their standard errors
# valid is a vector of logical flags indicating whether the corresponding
# measurements have already been rejected or not
chauvenet <- function(X,sX,valid,random.effects=TRUE){
    fit <- get.weightedmean(X,sX,random.effects=random.effects,
                            valid=valid)
    mu <- fit$mean[1]
    if (random.effects) sigma <- sqrt(fit$disp[1]^2 + sX^2)
    # max(1,mswd) is an ad hoc solution to avoid dealing with underdispersed datasets
    else sigma <- sqrt(fit$mean[2]^2 + max(1,fit$mswd)*sX^2)
    prob <- 2*(1-stats::pnorm(abs(X-mu),sd=sigma))
    iprob <- order(prob)
    npruned <- length(which(!valid))
    imin <- iprob[npruned+1]
    minp <- prob[imin]
    ns <- length(which(valid))
    if (ns*minp < 0.5) {
        valid[imin] <- FALSE # remove outlier
    } 
    valid
}
