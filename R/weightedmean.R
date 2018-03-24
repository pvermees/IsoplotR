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
#'     \code{UPb}, \code{PbPb}, \code{ArAr}, \code{ReOs}, \code{SmNd},
#'     \code{RbSr}, \code{LuHf}, \code{ThU}, \code{fissiontracks} or
#'     \code{UThHe}
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
#' @param rect.col the fill colour of the rectangles used to show the
#'     measurements or age estimates.
#' @param outlier.col if \code{detect.outliers=TRUE}, the outliers are
#'     given a different colour.
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output.
#' @param alpha the confidence limits of the error bars/rectangles.
#' @importFrom grDevices rgb
#' @rdname weightedmean
#' @export
weightedmean.default <- function(x,detect.outliers=TRUE,plot=TRUE,
                                 rect.col=rgb(0,1,0,0.5),
                                 outlier.col=rgb(0,1,1,0.5),
                                 sigdig=2,alpha=0.05,...){
    X <- x[,1]
    sX <- x[,2]
    valid <- !is.na(X) & !is.na(sX)
    nvalid <- count(valid)
    if (detect.outliers){
        while (TRUE){
            valid <- chauvenet(X,sX,valid)
            if (count(valid) < nvalid)
                nvalid <- count(valid)
            else
                break
        }
    }
    out <- get.weightedmean(X,sX,valid=valid,alpha=alpha)
    ns <- length(X)
    out$plotpar <-
        list(mean=list(x=c(0,ns+1),
                       y=rep(out$mean['x'],2)),
             rect=list(x=c(0,ns+1,ns+1,0),
                       y=c(rep(out$mean['x']+out$mean['ci[x]'],2),
                           rep(out$mean['x']-out$mean['ci[x]'],2))),
             dash1=list(x=c(0,ns+1),
                        y=rep(out$mean['x']+nfact(alpha)*out$disp['w'],2)),
             dash2=list(x=c(0,ns+1),
                        y=rep(out$mean['x']-nfact(alpha)*out$disp['w'],2))
             )
    if (plot){
        plot_weightedmean(X,sX,out,rect.col=rect.col,
                          outlier.col=outlier.col,sigdig=sigdig,
                          alpha=alpha,...)
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
weightedmean.UPb <- function(x,detect.outliers=TRUE,plot=TRUE,
                             rect.col=rgb(0,1,0,0.5),
                             outlier.col=rgb(0,1,1,0.5),
                             sigdig=2,type=4,cutoff.76=1100,
                             cutoff.disc=c(-15,5),alpha=0.05,
                             exterr=TRUE,common.Pb=0,...){
    if (common.Pb %in% c(1,2,3))
        X <- common.Pb.correction(x,option=common.Pb)
    else
        X <- x
    weightedmean_helper(X,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        type=type,cutoff.76=cutoff.76,
                        cutoff.disc=cutoff.disc,sigdig=sigdig,
                        alpha=alpha,exterr=exterr,units='Ma',...)
}
#' @rdname weightedmean
#' @export
weightedmean.PbPb <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,common.Pb=1,...){
    if (common.Pb %in% c(1,2,3))
        X <- common.Pb.correction(x,option=common.Pb)
    else
        X <- x
    weightedmean_helper(X,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        units='Ma',...)
}
#' @param i2i `isochron to intercept': calculates the initial (aka
#'     `inherited', `excess', or `common')
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar, \eqn{^{207}}Pb/\eqn{^{204}}Pb,
#'     \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd,
#'     \eqn{^{187}}Os/\eqn{^{188}}Os or \eqn{^{176}}Hf/\eqn{^{177}}Hf
#'     ratio from an isochron fit. Setting \code{i2i} to \code{FALSE}
#'     uses the default values stored in
#'     \code{settings('iratio',...)}. When applied to data of class
#'     \code{ThU}, setting \code{i2i} to \code{TRUE} applies a
#'     detrital Th-correction.
#' @rdname weightedmean
#' @export
weightedmean.ThU <- function(x,detect.outliers=TRUE,plot=TRUE,
                             rect.col=rgb(0,1,0,0.5),
                             outlier.col=rgb(0,1,1,0.5), sigdig=2,
                             alpha=0.05,i2i=TRUE,...){
    weightedmean_helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,i2i=i2i,units='ka',...)
}
#' @rdname weightedmean
#' @export
weightedmean.ArAr <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=FALSE,...){
    weightedmean_helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,units='Ma',...)
}
#' @rdname weightedmean
#' @export
weightedmean.ReOs <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=TRUE,...){
    weightedmean_helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,units='Ma',...)
}
#' @rdname weightedmean
#' @export
weightedmean.SmNd <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=TRUE,...){
    weightedmean_helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,units='Ma',...)
}
#' @rdname weightedmean
#' @export
weightedmean.RbSr <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=TRUE,...){
    weightedmean_helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,units='Ma',...)
}
#' @rdname weightedmean
#' @export
weightedmean.LuHf <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=TRUE,...){
    weightedmean_helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,units='Ma',...)
}
#' @rdname weightedmean
#' @export
weightedmean.UThHe <- function(x,detect.outliers=TRUE,plot=TRUE,
                               rect.col=rgb(0,1,0,0.5),
                               outlier.col=rgb(0,1,1,0.5),sigdig=2,
                               alpha=0.05,...){
    tt <- UThHe.age(x)
    fit <- weightedmean.default(tt,detect.outliers=detect.outliers,
                                alpha=alpha,plot=FALSE,...)
    if (plot){
        plot_weightedmean(tt[,1],tt[,2],fit,rect.col=rect.col,
                          outlier.col=outlier.col,sigdig=sigdig,
                          alpha=alpha,units='Ma')
    }
    invisible(fit)
}
#' @rdname weightedmean
#' @export
weightedmean.fissiontracks <- function(x,detect.outliers=TRUE,plot=TRUE,
                                       rect.col=rgb(0,1,0,0.5),
                                       outlier.col=rgb(0,1,1,0.5),
                                       sigdig=2,alpha=0.05,exterr=TRUE,...){
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
        plot_weightedmean(tt[,1],tt[,2],out,rect.col=rect.col,
                          outlier.col=outlier.col,sigdig=sigdig,
                          alpha=alpha,units='Ma')
    }
    invisible(out)
}
weightedmean_helper <- function(x,detect.outliers=TRUE,plot=TRUE,
                                rect.col=grDevices::rgb(0,1,0,0.5),
                                type=4,cutoff.76=1100,
                                cutoff.disc=c(-15,5),
                                outlier.col=grDevices::rgb(0,1,1,0.5),
                                sigdig=2,alpha=0.05,exterr=TRUE,
                                i2i=FALSE,common.Pb=1,units='',...){
    if (hasClass(x,'UPb')){
        tt <- filter.UPb.ages(x,type=type,cutoff.76=cutoff.76,
                              cutoff.disc=cutoff.disc,exterr=FALSE)
    } else if (hasClass(x,'PbPb')){
        tt <- PbPb.age(x,exterr=FALSE)
    } else if (hasClass(x,'ThU')){
        tt <- ThU.age(x,exterr=FALSE,i2i=i2i)
        exterr <- FALSE
    } else if (hasClass(x,'ArAr')){
        tt <- ArAr.age(x,jcu=FALSE,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'ReOs')){
        tt <- ReOs.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'SmNd')){
        tt <- SmNd.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'RbSr')){
        tt <- RbSr.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'LuHf')){
        tt <- LuHf.age(x,exterr=FALSE,i2i=i2i)
    }
    fit <- weightedmean.default(tt,detect.outliers=detect.outliers,
                                alpha=alpha,plot=FALSE,...)
    out <- fit
    if (exterr){
        out$mean[c('x','s[x]')] <-
            add.exterr(x,tt=fit$mean['x'],st=fit$mean['s[x]'],
                       cutoff.76=cutoff.76,type=type)
        out$mean['ci[x]'] <- nfact(alpha)*out$mean['s[x]']
    }
    if (plot){
        plot_weightedmean(tt[,1],tt[,2],out,rect.col=rect.col,
                          outlier.col=outlier.col,sigdig=sigdig,
                          alpha=alpha,units=units)
    }
    invisible(out)
}

get.weightedmean <- function(X,sX,valid=TRUE,alpha=0.05){
    x <- X[valid]
    sx <- sX[valid]
    out <- list()
    out$mean <- rep(NA,3)
    out$disp <- rep(NA,3)
    out$df <- length(x)-1 # degrees of freedom for the homogeneity test
    names(out$mean) <- c('x','s[x]','ci[x]')
    names(out$disp) <- c('w','ll','ul')    
    if (length(x)>1){
        fit <- continuous_mixture(x,sx)
        out$mean['x'] <- fit$mu[1]
        out$mean['s[x]'] <- fit$mu[2]
        out$mean['ci[x]'] <- nfact(alpha)*out$mean['s[x]']
        out$disp['w'] <- fit$sigma
        out$disp[c('ll','ul')] <-
            profile_LL_weightedmean_disp(fit,x,sx,alpha)
        SS <- sum(((x-out$mean['x'])/sx)^2)
        out$mswd <- SS/out$df
        out$p.value <- 1-stats::pchisq(SS,out$df)
        out$valid <- valid
    } else {
        out$mean <- x
        out$p.value <- 0
    }
    out
}

wtdmean.title <- function(fit,sigdig=2,units=''){
    rounded.mean <- roundit(fit$mean['x'],
                            fit$mean[c('s[x]','ci[x]')],
                            sigdig=sigdig)
    rounded.disp <- roundit(fit$disp['w'],
                            fit$disp[c('ll','ul')],
                            sigdig=sigdig)
    line1 <- substitute('mean ='~a%+-%b~'|'~c~u,
                        list(a=rounded.mean['x'],
                             b=rounded.mean['s[x]'],
                             c=rounded.mean['ci[x]'],
                             u=units))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,sigdig),
                             b=signif(fit$p.value,sigdig)))
    line3 <- substitute('dispersion ='~a+b/-c~u,
                        list(a=rounded.disp['w'],
                             b=rounded.disp['ul'],
                             c=rounded.disp['ll'],
                             u=units))
    graphics::mtext(line1,line=2)
    graphics::mtext(line2,line=1)
    graphics::mtext(line3,line=0)
}

plot_weightedmean <- function(X,sX,fit,
                              rect.col=grDevices::rgb(0,1,0,0.5),
                              outlier.col=grDevices::rgb(0,1,1,0.5),
                              sigdig=2,alpha=0.05,units=''){
    ns <- length(X)
    fact <- nfact(alpha)
    minX <- min(c(X-fact*sX,X-fact*fit$disp['w']),na.rm=TRUE)
    maxX <- max(c(X+fact*sX,X+fact*fit$disp['w']),na.rm=TRUE)
    graphics::plot(c(0,ns+1),c(minX,maxX),type='n',
                   axes=FALSE,xlab='N',ylab='')
    graphics::polygon(fit$plotpar$rect,col='gray80',border=NA)
    graphics::lines(fit$plotpar$mean)
    graphics::lines(fit$plotpar$dash1,lty=3)
    graphics::lines(fit$plotpar$dash2,lty=3)
    graphics::axis(side=1,at=1:ns)
    graphics::axis(side=2)
    for (i in 1:ns){
        if (fit$valid[i]) col <- rect.col
        else col <- outlier.col
        graphics::rect(xleft=i-0.4,ybottom=X[i]-fact*sX[i],
                       xright=i+0.4,ytop=X[i]+fact*sX[i],col=col)
    }
    graphics::title(wtdmean.title(fit,sigdig=sigdig,units=units))
}

# prune the data if necessary
# X and sX are some measurements and their standard errors
# valid is a vector of logical flags indicating whether the corresponding
# measurements have already been rejected or not
chauvenet <- function(X,sX,valid){
    fit <- get.weightedmean(X,sX,valid=valid)
    mu <- fit$mean[1]
    sigma <- sqrt(fit$disp[1]^2+sX^2)
    prob <- 2*(1-stats::pnorm(abs(X-mu),sd=sigma))
    minp <- min(prob[valid])
    imin <- which.min(sX/abs(X-mu))
    ns <- length(valid[valid])
    if (ns*minp < 0.5) valid[imin] <- FALSE # remove outlier
    valid
}
