#' Calculate the weighted mean age
#'
#' Models the data as a Normal distribution with two sources of
#' variance.  Estimates the mean and 'overdispersion' using the method
#' of Maximum Likelihood. Computes the MSWD of a Normal fit without
#' overdispersion. Implements Chauvenet's Criterion to detect and
#' reject outliers. Only propagates the analytical uncertainty
#' associated with decay constants and J-factors after computing the
#' weighted mean isotopic composition.
#' 
#' @param x a two column matrix of values (first column) and their
#'     standard errors (second column) OR an object of class
#'     \code{UPb}, \code{PbPb}, \code{ArAr}, \code{ReOs}, \code{SmNd},
#'     \code{RbSr}, \code{LuHf}, \code{ThU}, \code{fissiontracks} or
#'     \code{UThHe}
#' @param ... optional arguments
#' @return if \code{PLOT=FALSE}, returns a list with the following
#'     items:
#'
#' \describe{
#' \item{mean}{a two element vector with the weighted mean and its
#' standard error.}
#'
#' \item{disp}{a two element vector with the (over)dispersion and its
#' standard error.}
#'
#' \item{mswd}{the Mean Square of the Weighted Deviates
#' (a.k.a. `reduced Chi-square' statistic)}
#'
#' \item{p.value}{the p-value of a Chi-square test with n-1 degrees
#' of freedom, testing the null hypothesis that the underlying
#' population is not overdispersed.}
#'
#' \item{valid}{vector of logical flags indicating which steps are
#' included into the weighted mean calculation}
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
#' @rdname weightedmean
#' @export
weightedmean.default <- function(x,detect.outliers=TRUE,plot=TRUE,
                                 rect.col=rgb(0,1,0,0.5),
                                 outlier.col=rgb(0,1,1,0.5),
                                 sigdig=2,alpha=0.05,...){
    X <- x[,1]
    sX <- x[,2]
    ns <- length(X)
    valid <- !is.na(X) & !is.na(sX)
    if (detect.outliers){
        prob <- stats::pnorm(X,mean=mean(X,na.rm=TRUE),
                             sd=stats::sd(X,na.rm=TRUE))
        cutoff <- 0.5/ns
        valid <- valid & (prob > cutoff) & ((1-prob) > cutoff)
    }
    fit <- get.weightedmean(X,sX,valid)
    if (plot){
        plot.weightedmean(X,sX,fit,rect.col=rect.col,
                          outlier.col=outlier.col,sigdig=sigdig,
                          alpha=alpha,...)
    } else {
        return(fit)
    }
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
#' @param exterr propagate decay constant uncertainty?
#' @examples
#' ages <- c(251.9,251.59,251.47,251.35,251.1,251.04,250.79,250.73,251.22,228.43)
#' errs <- c(0.28,0.28,0.63,0.34,0.28,0.63,0.28,0.4,0.28,0.33)
#' weightedmean(cbind(ages,errs))
#' data(examples)
#' weightedmean(examples$ArAr)
#' @rdname weightedmean
#' @export
weightedmean.UPb <- function(x,detect.outliers=TRUE,plot=TRUE,
                             rect.col=rgb(0,1,0,0.5),
                             outlier.col=rgb(0,1,1,0.5),
                             sigdig=2,type=4,cutoff.76=1100,
                             cutoff.disc=c(-15,5),alpha=0.05,
                             exterr=TRUE,...){
    weightedmean.helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        type=type,cutoff.76=cutoff.76,
                        cutoff.disc=cutoff.disc,sigdig=sigdig,
                        alpha=alpha,exterr=exterr,...)
}
#' @param i2i `isochron to intercept': calculates the initial
#'     (aka `inherited', `excess', or `common') \eqn{^{40}}Ar/\eqn{^{36}}Ar,
#'     \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#'     \eqn{^{143}}Nd/\eqn{^{144}}Nd, \eqn{^{187}}Os/\eqn{^{188}}Os or
#'     \eqn{^{176}}Hf/\eqn{^{177}}Hf ratio from an isochron
#'     fit. Setting \code{i2i} to \code{FALSE} uses the default values
#'     stored in \code{settings('iratio',...)}  or zero (for the Pb-Pb
#'     method). When applied to data of class \code{ThU}, setting
#'     \code{i2i} to \code{TRUE} applies a detrital Th-correction.
#' @rdname weightedmean
#' @export
weightedmean.PbPb <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=FALSE,...){
    weightedmean.helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,...)
}
#' @rdname weightedmean
#' @export
weightedmean.ThU <- function(x,detect.outliers=TRUE,plot=TRUE,
                             rect.col=rgb(0,1,0,0.5),
                             outlier.col=rgb(0,1,1,0.5), sigdig=2,
                             alpha=0.05,i2i=TRUE,...){
    weightedmean.helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=FALSE,
                        i2i=i2i,...)
}
#' @rdname weightedmean
#' @export
weightedmean.ArAr <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=FALSE,...){
    weightedmean.helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,...)
}
#' @rdname weightedmean
#' @export
weightedmean.ReOs <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=TRUE,...){
    weightedmean.helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,...)
}
#' @rdname weightedmean
#' @export
weightedmean.SmNd <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=TRUE,...){
    weightedmean.helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,...)
}
#' @rdname weightedmean
#' @export
weightedmean.RbSr <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=TRUE,...){
    weightedmean.helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,...)
}
#' @rdname weightedmean
#' @export
weightedmean.LuHf <- function(x,detect.outliers=TRUE,plot=TRUE,
                              rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5), sigdig=2,
                              alpha=0.05,exterr=TRUE,i2i=TRUE,...){
    weightedmean.helper(x,detect.outliers=detect.outliers,plot=plot,
                        rect.col=rect.col,outlier.col=outlier.col,
                        sigdig=sigdig,alpha=alpha,exterr=exterr,
                        i2i=i2i,...)
}
#' @rdname weightedmean
#' @export
weightedmean.UThHe <- function(x,detect.outliers=TRUE,plot=TRUE,
                               rect.col=rgb(0,1,0,0.5),
                               outlier.col=rgb(0,1,1,0.5),
                               sigdig=2,alpha=0.05,exterr=TRUE,...){
    tt <- UThHe.age(x)
    fit <- weightedmean.default(tt,detect.outliers=detect.outliers,plot=FALSE,...)
    if (plot){
        plot.weightedmean(tt[,1],tt[,2],fit,rect.col=rect.col,
                          outlier.col=outlier.col,sigdig=sigdig,
                          alpha=alpha)
    } else {
        return(fit)
    }
}
#' @rdname weightedmean
#' @export
weightedmean.fissiontracks <- function(x,detect.outliers=TRUE,plot=TRUE,
                                       rect.col=rgb(0,1,0,0.5),
                                       outlier.col=rgb(0,1,1,0.5),
                                       sigdig=2,alpha=0.05,exterr=TRUE,...){
    tt <- fissiontrack.age(x,exterr=FALSE)
    # calculated weighted mean age ignoring zeta and rhoD uncertainties
    fit <- weightedmean.default(tt,detect.outliers=detect.outliers,plot=FALSE,...)
    if (exterr){
        stt2 <- (fit$mean[2]/fit$mean[1])^2
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
        fit$mean[2] <- fit$mean[1] *
            sqrt( stt2 + (rhoD[2]/rhoD[1])^2 + (zeta[2]/zeta[1])^2 )
    }
    if (plot){
        plot.weightedmean(tt[,1],tt[,2],fit,rect.col=rect.col,
                          outlier.col=outlier.col,sigdig=sigdig,
                          alpha=alpha)
    } else {
        return(fit)
    }
}
weightedmean.helper <- function(x,detect.outliers=TRUE,plot=TRUE,
                                rect.col=rgb(0,1,0,0.5),type=4,
                                cutoff.76=1100,cutoff.disc=c(-15,5),
                                outlier.col=rgb(0,1,1,0.5), sigdig=2,
                                alpha=0.05,exterr=TRUE,i2i=FALSE,...){
    if (hasClass(x,'UPb')){
        tt <- filter.UPb.ages(x,type=type,cutoff.76=cutoff.76,
                              cutoff.disc=cutoff.disc,exterr=FALSE)
    } else if (hasClass(x,'PbPb')){
        tt <- PbPb.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'ThU')){
        tt <- ThU.age(x,exterr=FALSE,i2i=i2i)
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
    fit <- weightedmean.default(tt,detect.outliers=detect.outliers,plot=FALSE,...)
    if (exterr){
        fit$mean <- add.exterr(x,tt=fit$mean[1],st=fit$mean[2],
                               cutoff.76=cutoff.76,type=type)
    }
    if (plot){
        plot.weightedmean(tt[,1],tt[,2],fit,rect.col=rect.col,
                          outlier.col=outlier.col,sigdig=sigdig,alpha=alpha)
    } else {
        return(fit)
    }
}

get.weightedmean <- function(X,sX,valid=TRUE){
    out <- list()
    X <- X[valid]
    sX <- sX[valid]
    if (length(X)>1){
        MZ <- c(mean(X,na.rm=TRUE),stats::sd(X,na.rm=TRUE))
        fit <- stats::optim(MZ,LL.weightedmean,X=X,sX=sX,
                            method='BFGS',hessian=TRUE)
        covmat <- solve(fit$hessian)
        out$mean <- c(fit$par[1],sqrt(covmat[1,1]))
        out$disp <- c(fit$par[2],sqrt(covmat[2,2]))
        df <- length(X)-1
        SS <- sum(((X-out$mean[1])/sX)^2)
        out$mswd <- SS/df
        out$p.value <- 1-pchisq(SS,df)
        out$valid <- valid
    } else {
        out$mean <- X
        out$p.value <- 0
    }
    out
}

LL.weightedmean <- function(MZ,X,sX){
    M <- MZ[1] # Mu (mean)
    Z <- MZ[2] # Zeta (overdispersion)
    out <- 0
    ns <- length(X)
    -sum(stats::dnorm(X,mean=M,sd=sqrt(sX^2+Z^2),log=TRUE))
}

wtdmean.title <- function(fit,sigdig=2){
    rounded.mean <- roundit(fit$mean[1],fit$mean[2],sigdig=sigdig)
    line1 <- substitute('mean ='~a%+-%b~' (1'~sigma~')',
                        list(a=rounded.mean[1], b=rounded.mean[2]))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,sigdig),
                             b=signif(fit$p.value,sigdig)))
    graphics::mtext(line1,line=2)
    graphics::mtext(line2,line=1)
    if (fit$p.value < 0.05){ # only show when the data are overdispersed
        rounded.disp <- roundit(fit$disp[1],fit$disp[2],sigdig=sigdig)
        line3 <- substitute('overdispersion ='~a%+-%b~' (1'~sigma~')',
                            list(a=rounded.disp[1], b=rounded.disp[2]))
        graphics::mtext(line3,line=0)
    }
}

plot.weightedmean <- function(X,sX,fit,rect.col=rgb(0,1,0,0.5),
                              outlier.col=rgb(0,1,1,0.5),sigdig=2,
                              alpha=0.05){
    ns <- length(X)
    fact <- stats::qnorm(1-alpha/2)
    minX <- min(X-fact*sX,na.rm=TRUE)
    maxX <- max(X+fact*sX,na.rm=TRUE)
    graphics::plot(c(0,ns+1),c(minX,maxX),type='n',axes=FALSE,xlab='N',ylab='')
    graphics::lines(c(0,ns+1),c(fit$mean[1],fit$mean[1]))
    graphics::axis(side=1,at=1:ns)
    graphics::axis(side=2)
    for (i in 1:ns){
        if (fit$valid[i]) col <- rect.col
        else col <- outlier.col
        graphics::rect(xleft=i-0.4,ybottom=X[i]-fact*sX[i],
                       xright=i+0.4,ytop=X[i]+fact*sX[i],col=col)
    }
    title(wtdmean.title(fit,sigdig=sigdig))
}
