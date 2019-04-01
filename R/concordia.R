#' Concordia diagram
#'
#' Plots U-Pb data on Wetherill and Tera-Wasserburg concordia
#' diagrams, calculate concordia ages and compositions, evaluates the
#' equivalence of multiple
#' (\eqn{^{206}}Pb/\eqn{^{238}}U-\eqn{^{207}}Pb/\eqn{^{235}}U or
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U)
#' compositions, computes the weighted mean isotopic composition and
#' the corresponding concordia age using the method of maximum
#' likelihood, computes the MSWD of equivalence and concordance and
#' their respective Chi-squared p-values. Performs linear regression
#' and computes the upper and lower intercept ages (for Wetherill) or
#' the lower intercept age and the \eqn{^{207}}Pb/\eqn{^{206}}Pb
#' intercept (for Tera-Wasserburg), taking into account error
#' correlations and decay constant uncertainties.
#'
#' @details
#' The concordia diagram is a graphical means of assessing the
#' internal consistency of U-Pb data. It sets out the measured
#' \eqn{^{206}}Pb/\eqn{^{238}}U- and
#' \eqn{^{207}}Pb/\eqn{^{235}}U-ratios against each other (`Wetherill'
#' diagram) or, equivalently, the \eqn{^{207}}Pb/\eqn{^{206}}Pb- and
#' \eqn{^{206}}Pb/\eqn{^{238}}U-ratios (`Tera-Wasserburg'
#' diagram). The space of concordant isotopic compositions is marked
#' by a curve, the `concordia line'. Isotopic ratio measurements are
#' shown as 100(1-\code{alpha})\% confidence ellipses. Concordant
#' samples plot near to, or overlap with, the concordia line. They
#' represent the pinnacle of geochronological robustness. Samples that
#' plot away from the concordia line but are aligned along a linear
#' trend form an isochron (or `discordia' line) that can be used
#' to infer the composition of the non-radiogenic
#' (`common') lead or to constrain the timing of prior lead loss.
#'
#' @param x an object of class \code{UPb}
#' 
#' @param tlim age limits of the concordia line
#' 
#' @param alpha probability cutoff for the error ellipses and
#'     confidence intervals
#' 
#' @param wetherill logical flag (\code{FALSE} for Tera-Wasserburg)
#' 
#' @param show.numbers logical flag (\code{TRUE} to show grain
#'     numbers)
#' 
#' @param levels a vector with \code{length(x)} values to be displayed
#'     as different background colours within the error ellipses.
#' 
#' @param clabel label for the colour legend (only used if
#'     \code{levels} is not \code{NA}.
#' 
#' @param ellipse.col a vector of two background colours for the error
#'     ellipses. If \code{levels=NA}, then only the first colour is
#'     used. If \code{levels} is a vector of numbers, then
#'     \code{ellipse.col} is used to construct a colour ramp.
#' 
#' @param concordia.col colour of the concordia line
#' 
#' @param exterr show decay constant uncertainty?
#' 
#' @param show.age one of either:
#'
#' \code{0}: plot the data without calculating an age
#'
#' \code{1}: fit a concordia composition and age
#'
#' \code{2}: fit a discordia line through the data using the maximum
#' likelihood algorithm of Ludwig (1998), which assumes that the
#' scatter of the data is solely due to the analytical
#' uncertainties. In this case, \code{IsoplotR} will either calculate
#' an upper and lower intercept age (for Wetherill concordia), or a
#' lower intercept age and common
#' (\eqn{^{207}}Pb/\eqn{^{206}}Pb)-ratio intercept (for
#' Tera-Wasserburg). If \code{mswd}>0, then the analytical
#' uncertainties are augmented by a factor \eqn{\sqrt{mswd}}.
#'
#' \code{3}: fit a discordia line ignoring the analytical uncertainties
#'
#' \code{4}: fit a discordia line using a modified maximum likelihood
#' algorithm that includes accounts for any overdispersion by adding a
#' geological (co)variance term.
#'
#' @param sigdig number of significant digits for the
#'     concordia/discordia age
#' 
#' @param common.Pb apply a common lead correction using one of three
#'     methods:
#'
#' \code{1}: use the Stacey-Kramers two-stage model to infer the initial
#' Pb-composition
#'
#' \code{2}: use the isochron intercept as the initial Pb-composition
#'
#' \code{3}: use the Pb-composition stored in
#' \code{settings('iratio','Pb207Pb206')} (if \code{x$format}<4) or
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')} (if \code{x$format}>3)
#' @param anchor control parameters to fix the intercept age or common
#'     Pb composition of the discordia fit. This is a two-element
#'     list.
#'

#' \itemize{
#'
#' \item The first element is a boolean flag indicating whether the
#' discordia line should be anchored. If this is \code{FALSE}, then
#' the second item is ignored and both the common Pb composition and
#' age are estimated.
#'
#' \item If the first element is \code{TRUE} and the second element is
#' \code{NA}, then the common Pb composition is fixed at the values
#' stored in \code{settings('iratio',...)}.
#'
#' item If the first element is \code{TRUE} and the second element is a
#' number, then the discordia line is forced to intersect the
#' concordia line at an age equal to that number.
#' }
#'
#' @param ticks either a scalar indicating the desired number of age
#'     ticks to be placed along the concordia line, OR a vector of
#'     tick ages.
#' @param hide vector with indices of aliquots that should be removed
#'     from the concordia diagram
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from concordia or discordia age calculation
#' @param omit.col colour that should be used for the omitted
#'     aliquots.
#' @param ... optional arguments to the generic \code{plot} function
#'
#' @return
#'
#' if \code{show.age=1}, returns a list with the following items:
#'
#' \describe{
#'
#' \item{x}{ a named vector with the (weighted mean) U-Pb composition }
#'
#' \item{cov}{ the covariance matrix of the (weighted mean) U-Pb composition }
#'
#' \item{mswd}{ a vector with three items (\code{equivalence},
#' \code{concordance} and \code{combined}) containing the MSWD (Mean
#' of the Squared Weighted Deviates, a.k.a the reduced Chi-squared
#' statistic) of isotopic equivalence, age concordance and combined
#' goodness of fit, respectively. }
#'
#' \item{p.value}{ a vector with three items (\code{equivalence},
#' \code{concordance} and \code{combined}) containing the p-value of
#' the Chi-square test for isotopic equivalence, age concordance and
#' combined goodness of fit, respectively. }
#'
#' \item{df}{ a three-element vector with the number of degrees of
#' freedom used for the \code{mswd} calculation.  These values are
#' useful when expanding the analytical uncertainties if
#' \code{mswd>1}.
#'
#' }
#'
#' \item{age}{a 4-element vector with:
#'
#' \code{t}: the concordia age (in Ma)
#'
#' \code{s[t]}: the estimated uncertainty of \code{t}
#'
#' \code{ci[t]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval of \code{t} for the appropriate degrees of freedom
#'
#' \code{disp[t]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{t} augmented by \eqn{\sqrt{mswd}} to account for
#' overdispersed datasets.}
#'
#' }
#'
#' if \code{show.age=2}, \code{3} or \code{4}, returns a list with the
#' following items:
#'
#' \describe{
#'
#' \item{model}{ the fitting model (\code{=show.age-1}).}
#'
#' \item{x}{ a two-element vector with the upper and lower intercept
#' ages (if \code{wetherill=TRUE}) or the lower intercept age and
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb intercept (if
#' \code{wetherill=FALSE}).}
#'
#' \item{cov}{ the covariance matrix of the elements in \code{x}.}
#'
#' \item{err}{ a \code{[2 x 2]} or \code{[3 x 2]} matrix with the
#' following rows:
#'
#' \code{s}: the estimated standard deviation for \code{x}
#'
#' \code{ci}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval of \code{x} for the appropriate degrees of freedom
#'
#' \code{disp[t]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{x} augmented by \eqn{\sqrt{mswd}} to account for
#' overdispersed datasets (only reported if \code{show.age=2}).  }
#'
#' \item{df}{ the degrees of freedom of the concordia fit (concordance
#' + equivalence)}
#'
#' \item{p.value}{ p-value of a Chi-square test for age homogeneity
#' (only reported if \code{ type=3}).}
#'
#' \item{mswd}{ mean square of the weighted deviates -- a
#' goodness-of-fit measure. \code{mswd > 1} indicates overdispersion
#' w.r.t the analytical uncertainties (not reported if
#' \code{show.age=3}).}
#'
#' \item{w}{ three-element vector with the standard deviation of the
#' (assumedly) Normal overdispersion term and the lower and upper
#' half-widths of its \eqn{100(1-\alpha)\%} confidence interval (only
#' important if \code{show.age=4}).}
#'
#' \item{n}{ the number of aliquots in the dataset }
#'
#' }
#'
#' @examples
#' data(examples)
#' concordia(examples$UPb,show.age=2)
#'
#' dev.new()
#' concordia(examples$UPb,wetherill=FALSE,
#'           xlim=c(24.9,25.4),ylim=c(0.0508,0.0518),
#'           ticks=249:254,exterr=TRUE)
#'
#' dev.new()
#' concordia(examples$UPb,wetherill=FALSE,show.age=2,anchor=list(TRUE,0))
#'
#' @references Ludwig, K.R., 1998. On the treatment of concordant
#'     uranium-lead ages. Geochimica et Cosmochimica Acta, 62(4),
#'     pp.665-676.
#' @export
concordia <- function(x=NULL,tlim=NULL,alpha=0.05,wetherill=TRUE,
                      show.numbers=FALSE,levels=NA,clabel="",
                      ellipse.col=c("#00FF0080","#FF000080"),
                      concordia.col='darksalmon',exterr=FALSE,
                      show.age=0,sigdig=2,common.Pb=0,ticks=5,
                      anchor=list(FALSE,NA),hide=NULL,omit=NULL,
                      omit.col=NA,...){
    if (is.null(x)){
        emptyconcordia(tlim=tlim,alpha=alpha,wetherill=wetherill,exterr=exterr,
                       concordia.col=concordia.col,ticks=ticks,...)
        return(invisible(NULL))
    }
    ns <- length(x)
    plotit <- (1:ns)%ni%hide
    calcit <- (1:ns)%ni%c(hide,omit)
    if (common.Pb<1) X <- x
    else X <- Pb0corr(x,option=common.Pb,omit=unique(c(hide,omit)))
    X2plot <- subset(X,subset=plotit)
    lims <- prepare.concordia.line(X2plot,tlim=tlim,wetherill=wetherill,...)
    fit <- NULL
    if (show.age>1){
        x2calc <- subset(x,subset=calcit)
        fit <- concordia.intersection.ludwig(x2calc,wetherill=wetherill,exterr=exterr,
                                             alpha=alpha,model=(show.age-1),anchor=anchor)
        discordia.line(fit,wetherill=wetherill,d=x$d)
        fit$n <- length(x2calc)
        graphics::title(discordia.title(fit,wetherill=wetherill,sigdig=sigdig))
    }
    plot.concordia.line(X2plot,lims=lims,wetherill=wetherill,col=concordia.col,
                        alpha=alpha,exterr=exterr,ticks=ticks)
    y <- data2york(X,option=(2-1*wetherill))
    scatterplot(y,alpha=alpha,show.numbers=show.numbers,
                show.ellipses=1*(show.age!=3),levels=levels,
                clabel=clabel,ellipse.col=ellipse.col,add=TRUE,
                hide=hide,omit=omit,omit.col=omit.col,addcolourbar=FALSE,...)
    if (show.age==1){
        X2calc <- subset(X,subset=calcit)
        fit <- concordia.age(X2calc,wetherill=wetherill,exterr=exterr,alpha=alpha)
        ell <- ellipse(fit$x[1],fit$x[2],fit$cov)
        graphics::polygon(ell,col='white')
        fit$n <- length(X2calc)
        graphics::title(concordia.title(fit,sigdig=sigdig))
    }
    colourbar(z=levels[calcit],col=ellipse.col,clabel=clabel)
    invisible(fit)
}

# helper function for plot.concordia
plot.concordia.line <- function(x,lims,wetherill=TRUE,col='darksalmon',
                                alpha=0.05,exterr=TRUE,ticks=5){
    if (all(is.null(x$x))){
        m <- 0
        M <- 4500
    } else {
        range.t <- range(lims$t)
        m <- max(0.8*lims$t[1],lims$t[1]-range.t/20)
        M <- min(1.2*lims$t[2],lims$t[2]+range.t/20)
    }
    nn <- 30 # number of segments into which the concordia line is divided
    tt <- cseq(m,M,wetherill=wetherill,n=nn)
    conc <- matrix(0,nn,2)
    colnames(conc) <- c('x','y')
    for (i in 1:nn){ # build the concordia line
        xy <- age_to_concordia_ratios(tt[i],wetherill=wetherill,exterr=exterr,d=x$d)
        if (exterr){ # show decay constant uncertainty
            if (i > 1) oldell <- ell
            ell <- ellipse(xy$x[1],xy$x[2],xy$cov,alpha=alpha)
            if (i > 1){
                xycd <- rbind(oldell,ell)
                ii <- grDevices::chull(xycd)
                graphics::polygon(xycd[ii,],col=col,border=NA)
            }
        }
        conc[i,] <- xy$x
    }
    if (length(ticks)<2)
        ticks <- prettier(lims$t,wetherill=wetherill,n=ticks)
    graphics::lines(conc[,'x'],conc[,'y'],col=col,lwd=2)
    for (i in 1:length(ticks)){
        xy <- age_to_concordia_ratios(ticks[i],wetherill=wetherill,exterr=exterr,d=x$d)
        if (exterr){ # show ticks as ellipse
            ell <- ellipse(xy$x[1],xy$x[2],xy$cov,alpha=alpha)
            graphics::polygon(ell,col='white')
        } else {
            graphics::points(xy$x[1],xy$x[2],pch=21,bg='white')
        }
        pos <- 2
        if ((wetherill  & diff(range(conc[,'x']))<0.05) |
            (!wetherill & diff(range(conc[,'x']))<2.5) & exterr){ pos <- NULL }
        graphics::text(xy$x[1],xy$x[2],as.character(ticks[i]),pos=pos)
    }
    graphics::box()
}
# helper function for plot.concordia
prepare.concordia.line <- function(x,tlim,wetherill=TRUE,...){
    lims <- get.concordia.limits(x,tlim=tlim,wetherill=wetherill,...)
    if (wetherill){
        x.lab <- expression(paste(""^"207","Pb/"^"235","U"))
        y.lab <- expression(paste(""^"206","Pb/"^"238","U"))
    } else {
        x.lab <- expression(paste(""^"238","U/"^"206","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"206","Pb"))
    }
    graphics::plot(lims$x,lims$y,type='n',xlab=x.lab,ylab=y.lab,bty='n',...)
    lims
}
cseq <- function(m,M,wetherill=TRUE,n=50){
    if (wetherill){
        return(seq(m,M,length.out=n))
    }
    if (m>0){
        out <- exp(seq(0,log(M/m),length.out=n))*m
    } else {
        out <- exp(seq(0,log(M+1-m),length.out=n))+1-m
    }
    out[out<=0] <- min(out[out>0])/10
    out
}
prettier <- function(x,wetherill=TRUE,n=5){
    pilot <- pretty(x,n=n)
    if (wetherill){
        return(pilot)
    }
    m <- min(x)
    M <- max(x)
    if (M/m < 50){ # linear spacing if spanning less than 1 order of magnitude
        out <- pilot
    } else { # log spaced if spanning more than 1 order of magnitude
        out <- cseq(m,M,wetherill=wetherill,n=n)
        init <- out
        out[1] <- pilot[1]
        out[n] <- utils::tail(pilot,1)
        for (i in 2:(n-1)){ # prettify
            out[i] <- pretty(init[(i-1):i],n=2)[2]
        }
    }
    out[out<=0] <- min(out[out>0])/2
    out
}
age_to_concordia_ratios <- function(tt,wetherill=TRUE,exterr=FALSE,d=diseq()){
    if (wetherill) return(age_to_wetherill_ratios(tt,exterr=exterr,d=d))
    else return(age_to_terawasserburg_ratios(tt,exterr=exterr,d=d))
}
get.concordia.limits <- function(x,tlim=NULL,wetherill=FALSE,...){
    out <- list()
    args <- list(...)
    xset <- ('xlim' %in% names(args))
    yset <- ('ylim' %in% names(args))
    if (xset) {
        out$x <- args$xlim
        minx <- args$xlim[1]
        maxx <- args$xlim[2]
    } else {
        out$x <- c(0,0)
    }
    if (yset) {
        out$y <- args$ylim
        miny <- args$ylim[1]
        maxy <- args$ylim[2]
    } else {
        out$y <- c(0,0)
    }
    if (is.null(tlim)) out$t <- c(0,0)
    else out$t <- tlim
    nse <- 3 # number of standard errors used for buffer
    if (!is.null(tlim) && wetherill){
        if (!xset) out$x <- age_to_Pb207U235_ratio(tlim,d=x$d)[,'75']
        if (!yset) out$y <- age_to_Pb206U238_ratio(tlim,d=x$d)[,'68']
    } else if (!is.null(tlim) && !wetherill){
        if (tlim[1] <= 0){
            U238Pb206 <- get.U238Pb206.ratios(x)
            if (xset) maxx <- out$x[2]
            else maxx <- max(U238Pb206[,1]+nse*U238Pb206[,2],na.rm=TRUE)
            out$t[1] <- get.Pb206U238.age(1/maxx,d=x$d)[1]
        }
        if (!xset) out$x <- age_to_U238Pb206_ratio(out$t,d=x$d)[,'86']
        if (!yset) out$y <- age_to_Pb207Pb206_ratio(out$t,d=x$d)[,'76']
    } else if (is.null(tlim) && wetherill) {
        if (!xset){
            Pb207U235 <- get.Pb207U235.ratios(x)
            minx <- min(Pb207U235[,1]-nse*Pb207U235[,2],na.rm=TRUE)
            maxx <- max(Pb207U235[,1]+nse*Pb207U235[,2],na.rm=TRUE)
        }
        if (!yset){
            Pb206U238 <- get.Pb206U238.ratios(x)
            miny <- min(Pb206U238[,1]-nse*Pb206U238[,2],na.rm=TRUE)
            maxy <- max(Pb206U238[,1]+nse*Pb206U238[,2],na.rm=TRUE)
        }
        out$t[1] <- get.Pb206U238.age(miny,d=x$d)[1]
        out$t[2] <- get.Pb207U235.age(maxx,d=x$d)[1]
        if (!xset){
            minx <- min(minx,age_to_Pb207U235_ratio(out$t[1],d=x$d)[,'75'])
            maxx <- max(maxx,age_to_Pb207U235_ratio(out$t[2],d=x$d)[,'75'])
        }
        if (!yset){
            miny <- min(miny,age_to_Pb206U238_ratio(out$t[1],d=x$d)[,'68'])
            maxy <- max(maxy,age_to_Pb206U238_ratio(out$t[2],d=x$d)[,'68'])
        }
        out$x <- c(minx,maxx)
        out$y <- c(miny,maxy)
    } else if (is.null(tlim) && !wetherill){
        U238Pb206 <- get.U238Pb206.ratios(x)
        Pb207Pb206 <- get.Pb207Pb206.ratios(x)
        if (!xset){
            minx <- min(U238Pb206[,1]-nse*U238Pb206[,2],na.rm=TRUE)
            maxx <- max(U238Pb206[,1]+nse*U238Pb206[,2],na.rm=TRUE)
        }
        if (!yset){
            miny <- min(Pb207Pb206[,1]-nse*Pb207Pb206[,2],na.rm=TRUE)
            maxy <- max(Pb207Pb206[,1]+nse*Pb207Pb206[,2],na.rm=TRUE)
        }
        out$t[1] <- get.Pb206U238.age(1/maxx,d=x$d)[1]
        out$t[2] <- get.Pb207Pb206.age(maxy,d=x$d)[1]
        if (!xset) minx <- min(minx,age_to_U238Pb206_ratio(out$t[2],d=x$d)[,'86'])
        if (!yset) miny <- min(miny,age_to_Pb207Pb206_ratio(out$t[1],d=x$d)[,'76'])
        out$x <- c(minx,maxx)
        out$y <- c(miny,maxy)
    }
    out
}

# this would be much easier in unicode but that doesn't render in PDF:
concordia.title <- function(fit,sigdig=2,alpha=0.05,...){
    rounded.age <- roundit(fit$age[1],fit$age[-1],sigdig=sigdig)
    expr1 <- expression('concordia age ='~a%+-%b~'|'~c~'Ma (n='*n*')')
    list1 <- list(a=rounded.age[1],b=rounded.age[2],c=rounded.age[3],n=fit$n)
    if (fit$mswd['combined']>1){
        expr1 <- expression('concordia age ='~a%+-%b~'|'~c~'|'~d~'Ma (n='*n*')')
        list1$d <- rounded.age[4]
    }
    line1 <- do.call('substitute',list(eval(expr1),list1))
    line2 <- substitute('MSWD ='~a~'|'~b~'|'~c~
                            ', p('*chi^2*') ='~d~'|'~e~'|'~f,
                        list(a=signif(fit$mswd['concordance'],2),
                             b=signif(fit$mswd['equivalence'],2),
                             c=signif(fit$mswd['combined'],2),
                             d=signif(fit$p.value['concordance'],2),
                             e=signif(fit$p.value['equivalence'],2),
                             f=signif(fit$p.value['combined'],2)))
    mymtext(line1,line=1,...)
    mymtext(line2,line=0,...)
}

concordia.age <- function(x,...){ UseMethod("concordia.age",x) }
concordia.age.default <- function(x,...){
    stop("no default method implemented for concordia.age()")
}
concordia.age.UPb <- function(x,i=NA,wetherill=TRUE,exterr=TRUE,alpha=0.05,...){
    if (is.na(i)){
        ccw <- concordia.comp(x,wetherill=TRUE)
        cct <- concordia.comp(x,wetherill=FALSE)
    } else {
        ccw <- wetherill(x,i)
        cct <- tera.wasserburg(x,i)
    }
    tt <- concordia.age(ccw,d=x$d,exterr=exterr)
    out <- list()
    if (is.na(i)){ # these calculations are only relevant to weighted means
        out <- c(out,mswd.concordia(x,ccw,tt[1],exterr=exterr))
        out$age <- rep(NA,4)
        names(out$age) <- c('t','s[t]','ci[t]','disp[t]')
        tfact <- stats::qt(1-alpha/2,out$df['combined'])
        out$age[c('t','s[t]')] <- tt
        out$age['ci[t]'] <- tfact*out$age['s[t]']
        if (out$mswd['combined']>1)
            out$age['disp[t]'] <- tfact*out$mswd['combined']*out$age['s[t]']
        if (wetherill) cc <- ccw
        else cc <- cct
        out$x <- cc$x
        out$cov <- cc$cov
    } else {
        out$age <- tt
    }
    out
}
concordia.age.wetherill <- function(x,d=diseq(),exterr=FALSE,...){
    fit <- stats::optimise(LL.concordia.age,interval=c(0,5000),
                           exterr=FALSE,ccw=x,d=d)
    tt <- fit$minimum
    hess <- stats::optimHess(tt,fn=LL.concordia.age,
                             ccw=x,exterr=exterr,d=d)
    if (det(hess)>1e-15)
        tt.err <- as.numeric(sqrt(solve(hess)))
    else
        tt.err <- .Machine$double.ulp.digits
    out <- c(tt,tt.err)
    names(out) <- c('t.conc','s[t.conc]')
    out
}

# x has class 'UPb'
concordia.comp <- function(x,wetherill=TRUE){
    X <- flat.UPb.table(x,wetherill=wetherill)
    out <- wtdmean2D(X)
    cnames <- colnames(X)[c(1,3)]
    names(out$x) <- cnames
    colnames(out$cov) <- cnames
    if (wetherill) class(out) <- 'wetherill'
    else class(out) <- 'terawasserburg'
    out
}

# x is an object of class "terawasserburg"
initial.concordia.age <- function(x,d=diseq()){
    e <- eigen(x$cov)
    v <- e$vectors[,1]
    if (v[1]==0) return(get.Pb207Pb206.age(x$x['Pb207Pb206'],0,d=d)[1])
    if (v[2]==0) return(get.Pb206U238.age(1/x$x['U238Pb206'],0,d=d)[1])
    b <- v[1]/v[2] # slope of the ellipse
    x0 <- x$x['U238Pb206']
    y0 <- x$x['Pb207Pb206']
    a <- y0 - b*x0
    fit <- concordia.intersection.ab(a,b,exterr=FALSE,d=d)
    fit[1]
}

mswd.concordia <- function(x,ccw,tt,exterr=TRUE){
    SS.equivalence <- 
        LL.concordia.comp(mu=ccw$x,x=x,wetherill=TRUE,mswd=TRUE)
    SS.concordance <- 
        LL.concordia.age(tt=tt,ccw=ccw,mswd=TRUE,exterr=exterr,d=x$d)
    df.equivalence <- 2*length(x)-2
    df.concordance <- 1
    mswd <- rep(0,3)
    p.value <- rep(0,3)
    df <- rep(0,3)
    labels <- c('equivalence','concordance','combined')
    names(mswd) <- labels
    names(p.value) <- labels
    names(df) <- labels
    mswd['equivalence'] <-
        SS.equivalence/df.equivalence
    mswd['concordance'] <-
        SS.concordance/df.concordance
    mswd['combined'] <-
        (SS.equivalence+SS.concordance)/(df.equivalence+df.concordance)
    p.value['equivalence'] <-
        1-stats::pchisq(SS.equivalence,df.equivalence)
    p.value['concordance'] <-
        1-stats::pchisq(SS.concordance,df.concordance)
    p.value['combined'] <-
        1-stats::pchisq(SS.equivalence+SS.concordance,
                        df.equivalence+df.concordance)
    df['equivalence'] <- df.equivalence
    df['concordance'] <- df.concordance
    df['combined'] <- df.equivalence + df.concordance
    list(mswd=mswd,p.value=p.value,df=df)
}

LL.concordia.comp <- function(mu,x,wetherill=TRUE,mswd=FALSE,...){
    out <- 0
    for (i in 1:length(x)){
        if (wetherill) xi <- wetherill(x,i)
        else xi <- tera.wasserburg(x,i)
        X <- matrix(xi$x[1:2]-mu[1:2],1,2)
        covmat <- xi$cov[1:2,1:2]
        if (mswd) out <- out + get.concordia.SS(X,covmat)
        else out <- out + LL.norm(X,covmat)
    }
    out
}

LL.concordia.age <- function(tt,ccw,mswd=FALSE,exterr=TRUE,d=diseq()){
    y <- age_to_wetherill_ratios(tt,d=d)
    dx <- matrix(ccw$x[1:2]-y$x,1,2)
    covmat <- ccw$cov[1:2,1:2]
    if (exterr){
        l8 <- settings('lambda','U238')
        l5 <- settings('lambda','U235')
        P235 <- tt*exp(l5[1]*tt)
        P238 <- tt*exp(l8[1]*tt)
        E <- diag(c(P235*l5[2],P238*l8[2]))^2
        covmat <- ccw$cov[1:2,1:2] + E
    }
    if (mswd) out <- get.concordia.SS(dx,covmat)
    else out <- LL.norm(dx,covmat)
    out
}

get.concordia.SS <- function(x,covmat){
    x %*% solve(covmat) %*% t(x)
}

emptyconcordia <- function(tlim=NULL,alpha=0.05,wetherill=TRUE,exterr=TRUE,
                           concordia.col='darksalmon',ticks=5,...){
    if (is.null(tlim) && wetherill) tlim <- c(1,3500)
    else if (is.null(tlim)) tlim <- c(100,3500)
    dat <- list()
    class(dat) <- 'UPb'
    dat$d <- diseq()
    if (wetherill){
        dat$x <- rbind(c(age_to_Pb207U235_ratio(tlim[1]),age_to_Pb206U238_ratio(tlim[1]),0),
                       c(age_to_Pb207U235_ratio(tlim[1]),age_to_Pb206U238_ratio(tlim[1]),0))
        dat$format <- 1
    } else {
        dat$x <- rbind(c(age_to_U238Pb206_ratio(tlim[1]),age_to_Pb207Pb206_ratio(tlim[1])),
                       c(age_to_U238Pb206_ratio(tlim[1]),age_to_Pb207Pb206_ratio(tlim[1])))
        dat$format <- 2
    }
    lims <- prepare.concordia.line(x=dat,tlim=tlim,wetherill=wetherill,...)
    plot.concordia.line(x=dat,lims=lims,wetherill=wetherill,
                        col=concordia.col,alpha=alpha,exterr=exterr,ticks=ticks)
}
