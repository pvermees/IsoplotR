#' @title Concordia diagram
#'
#' @description
#' Plots U-Pb data on Wetherill, Tera-Wasserburg or U-Th-Pb concordia
#' diagrams, calculates concordia ages and compositions, evaluates the
#' equivalence of multiple
#' (\eqn{^{206}}Pb/\eqn{^{238}}U-\eqn{^{207}}Pb/\eqn{^{235}}U,
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U, or
#' \eqn{^{208}}Th/\eqn{^{232}}Th-\eqn{^{206}}Pb/\eqn{^{238}}U)
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
#' diagram); or, equivalently, the \eqn{^{207}}Pb/\eqn{^{206}}Pb- and
#' \eqn{^{206}}Pb/\eqn{^{238}}U-ratios (`Tera-Wasserburg'
#' diagram). Alternatively, for data format 7 and 8, it is also
#' possible to plot \eqn{^{208}}Pb/\eqn{^{232}}Th against the
#' \eqn{^{206}}Pb/\eqn{^{238}}U.  The space of concordant isotopic
#' compositions is marked by a curve, the `concordia line'. Isotopic
#' ratio measurements are shown as 100(1-\code{alpha})\% confidence
#' ellipses. Concordant samples plot near to, or overlap with, the
#' concordia line. They represent the pinnacle of geochronological
#' robustness. Samples that plot away from the concordia line but are
#' aligned along a linear trend form an isochron (or `discordia' line)
#' that can be used to infer the composition of the non-radiogenic
#' (`common') lead or to constrain the timing of prior lead loss.
#'
#' @param x an object of class \code{UPb}
#' 
#' @param tlim age limits of the concordia line
#' 
#' @param alpha probability cutoff for the error ellipses and
#'     confidence intervals
#' 
#' @param type one of
#'
#' \code{1}: Wetherill -- \eqn{{}^{206}}Pb/\eqn{{}^{238}}U
#' vs. \eqn{{}^{207}}Pb/\eqn{{}^{235}}U
#'
#' \code{2}: Tera-Wasserburg -- \eqn{{}^{207}}Pb/\eqn{{}^{206}}Pb
#' vs. \eqn{{}^{238}}U/\eqn{{}^{206}}Pb
#'
#' \code{3}: U-Th-Pb concordia -- \eqn{{}^{208}}Pb/\eqn{{}^{232}}Th
#' vs. \eqn{{}^{206}}Pb/\eqn{{}^{238}}U (only available if
#' \code{x$format=7} or \code{8})
#' 
#' @param show.numbers logical flag (\code{TRUE} to show grain
#'     numbers)
#' 
#' @param levels a vector with \code{length(x)} values to be displayed
#'     as different background colours within the error ellipses.
#' 
#' @param clabel label for the colour legend (only used if
#'     \code{levels} is not \code{NA}).
#' 
#' @param ellipse.fill
#' Fill colour for the error ellipses. This can either be a single
#' colour or multiple colours to form a colour ramp. Examples:
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
#' 
#' @param concordia.col colour of the concordia line
#' 
#' @param exterr show decay constant uncertainties?
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
#' @param common.Pb common lead correction:
#'
#' \code{0}:none
#'
#' \code{1}: use the Pb-composition stored in
#' 
#' \code{settings('iratio','Pb207Pb206')} (if \code{x$format<4});
#' 
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')} (if \code{3<x$format<7}); or
#' 
#' \code{settings('iratio','Pb208Pb206')} and
#' \code{settings('iratio','Pb208Pb207')} (if \code{x$format>6}).
#' 
#' \code{2}: use the isochron intercept as the initial Pb-composition
#'
#' \code{3}: use the Stacey-Kramers two-stage model to infer the initial
#' Pb-composition.
#'
#' @param anchor
#' control parameters to fix the intercept age or common Pb
#' composition of the isochron fit. This can be a scalar or a vector.
#'
#' If \code{anchor[1]=0}: do not anchor the isochron.
#'
#' If \code{anchor[1]=1}: fix the common Pb composition at the values
#' stored in \code{settings('iratio',...)}.
#'
#' If \code{anchor[1]=2}: force the isochron line to intersect the
#' concordia line at an age equal to \code{anchor[2]}.
#'
#' @param ticks either a scalar indicating the desired number of age
#'     ticks to be placed along the concordia line, OR a vector of
#'     tick ages.
#' @param hide vector with indices of aliquots that should be removed
#'     from the concordia diagram
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from concordia or discordia age calculation
#' @param omit.fill fill colour that should be used for the omitted
#'     aliquots.
#' @param omit.stroke stroke colour that should be used for the
#'     omitted aliquots.
#' @param ... optional arguments to the generic \code{plot} function
#'
#' @return
#'
#' If \code{show.age=1}, returns a list with the following items:
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
#' freedom used for the \code{mswd} calculation. }
#'
#' \item{age}{a 4-element vector with:\cr
#' 
#' \code{t}: the concordia age (in Ma)\cr
#' 
#' \code{s[t]}: the estimated uncertainty of \code{t}\cr
#' 
#' \code{ci[t]}: the \eqn{100(1-\alpha)\%} confidence interval of
#' 
#' \code{t} for the appropriate degrees of freedom\cr
#' 
#' \code{disp[t]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{t} augmented by \eqn{\sqrt{mswd}} to account for
#' overdispersed datasets.
#' }
#'
#' }
#'
#' If \code{show.age=2}, \code{3} or \code{4}, returns a list with the
#' following items:
#'
#' \describe{
#'
#' \item{model}{the fitting model (\code{=show.age-1}).}
#'
#' \item{par}{a vector with the upper and lower intercept ages (if
#' \code{type=1}) or the lower intercept age and common Pb
#' intercept(s) (if \code{type=2}). If \code{show.age=4}, includes an
#' overdispersion term as well.}
#'
#' \item{cov}{ the covariance matrix of the elements in \code{par}.}
#'
#' \item{logpar}{the logarithm of \code{par}}
#'
#' \item{logcov}{the logarithm of \code{cov}}
#'
#' \item{err}{ a matrix with the following rows:
#'
#' \code{s}: the estimated standard deviation for \code{x}
#'
#' \code{ci}: the \eqn{100(1-\alpha)\%} confidence interval of
#' \code{x} for the appropriate degrees of freedom
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
#' \item{n}{ the number of aliquots in the dataset }
#'
#' }
#'
#' @examples
#' attach(examples)
#' concordia(UPb,show.age=2)
#'
#' dev.new()
#' concordia(UPb,type=2,xlim=c(24.9,25.4),
#'           ylim=c(0.0508,0.0518),ticks=249:254,exterr=TRUE)
#'
#' dev.new()
#' concordia(UPb,show.age=2,anchor=c(2,260))
#'
#' @references Ludwig, K.R., 1998. On the treatment of concordant
#'     uranium-lead ages. Geochimica et Cosmochimica Acta, 62(4),
#'     pp.665-676.
#' @export
concordia <- function(x=NULL,tlim=NULL,alpha=0.05,type=1,
                      show.numbers=FALSE,levels=NA,clabel="",
                      ellipse.fill=c("#00FF0080","#FF000080"),
                      ellipse.stroke='black',
                      concordia.col='darksalmon',exterr=FALSE,
                      show.age=0,sigdig=2,common.Pb=0,
                      ticks=5,anchor=0,hide=NULL,omit=NULL,
                      omit.fill=NA,omit.stroke='grey',...){    
    if (is.null(x)){
        emptyconcordia(tlim=tlim,alpha=alpha,type=type,exterr=exterr,
                       concordia.col=concordia.col,ticks=ticks,...)
        return(invisible(NULL))
    }
    ns <- length(x)
    plotit <- (1:ns)%ni%hide
    calcit <- (1:ns)%ni%c(hide,omit)
    x2calc <- subset(x,subset=calcit)
    if (common.Pb<1) X <- x
    else X <- Pb0corr(x,option=common.Pb,omit4c=unique(c(hide,omit)))
    X2plot <- subset(X,subset=plotit)
    fit <- NULL
    lims <- prepare.concordia.line(x=X2plot,tlim=tlim,type=type,...)
    if (show.age>1){
        wetherill <- (type==1)
        fit <- concordia.intersection.ludwig(x2calc,wetherill=wetherill,
                                             exterr=exterr,alpha=alpha,
                                             model=(show.age-1),anchor=anchor)
        wethertit <- wetherill & !measured.disequilibrium(x2calc$d)
        if (measured.disequilibrium(x2calc$d))
            x2calc$d <- replace.impossible.diseq(tt=fit$par[1],d=x2calc$d)
        fit$n <- length(x2calc)
        discordia.line(fit,wetherill=wetherill,d=mediand(x2calc$d))
        graphics::title(discordia.title(fit,wetherill=wethertit,sigdig=sigdig))
    }
    plot.concordia.line(x2calc,lims=lims,type=type,col=concordia.col,
                        alpha=alpha,exterr=exterr,ticks=ticks)
    if (type==1) y <- data2york(X,option=1)
    else if (type==2) y <- data2york(X,option=2)
    else if (x$format%in%c(7,8) & type==3) y <- data2york(X,option=5)
    else stop('Concordia type incompatible with this input format.')
    scatterplot(y,alpha=alpha,show.numbers=show.numbers,
                show.ellipses=1*(show.age!=3),levels=levels,
                clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,add=TRUE,
                hide=hide,omit=omit,omit.fill=omit.fill,
                omit.stroke=omit.stroke,addcolourbar=FALSE,...)
    if (show.age==1){
        X2calc <- subset(X,subset=calcit)
        fit <- concordia.age(X2calc,type=type,exterr=exterr,alpha=alpha)
        ell <- ellipse(fit$x[1],fit$x[2],fit$cov)
        graphics::polygon(ell,col='white')
        fit$n <- length(X2calc)
        graphics::title(concordia.title(fit,sigdig=sigdig))
    }
    # must be added to the end because otherwise R doesn't
    # add the concordia ellipse to the scatterplot
    colourbar(z=levels[calcit],fill=ellipse.fill,
              stroke=ellipse.stroke,clabel=clabel)
    invisible(fit)
}

# helper function for plot.concordia
plot.concordia.line <- function(x,lims,type=1,col='darksalmon',
                                alpha=0.05,exterr=TRUE,ticks=5){
    if (length(ticks)<2)
        ticks <- prettier(lims$t,type=type,n=ticks)
    m <- min(lims$t[1],ticks[1])
    M <- max(lims$t[2],utils::tail(ticks,1))
    nn <- 30 # number of segments into which the concordia line is divided
    tt <- cseq(m,M,type=type,n=nn)
    conc <- matrix(0,nn,2)
    colnames(conc) <- c('x','y')
    for (i in 1:nn){ # build the concordia line
        xy <- age_to_concordia_ratios(tt[i],type=type,
                                      exterr=exterr,d=mediand(x$d))
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
    graphics::lines(conc[,'x'],conc[,'y'],col=col,lwd=2)
    dx <- diff(graphics::par('usr')[1:2])
    if (exterr & ((type==1 & dx<0.03) | (type==2 & dx<3) | (type==3 & dx<0.005)))
    { pos <- NULL } else { pos <- 2 }
    for (i in 1:length(ticks)){
        xy <- age_to_concordia_ratios(ticks[i],type=type,
                                      exterr=exterr,d=mediand(x$d))
        if (exterr){ # show ticks as ellipse
            ell <- ellipse(xy$x[1],xy$x[2],xy$cov,alpha=alpha)
            graphics::polygon(ell,col='white')
        } else {
            graphics::points(xy$x[1],xy$x[2],pch=21,bg='white')
        }
        graphics::text(xy$x[1],xy$x[2],as.character(ticks[i]),pos=pos)
    }
    graphics::box()
}
# helper function for plot.concordia
prepare.concordia.line <- function(x,tlim,type=1,...){
    out <- get.concordia.limits(x,tlim=tlim,type=type,...)
    if (type==1){
        y.lab <- expression(paste(""^"206","Pb/"^"238","U"))
        x.lab <- expression(paste(""^"207","Pb/"^"235","U"))
    } else if (type==2){
        x.lab <- expression(paste(""^"238","U/"^"206","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"206","Pb"))
    } else if (x$format>6 & type==3){
        x.lab <- expression(paste(""^"206","Pb/"^"238","U"))
        y.lab <- expression(paste(""^"208","Pb/"^"232","Th"))
    } else {
        stop('Incorrect input format.')
    }
    graphics::plot(out$x,out$y,type='n',xlab=x.lab,ylab=y.lab,bty='n',...)
    out
}
# concordia sequence
cseq <- function(m,M,type=1,n=50){
    if (type==1){
        out <- seq(m,M,length.out=n)
    } else if (m>0){
        out <- exp(seq(log(m),log(M),length.out=n))
    } else {
        out <- exp(seq(0,log(M+1),length.out=n))-1
    }
    out[out<=0] <- min(out[out>0])/10
    out
}
prettier <- function(x,type=1,n=5){
    if (type%in%c(1,3)){
        out <- pretty(x,n=n)
    } else {
        m <- min(x)
        M <- max(x)
        if (M/m<20){ # linear spacing if TW spans less than 1 order of magnitude
            out <- pretty(x,n=n)
            if (out[1]<=0) out[1] <- out[2]/10
        } else { # log spacing if TW spans more than 1 order of magnitude
            mexp <- log10(floor(10^(log10(m)%%1))) + floor(log10(m))
            Mexp <- log10(ceiling(10^(log10(M)%%1))) + floor(log10(M))
            out <- c(c(1,2,5) %o% 10^(mexp:Mexp))
        }
    }
    out
}
age_to_concordia_ratios <- function(tt,type=1,exterr=FALSE,d=diseq()){
    if (type==1)
        return(age_to_wetherill_ratios(tt,exterr=exterr,d=d))
    else if (type==2)
        return(age_to_terawasserburg_ratios(tt,exterr=exterr,d=d))
    else if (type==3)
        return(age_to_cottle_ratios(tt,exterr=exterr,d=d))
    else
        stop('Invalid concordia type.')
}
get.concordia.limits <- function(x,tlim=NULL,type=1,xlim,ylim,...){
    out <- list()
    if (missing(xlim)) {
        xset <- FALSE
        out$x <- c(NA,NA)
    } else {
        xset <- TRUE
        out$x <- xlim
        minx <- xlim[1]
        maxx <- xlim[2]        
    }
    if (missing(ylim)) {
        yset <- FALSE
        out$y <- c(NA,NA)
    } else {
        yset <- TRUE
        out$y <- ylim
        miny <- ylim[1]
        maxy <- ylim[2]
    }
    md <- mediand(x$d)
    nse <- 3 # number of standard errors used for buffer
    if (is.null(tlim)) out$t <- c(0,0)
    else out$t <- tlim
    if (measured.disequilibrium(x$d)){
        if (is.null(tlim)) out$t[2] <- meas.diseq.maxt(md)
        if (type==1){
            if (!xset){
                Pb7U5 <- get.Pb207U235.ratios(x)
                Pb7U5t <- age_to_Pb207U235_ratio(out$t,d=x$d)[,'75']
                minx <- min(Pb7U5[,1]-nse*Pb7U5[,2],Pb7U5t,na.rm=TRUE)
                maxx <- max(Pb7U5[,1]+nse*Pb7U5[,2],Pb7U5t,na.rm=TRUE)
            }
            Pb6U8t <- age_to_Pb206U238_ratio(out$t,d=x$d)[,'68']
            if (!yset){
                Pb6U8 <- get.Pb206U238.ratios(x)
                miny <- min(Pb6U8[,1]-nse*Pb6U8[,2],Pb6U8t,na.rm=TRUE)
                maxy <- max(Pb6U8[,1]+nse*Pb6U8[,2],na.rm=TRUE)
            }
            if (is.null(tlim) & maxy<Pb6U8t[2])
                out$t[2] <- get.Pb206U238.age(maxy,d=md)[1]
            out$x <- c(minx,maxx)
            out$y <- c(miny,maxy)
        } else if (type==2){
            if (is.null(tlim)) out$t[1] <- out$t[2]/2
            U8Pb6 <- get.U238Pb206.ratios(x)
            U8Pb6t <- age_to_U238Pb206_ratio(out$t,d=x$d)[,'86']
            if (!xset){
                minx <- min(U8Pb6[,1]-nse*U8Pb6[,2],U8Pb6t,na.rm=TRUE)
                maxx <- max(U8Pb6[,1]+nse*U8Pb6[,2],U8Pb6t,na.rm=TRUE)
            }
            if (is.null(tlim) & maxx>U8Pb6t[1])
                out$t[1] <- get.Pb206U238.age(1/maxx,d=md)[1]
            Pb76 <- get.Pb207Pb206.ratios(x)
            if (!yset){
                Pb76t <- age_to_Pb207Pb206_ratio(out$t,d=x$d)[,'76']
                miny <- min(Pb76[,1]-nse*Pb76[,2],Pb76t,na.rm=TRUE)
                maxy <- max(Pb76[,1]+nse*Pb76[,2],Pb76t,na.rm=TRUE)
            }
            out$x <- c(minx,maxx)
            out$y <- c(miny,maxy)
        } else if (type==3){
            Pb6U8t <- age_to_Pb206U238_ratio(out$t,d=x$d)[,'68']
            if (!xset){
                Pb6U8 <- get.Pb206U238.ratios(x)
                minx <- min(Pb6U8[,1]-nse*Pb6U8[,2],Pb6U8t,na.rm=TRUE)
                maxx <- max(Pb6U8[,1]+nse*Pb6U8[,2],na.rm=TRUE)
            }
            if (is.null(tlim) & maxx>Pb6U8t[1])
                out$t[2] <- get.Pb206U238.age(1/maxx,d=md)[1]
            if (!yset){
                Pb8Th2 <- get.Pb208Th232.ratios(x)
                Pb8Th2t <- age_to_Pb208Th232_ratio(out$t)[,'82']
                miny <- min(Pb8Th2[,1]-nse*Pb8Th2[,2],Pb8Th2t,na.rm=TRUE)
                maxy <- max(Pb8Th2[,1]+nse*Pb8Th2[,2],Pb8Th2t,na.rm=TRUE)
            }
            out$x <- c(minx,maxx)
            out$y <- c(miny,maxy)
        }
    } else {
        if (!is.null(tlim) & type==1){
            if (!xset) out$x <- age_to_Pb207U235_ratio(tlim,d=x$d)[,'75']
            if (!yset) out$y <- age_to_Pb206U238_ratio(tlim,d=x$d)[,'68']
        } else if (!is.null(tlim) & type==2){
            if (tlim[1] <= 0){
                U238Pb206 <- get.U238Pb206.ratios(x)
                if (xset) maxx <- out$x[2]
                else maxx <- max(U238Pb206[,1]+nse*U238Pb206[,2],na.rm=TRUE)
                out$t[1] <- get.Pb206U238.age(1/maxx,d=md)[1]
            }
            if (!xset) out$x <- age_to_U238Pb206_ratio(out$t,d=x$d)[,'86']
            if (!yset) out$y <- age_to_Pb207Pb206_ratio(out$t,d=x$d)[,'76']
        } else if (!is.null(tlim) & type==3){
            if (!xset) out$x <- age_to_Pb206U238_ratio(tlim,d=x$d)[,'68']
            if (!yset) out$y <- age_to_Pb208Th232_ratio(tlim)[,'82']
        } else if (is.null(tlim) & type==1) {
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
            out$t[1] <- get.Pb206U238.age(miny,d=md)[1]
            out$t[2] <- get.Pb207U235.age(maxx,d=md)[1]
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
        } else if (is.null(tlim) & type==2){
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
            out$t[1] <- get.Pb206U238.age(1/maxx,d=mediand(md))[1]
            out$t[2] <- get.Pb207Pb206.age(maxy,d=mediand(md),
                                           interval=c(out$t[1],10000))[1]
            if (!xset)
                minx <- min(minx,age_to_U238Pb206_ratio(out$t[2],d=md)[,'86'])
            if (!yset)
                miny <- min(miny,age_to_Pb207Pb206_ratio(out$t[1],d=md)[,'76'])
            out$x <- c(minx,maxx)
            out$y <- c(miny,maxy)
        } else if (is.null(tlim) & type==3){
            if (!xset){
                Pb206U238 <- get.Pb206U238.ratios(x)
                minx <- min(Pb206U238[,1]-nse*Pb206U238[,2],na.rm=TRUE)
                maxx <- max(Pb206U238[,1]+nse*Pb206U238[,2],na.rm=TRUE)
            }
            if (!yset){
                Pb208Th232 <- get.Pb208Th232.ratios(x)
                miny <- min(Pb208Th232[,1]-nse*Pb208Th232[,2],na.rm=TRUE)
                maxy <- max(Pb208Th232[,1]+nse*Pb208Th232[,2],na.rm=TRUE)
            }
            out$t[1] <- get.Pb206U238.age(minx,d=md)[1]
            out$t[2] <- get.Pb208Th232.age(maxy)[1]
            if (!xset){
                minx <- min(minx,age_to_Pb206U238_ratio(out$t[1],d=md)[,'68'])
                maxx <- max(maxx,age_to_Pb206U238_ratio(out$t[2],d=md)[,'68'])
            }
            if (!yset){
                miny <- min(miny,age_to_Pb208Th232_ratio(out$t[1])[,'82'])
                maxy <- max(maxy,age_to_Pb208Th232_ratio(out$t[2])[,'82'])
            }
            out$x <- c(minx,maxx)
            out$y <- c(miny,maxy)
        }
    }
    out
}

# this would be much easier in unicode but that doesn't render in PDF:
concordia.title <- function(fit,sigdig=2,alpha=0.05,...){
    rounded.age <- roundit(fit$age[1],fit$age[-1],sigdig=sigdig,text=TRUE)
    expr1 <- expression('concordia age ='~a%+-%b~'|'~c~'Ma (n='*n*')')
    list1 <- list(a=rounded.age[1],b=rounded.age[2],c=rounded.age[3],n=fit$n)
    if (fit$p.value['combined']<alpha){
        expr1 <- expression('concordia age ='~a%+-%b~'|'~c~'|'~d~'Ma (n='*n*')')
        list1$d <- rounded.age[4]
    }
    line1 <- do.call(substitute,list(eval(expr1),list1))
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

concordia.age <- function(x,i=NA,type=1,exterr=TRUE,alpha=0.05,...){
    if (is.na(i)){
        cc <- concordia.comp(x,type=type)
        if (type==3){
            cc4age <- cc
            type4age <- 3
        } else { # use Wetherill
            cc4age <- concordia.comp(x,type=1)
            type4age <- 1
        }
    } else {
        cc <- wetherill(x,i)
        cc4age <- cc
        type4age <- 1
    }
    tt <- concordia_age_helper(cc4age,d=mediand(x$d),type=type4age,exterr=exterr)
    out <- list()
    if (is.na(i)){ # these calculations are only relevant to weighted means
        out <- c(out,mswd.concordia(x,cc4age,type=type4age,tt=tt[1],exterr=exterr))
        out$age <- rep(NA,4)
        names(out$age) <- c('t','s[t]','ci[t]','disp[t]')
        out$age[c('t','s[t]')] <- tt
        out$age['ci[t]'] <- ntfact(alpha)*out$age['s[t]']
        mswd <- list(mswd=out$mswd['combined'],model=1,
                     p.value=out$p.value['combined'],
                     df=out$df['combined'],alpha=alpha)
        if (inflate(mswd)){
            out$age['disp[t]'] <- ntfact(alpha,mswd)*out$age['s[t]']
        }
        out$x <- cc$x
        out$cov <- cc$cov
    } else {
        out$age <- tt
    }
    out
}
concordia_age_helper <- function(cc,d=diseq(),type=1,exterr=FALSE,...){
    if (type==1) Pb206U238 <- cc$x['Pb206U238']
    else if (type==2) Pb206U238 <- 1/cc$x['U238Pb206']
    else if (type==3) Pb206U238 <- cc$x['Pb206U238']
    else stop('Incorrect concordia type.')
    midpoint <- get.Pb206U238.age(x=Pb206U238,d=d)[1]
    fit1 <- stats::optimise(LL.concordia.age,interval=c(0,midpoint),
                            exterr=FALSE,cc=cc,type=type,d=d)
    fit2 <- stats::optimise(LL.concordia.age,interval=c(midpoint,5000),
                            exterr=FALSE,cc=cc,type=type,d=d)
    if (fit1$objective<fit2$objective){
        tt <- fit1$minimum
    } else {
        tt <- fit2$minimum
    }
    hess <- stats::optimHess(tt,fn=LL.concordia.age,cc=cc,
                             exterr=exterr,type=type,d=d)
    if (det(hess)>1e-15)
        tt.err <- as.numeric(sqrt(solve(hess)))
    else
        tt.err <- .Machine$double.ulp.digits
    out <- c(tt,tt.err)
    names(out) <- c('t.conc','s[t.conc]')
    out
}

# x has class 'UPb'
concordia.comp <- function(x,type=1){
    if (type==1){
        X <- data2york(x,option=1)
        colnames(X) <- c('Pb207U235','errPb207U235',
                         'Pb206U238','errPb206U238','rhoXY')
    } else if (type==2){
        X <- data2york(x,option=2)
        colnames(X) <- c('U238Pb206','errU238Pb206',
                         'Pb207Pb206','errPb207Pb206','rhoXY')
    } else if (type==3){
        X <- data2york(x,option=5)
        colnames(X) <- c('Pb206U238','errPb206U238',
                         'Pb208Th232','errPb208Th232','rhoXY')
    } else {
        stop('Incorrect concordia type.')
    }
    out <- wtdmean2D(X)
    cnames <- colnames(X)[c(1,3)]
    names(out$x) <- cnames
    colnames(out$cov) <- cnames
    rownames(out$cov) <- cnames
    out
}

mswd.concordia <- function(x,cc,type=1,tt=0,exterr=TRUE){
    SS.equivalence <- 
        LL.concordia.comp(mu=cc$x,x=x,type=type,mswd=TRUE)
    SS.concordance <- 
        LL.concordia.age(tt=tt,cc=cc,type=type,exterr=exterr,
                         d=mediand(x$d),mswd=TRUE)
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

LL.concordia.comp <- function(mu,x,type=1,mswd=FALSE,...){
    out <- 0
    for (i in 1:length(x)){
        if (type==1){
            xi <- wetherill(x,i)
            j <- c(1,2)
        } else if (type==2){
            xi <- tera.wasserburg(x,i)
            j <- c(1,2)
        } else if (type==3){
            xi <- wetherill(x,i)
            j <- c(2,3)
        } else {
            stop('Incorrect concordia type.')
        }
        X <- matrix(xi$x[j]-mu,1,2)
        covmat <- xi$cov[j,j]
        if (mswd) out <- out + stats::mahalanobis(x=X,center=FALSE,cov=covmat)
        else out <- out + LL.norm(X,covmat)
    }
    out
}

LL.concordia.age <- function(tt,cc,type=1,exterr=TRUE,d=diseq(),mswd=FALSE){
    if (type==1){
        y <- age_to_wetherill_ratios(tt,d=d)
        cols <- c('Pb207U235','Pb206U238')
    } else if (type==2){
        y <- age_to_terawasserburg_ratios(tt,d=d)
        cols <- c('U238Pb206','Pb207Pb206')
    } else if (type==3){
        y <- age_to_cottle_ratios(tt,d=d)
        cols <- c('Pb206U238','Pb208Th232')
    } else {
        stop('Incorrect concordia type.')
    }
    dx <- matrix(cc$x[cols]-y$x[cols],1,2)
    covmat <- cc$cov[cols,cols]
    if (exterr){
        l5 <- settings('lambda','U235')
        l8 <- settings('lambda','U238')
        l2 <- settings('lambda','Th232')
        U <- settings('iratio','U238U235')
        Lcov <- diag(c(l5[2],l8[2],l2[2]))^2
        J <- matrix(0,2,3)
        D <- mclean(tt=tt,d=d,exterr=TRUE)
        if (type==1){
            J[1,1] <- D$dPb207U235dl35
            J[2,2] <- D$dPb206U238dl38
        } else if (type==2){
            J[1,2] <- -D$dPb206U238dl38/D$Pb206U238^2
            J[2,1] <- D$dPb207U235dl35/(U*D$Pb206U238)
            J[2,2] <- -D$dPb206U238dl38/(U*D$Pb206U238^2)
        } else { # type == 3
            J[1,2] <- D$dPb206U238dl38
            J[2,3] <- tt*exp(l2[1]*tt)
        }
        E <- J %*% Lcov %*% t(J)
        covmat <- covmat + E
    }
    if (mswd) out <- stats::mahalanobis(x=dx,center=FALSE,cov=covmat)
    else out <- LL.norm(dx,covmat)
    out
}

emptyconcordia <- function(tlim=NULL,alpha=0.05,type=1,exterr=TRUE,
                           concordia.col='darksalmon',ticks=5,...){
    if (is.null(tlim)){
        if (type%in%c(1,3)) tlim <- c(1,4500)
        else tlim <- c(100,4500)
    } 
    dat <- list()
    class(dat) <- 'UPb'
    dat$d <- diseq()
    if (type==1){
        dat$x <- rbind(c(age_to_Pb207U235_ratio(tlim[1]),0,
                         age_to_Pb206U238_ratio(tlim[1]),0,0),
                       c(age_to_Pb207U235_ratio(tlim[2]),0,
                         age_to_Pb206U238_ratio(tlim[2]),0,0))
        dat$format <- 1
    } else if (type==2){
        dat$x <- rbind(c(age_to_U238Pb206_ratio(tlim[1]),0,
                         age_to_Pb207Pb206_ratio(tlim[1]),0,0),
                       c(age_to_U238Pb206_ratio(tlim[2]),0,
                         age_to_Pb207Pb206_ratio(tlim[2]),0,0))
        dat$format <- 2
    } else if (type==3){
        dat$x <- rbind(c(0,0,age_to_Pb206U238_ratio(tlim[1]),0,
                         age_to_Pb208Th232_ratio(tlim[1]),0,rep(0,8)),
                       c(0,0,age_to_Pb206U238_ratio(tlim[2]),0,
                         age_to_Pb208Th232_ratio(tlim[2]),0,rep(0,8)))
        dat$format <- 7
    } else {
        stop('Invalid concordia type.')
    }
    lims <- prepare.concordia.line(x=dat,tlim=tlim,type=type,...)
    plot.concordia.line(x=dat,lims=lims,type=type,col=concordia.col,
                        alpha=alpha,exterr=exterr,ticks=ticks)
}
