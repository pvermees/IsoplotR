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
#' @param x an object of class \code{UPb}
#' @param tlim age limits of the concordia line
#' @param alpha confidence cutoff for the error ellipses
#' @param wetherill logical flag (\code{FALSE} for Tera-Wasserburg)
#' @param show.numbers logical flag (\code{TRUE} to show grain
#'     numbers)
#' @param levels a vector with additional values to be displayed as
#'     different background colours within the error ellipses.
#' @param ellipse.col a vector of two background colours for the error
#'     ellipses. If \code{levels=NA}, then only the first colour will
#'     be used. If \code{levels} is a vector of numbers, then
#'     \code{ellipse.col} is used to construct a colour ramp.
#' @param concordia.col colour of the concordia line
#' @param exterr show decay constant uncertainty?
#' @param show.age one of either:
#'
#' \code{0}: just plot the data but don't calculate the age
#'
#' \code{1}: calculate the concordia age
#'
#' \code{2}: fit a discordia line
#'
#' @param sigdig number of significant digits for the
#'     concordia/discordia age
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
#' @param ticks an optional vector of age ticks to be added to the
#'     concordia line.
#' @param ... optional arguments to the generic \code{plot} function
#' @examples
#' data(examples) 
#' concordia(examples$UPb)
#' @references Ludwig, K.R., 1998. On the treatment of concordant
#'     uranium-lead ages. Geochimica et Cosmochimica Acta, 62(4),
#'     pp.665-676.
#' @export
concordia <- function(x,tlim=NULL,alpha=0.05,wetherill=TRUE,
                      show.numbers=FALSE,levels=NA,
                      ellipse.col=c("#00FF0080","#FF000080"),
                      concordia.col='darksalmon',exterr=TRUE,
                      show.age=0,sigdig=2,common.Pb=0,ticks=NULL,...){
    if (common.Pb>0) X <- common.Pb.correction(x,option=common.Pb)
    else X <- x
    concordia.line(X,tlim=tlim,wetherill=wetherill,col=concordia.col,
                   alpha=alpha,exterr=exterr,ticks=ticks,...)
    if (show.age==2){
        fit <- concordia.intersection.ludwig(x,wetherill=wetherill,exterr=exterr)
        discordia.plot(fit,wetherill=wetherill)
        graphics::title(discordia.title(fit,wetherill=wetherill,sigdig=sigdig))
    }
    ns <- length(x)
    ellipse.cols <- set.ellipse.colours(ns=ns,levels=levels,col=ellipse.col)
    for (i in 1:ns){
        if (wetherill) xyc <- wetherill(X,i)
        else xyc <- tera.wasserburg(X,i)
        x0 <- xyc$x[1]
        y0 <- xyc$x[2]
        covmat <- xyc$cov
        ell <- ellipse(x0,y0,covmat,alpha=alpha)
        graphics::polygon(ell,col=ellipse.cols[i])
        if (show.numbers) graphics::text(x0,y0,i)
        else graphics::points(x0,y0,pch=19,cex=0.25)
    }
    colourbar(z=levels,col=ellipse.col)
    if (show.age==1){
        fit <- concordia.age(X,wetherill=wetherill,exterr=exterr)
        ell <- ellipse(fit$x[1],fit$x[2],fit$cov)
        graphics::polygon(ell,col='white')
        graphics::title(concordia.title(fit,sigdig=sigdig))
    }
}

# helper function for plot.concordia
concordia.line <- function(x,tlim,wetherill,col,alpha=0.05,
                           exterr=TRUE,ticks=NULL,...){
    lims <- get.concordia.limits(x,tlim=tlim,wetherill=wetherill,...)
    if (wetherill){
        x.lab <- expression(paste(""^"207","Pb/"^"235","U"))
        y.lab <- expression(paste(""^"206","Pb/"^"238","U"))
    } else {
        x.lab <- expression(paste(""^"238","U/"^"206","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"206","Pb"))
    }
    graphics::plot(lims$x,lims$y,type='n',xlab=x.lab,ylab=y.lab,...)
    range.t <- range(lims$t)
    m <- max(0.8*lims$t[1],lims$t[1]-range.t/20)
    M <- min(1.2*lims$t[2],lims$t[2]+range.t/20)
    nn <- 30 # number of segments into which the concordia line is divided
    tt <- prettier(c(m,M),wetherill=wetherill,n=nn)
    concordia <- matrix(0,nn,2)
    colnames(concordia) <- c('x','y')
    for (i in 1:nn){ # build the concordia line
        xy <- age_to_concordia_ratios(tt[i],wetherill=wetherill,exterr=exterr)
        if (exterr){ # show decay constant uncertainty
            if (i > 1) oldell <- ell
            ell <- ellipse(xy$x[1],xy$x[2],xy$cov,alpha=alpha)
            if (i > 1){
                xycd <- rbind(oldell,ell)
                ii <- grDevices::chull(xycd)
                graphics::polygon(xycd[ii,],col=col,border=NA)
            }
        }
        concordia[i,] <- xy$x
    }
    graphics::lines(concordia[,'x'],concordia[,'y'],col=col,lwd=2)
    # prepare and plot ticks
    if (is.null(ticks)) ticks <- pretty(tt)
    for (i in 1:length(ticks)){
        xy <- age_to_concordia_ratios(ticks[i],wetherill=wetherill,exterr=exterr)
        if (exterr){ # show ticks as ellipse
            ell <- ellipse(xy$x[1],xy$x[2],xy$cov,alpha=alpha)
            graphics::polygon(ell,col='white')
        } else {
            graphics::points(xy$x[1],xy$x[2],pch=21,bg='white')
        }
        pos <- 2
        if ((wetherill  & diff(range(concordia[,'x']))<0.05) |
            (!wetherill & diff(range(concordia[,'x']))<2.5) & exterr){ pos <- NULL }
        graphics::text(xy$x[1],xy$x[2],as.character(ticks[i]),pos=pos)
    }
}
prettier <- function(x,wetherill=TRUE,n=20){
    m <- min(x)
    M <- max(x)
    out <- pretty(x,n=n,min.n=n)
    ntocull <- length(out)-n
    icull <- seq(from=2,to=2*ntocull,by=2)
    out <- out[-icull]
    if (wetherill){
        out[out<m] <- m
        out[out<0] <- 0
    } else {
        out[1] <- m
        out[out<=0] <- 1e-10
    }
    out[n] <- M
    out
}
age_to_concordia_ratios <- function(tt,wetherill=TRUE,exterr=FALSE){
    if (wetherill) return(age_to_wetherill_ratios(tt,exterr=exterr))
    else return(age_to_terawasserburg_ratios(tt,exterr=exterr))
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
        if (!xset) out$x <- age_to_Pb207U235_ratio(tlim)[,'75']
        if (!yset) out$y <- age_to_Pb206U238_ratio(tlim)[,'68']
    } else if (!is.null(tlim) && !wetherill){
        if (tlim[1] <= 0){
            U238Pb206 <- get.U238Pb206.ratios(x)
            if (xset) maxx <- out$x[2]
            else maxx <- max(U238Pb206[,1]+nse*U238Pb206[,2],na.rm=TRUE)
            out$t[1] <- get.Pb206U238.age(1/maxx)[1]
        }
        if (!xset) out$x <- age_to_U238Pb206_ratio(out$t)[,'86']
        if (!yset) out$y <- age_to_Pb207Pb206_ratio(out$t)[,'76']
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
        out$t[1] <- get.Pb206U238.age(miny)[1]
        out$t[2] <- get.Pb207U235.age(maxx)[1]
        if (!xset){
            minx <- min(minx,age_to_Pb207U235_ratio(out$t[1])[,'75'])
            maxx <- max(maxx,age_to_Pb207U235_ratio(out$t[2])[,'75'])
        }
        if (!yset){
            miny <- min(miny,age_to_Pb206U238_ratio(out$t[1])[,'68'])
            maxy <- max(maxy,age_to_Pb206U238_ratio(out$t[2])[,'68'])
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
        out$t[1] <- get.Pb206U238.age(1/maxx)[1]
        out$t[2] <- get.Pb207Pb206.age(maxy)[1]
        if (!xset) minx <- min(minx,age_to_U238Pb206_ratio(out$t[2])[,'86'])
        if (!yset) miny <- min(miny,age_to_Pb207Pb206_ratio(out$t[1])[,'76'])
        out$x <- c(minx,maxx)
        out$y <- c(miny,maxy)
    }
    out
}

concordia.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age[1],fit$age[2],sigdig=sigdig)
    line1 <- substitute('concordia age ='~a%+-%b~'[Ma] (1'~sigma~')',
                        list(a=rounded.age[1], b=rounded.age[2]))
    line2 <- substitute('MSWD (concordance) ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd$concordance,2),
                             b=signif(fit$p.value$concordance,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}

concordia.age <- function(x,...){ UseMethod("concordia.age",x) }
concordia.age.default <- function(x,...){
    stop("no default method implemented for concordia.age()")
}
concordia.age.UPb <- function(x,i=NA,wetherill=TRUE,exterr=TRUE,...){
    if (is.na(i)){
        ccw <- concordia.comp(x,wetherill=TRUE)
        cct <- concordia.comp(x,wetherill=FALSE)
    } else {
        ccw <- wetherill(x,i)
        cct <- tera.wasserburg(x,i)
    }
    t.init <- initial.concordia.age(cct)
    out <- concordia.age(ccw,t.init=t.init,exterr=exterr)
    mswd <- mswd.concordia(x,ccw,out$age[1],exterr=exterr)
    if (wetherill) cc <- ccw
    else cc <- cct
    out$x <- cc$x
    out$cov <- cc$cov
    out$mswd <- mswd$mswd
    out$p.value <- mswd$p.value
    out
}
concordia.age.wetherill <- function(x,t.init,exterr=TRUE,...){
    out <- list()
    out$age <- tryCatch({
        fit <- stats::optim(par=t.init,fn=LL.concordia.age,
                            ccw=x,exterr=exterr,
                            method="BFGS",hessian=TRUE)
        tt <- fit$par
        tt.err <- as.numeric(sqrt(solve(fit$hessian)))
        c(tt,tt.err)
    }, warning = function(w){
        c(t.init,1)
    }, error = function(e){
        c(t.init,1)
    })
    names(out$age) <- c('t.conc','s[t.conc]')
    out
}

# x has class 'UPb'
concordia.comp <- function(x,wetherill=TRUE){
    if (wetherill) out <- wetherill(x,1)
    else out <- tera.wasserburg(x,1)
    fit.comp <- stats::optim(out$x[1:2], LL.concordia.comp,
                             x=x,wetherill=wetherill,
                             method="BFGS",hessian=TRUE)
    out$x <- fit.comp$par
    out$cov <- solve(fit.comp$hessian)
    out
}

# x is an object of class "terawasserburg"
initial.concordia.age <- function(x){
    e <- eigen(x$cov)
    v <- e$vectors[,1]
    if (v[1]==0) return(get.Pb207Pb206.age(x$x['Pb207Pb206'],0)[1])
    if (v[2]==0) return(get.Pb206U238.age(1/x$x['U238Pb206'],0)[1])
    b <- v[1]/v[2] # slope of the ellipse
    x0 <- x$x['U238Pb206']
    y0 <- x$x['Pb207Pb206']
    a <- y0 - b*x0
    fit <- concordia.intersection.york.ab(a,b,exterr=FALSE)
    fit$x[1]
}

mswd.concordia <- function(x,ccw,tt,exterr=TRUE){
    out <- list()
    SS.equivalence <- LL.concordia.comp(mu=ccw$x,x=x,wetherill=TRUE,mswd=TRUE)
    SS.concordance <- LL.concordia.age(tt=tt,ccw=ccw,mswd=TRUE,exterr=exterr)
    df.equivalence <- 2*length(x)-2
    df.concordance <- 1
    out$mswd <- list(equivalence = SS.equivalence/df.equivalence,
                     concordance = SS.concordance/df.concordance)
    out$p.value <- list(equivalence = 1-stats::pchisq(SS.equivalence,df.equivalence),
                        concordance = 1-stats::pchisq(SS.concordance,df.concordance))
    out
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

LL.concordia.age <- function(tt,ccw,mswd=FALSE,exterr=TRUE){
    y <- age_to_wetherill_ratios(tt)
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
