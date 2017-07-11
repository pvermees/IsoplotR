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
#' @param limits age limits of the concordia line
#' @param alpha confidence cutoff for the error ellipses
#' @param wetherill logical flag (\code{FALSE} for Tera-Wasserburg)
#' @param show.numbers logical flag (\code{TRUE} to show grain numbers)
#' @param ellipse.col background colour of the error ellipses
#' @param concordia.col colour of the concordia line
#' @param exterr show decay constant uncertainty?
#' @param show.age one of either
#'
#' \code{1}: don't show the age
#'
#' \code{2}: calculate the concordia age
#'
#' \code{3}: fit a discordia line
#'
#' @param sigdig number of significant digits for the
#'     concordia/discordia age
#' @importFrom grDevices rgb
#' @importFrom graphics polygon title points text
#' @importFrom stats pchisq
#' @examples
#' data(examples) 
#' concordia(examples$UPb)
#' @references Ludwig, K.R., 1998. On the treatment of concordant
#'     uranium-lead ages. Geochimica et Cosmochimica Acta, 62(4),
#'     pp.665-676.
#' @export
concordia <- function(x,limits=NULL,alpha=0.05,wetherill=TRUE,show.numbers=FALSE,
                      ellipse.col=rgb(0,1,0,0.5),concordia.col='darksalmon',
                      exterr=TRUE,show.age=1,sigdig=2){
    concordia.line(x,limits=limits,wetherill=wetherill,
                   col=concordia.col,alpha=alpha,exterr=exterr)
    if (show.age==3){
        fit <- discordia.age(x,wetherill=wetherill,exterr=exterr)
        discordia.plot(fit,wetherill=wetherill)
        title(discordia.title(fit,wetherill=wetherill,sigdig=sigdig))
    }
    ns <- length(x)
    for (i in 1:ns){
        if (wetherill) xyc <- wetherill(x,i)
        else xyc <- tera.wasserburg(x,i)
        x0 <- xyc$x[1]
        y0 <- xyc$x[2]
        covmat <- xyc$cov
        ell <- ellipse(x0,y0,covmat,alpha=alpha)
        polygon(ell,col=ellipse.col)
        points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) { text(x0,y0,i) }
    }
    if (show.age==2){
        fit <- concordia.age(x,wetherill=wetherill,exterr=exterr)
        ell <- ellipse(fit$x[1],fit$x[2],fit$cov)
        polygon(ell,col='white')
        title(concordia.title(fit,sigdig=sigdig))
    }
}

# helper function for plot.concordia
concordia.line <- function(x,limits,wetherill,col,alpha=0.05,exterr=TRUE){
    lims <- get.concordia.limits(x,limits=limits,wetherill=wetherill)
    if (wetherill){
        x.lab <- expression(paste(""^"207","Pb/"^"235","U"))
        y.lab <- expression(paste(""^"206","Pb/"^"238","U"))
    } else {
        x.lab <- expression(paste(""^"238","U/"^"206","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"206","Pb"))
    }
    graphics::plot(lims$x,lims$y,type='n',xlab=x.lab,ylab=y.lab)
    range.t <- range(lims$t)
    m <- max(0.8*lims$t[1],lims$t[1]-range.t/20)
    M <- min(1.2*lims$t[2],lims$t[2]+range.t/20)
    nn <- 30 # number of segments into which the concordia line is divided
    if (wetherill) tt <- seq(from=m,to=M,length.out=nn)
    else tt <- exp(seq(from=log(m),to=log(M),length.out=nn))
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
    ticks <- pretty(tt)
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
age_to_concordia_ratios <- function(tt,wetherill=TRUE,exterr=FALSE){
    if (wetherill) return(age_to_wetherill_ratios(tt,exterr=exterr))
    else return(age_to_terawasserburg_ratios(tt,exterr=exterr))
}
get.concordia.limits <- function(x,limits=NULL,wetherill=FALSE){
    out <- list()
    out$x <- c(0,0)
    out$y <- c(0,0)
    if (is.null(limits)) out$t <- c(0,0)
    else out$t <- limits
    nse <- 3 # number of standard errors used for buffer
    if (!is.null(limits) && wetherill){
        out$x <- age_to_Pb207U235_ratio(limits)[,'75']
        out$y <- age_to_Pb206U238_ratio(limits)[,'68']
    } else if (!is.null(limits) && !wetherill){
        if (limits[1] <= 0){
            U238Pb206 <- get.U238Pb206.ratios(x)
            maxx <- max(U238Pb206[,1]+nse*U238Pb206[,2],na.rm=TRUE)
            out$t[1] <- get.Pb206U238.age(1/maxx)[1]
        }
        out$x <- age_to_U238Pb206_ratio(out$t)[,'86']
        out$y <- age_to_Pb207Pb206_ratio(out$t)[,'76']
    } else if (is.null(limits) && wetherill) {
        Pb207U235 <- get.Pb207U235.ratios(x)
        out$x[1] <- min(Pb207U235[,1]-nse*Pb207U235[,2],na.rm=TRUE)
        out$x[2] <- max(Pb207U235[,1]+nse*Pb207U235[,2],na.rm=TRUE)
        Pb206U238 <- get.Pb206U238.ratios(x)
        out$y[1] <- min(Pb206U238[,1]-nse*Pb206U238[,2],na.rm=TRUE)
        out$y[2] <- max(Pb206U238[,1]+nse*Pb206U238[,2],na.rm=TRUE)
        out$t[1] <- get.Pb206U238.age(out$y[1])[1]
        out$t[2] <- get.Pb207U235.age(out$x[2])[1]
    } else if (is.null(limits) && !wetherill){
        U238Pb206 <- get.U238Pb206.ratios(x)
        Pb207Pb206 <- get.Pb207Pb206.ratios(x)
        out$x[1] <- min(U238Pb206[,1]-nse*U238Pb206[,2],na.rm=TRUE)
        out$x[2] <- max(U238Pb206[,1]+nse*U238Pb206[,2],na.rm=TRUE)
        out$y[1] <- min(Pb207Pb206[,1]-nse*Pb207Pb206[,2],na.rm=TRUE)
        out$y[2] <- max(Pb207Pb206[,1]+nse*Pb207Pb206[,2],na.rm=TRUE)
        out$t[1] <- min(get.Pb206U238.age(1/out$x[2])[1],
                        get.Pb207Pb206.age(out$y[1])[1],na.rm=TRUE)
        out$t[2] <- max(get.Pb206U238.age(1/out$x[1])[1],
                        get.Pb207Pb206.age(out$y[2])[1],na.rm=TRUE)
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
concordia.age.default <- function(x,exterr=TRUE,...){
    out <- x
    t.init <- initial.concordia.age(x)
    fit.age <- stats::optim(par=t.init, fn=LL.concordia.age,
                            x=x, exterr=exterr,
                            method="BFGS", hessian=TRUE)
    tt <- fit.age$par
    tt.err <- as.numeric(sqrt(solve(fit.age$hessian)))
    out$age <- c(tt,tt.err)
    names(out$age) <- c('t.conc','s[t.conc]')
    out
}
concordia.age.UPb <- function(x,wetherill=TRUE,exterr=TRUE,...){
    cc <- concordia.comp(x,wetherill=wetherill,exterr=exterr)
    out <- concordia.age(cc)
    mswd <- mswd.concordia(x,cc,out$age[1])
    out$mswd <- mswd$mswd
    out$p.value <- mswd$p.value
    out
}

# x has class 'UPb'
concordia.comp <- function(x,wetherill=TRUE,exterr=TRUE){
    if (wetherill) out <- wetherill(x,1,exterr=exterr)
    else out <- tera.wasserburg(x,1,exterr=exterr)
    fit.comp <- stats::optim(out$x, LL.concordia.comp.default,
                             x=x, wetherill=wetherill,
                             method="BFGS", hessian=TRUE)
    out$x <- fit.comp$par
    out$cov <- solve(fit.comp$hessian)
    out
}

initial.concordia.age <- function(x,...){ UseMethod("initial.concordia.age",x) }
initial.concordia.age.default <- function(x,...){ stop('No default method.') }
initial.concordia.age.wetherill <- function(x,...){
    ages <- c(get.Pb207U235.age(x),get.Pb206U238.age(x))
    mean(ages)
}
initial.concordia.age.terawasserburg <- function(x,...){
    ages <- c(get.Pb206U238.age(x),get.Pb207Pb206.age(x))
    mean(ages)
}

mswd.concordia <- function(x,cc,tt){
    out <- list()
    SS.equivalence <- LL.concordia.comp(cc,x,mswd=TRUE)
    SS.concordance <- LL.concordia.age(tt,cc,mswd=TRUE)
    df.equivalence <- 2*length(x)-1
    df.concordance <- 1
    out$mswd <- list(equivalence = SS.equivalence/df.equivalence,
                     concordance = SS.concordance/df.concordance)
    out$p.value <- list(equivalence = 1-pchisq(SS.equivalence,df.equivalence),
                        concordance = 1-pchisq(SS.concordance,df.concordance))
    out
}

# cc = two element vector
LL.concordia.comp <- function(mu,...){ UseMethod("LL.concordia.comp",mu) }
LL.concordia.comp.default <- function(mu,x,wetherill=TRUE,mswd=FALSE,...){
    out <- 0
    for (i in 1:length(x)){
        if (wetherill) xi <- wetherill(x,i)
        else xi <- tera.wasserburg(x,i)
        X <- matrix(xi$x-mu,1,2)
        covmat <- xi$cov
        if (mswd) out <- out + get.concordia.SS(X,covmat)
        else out <- out + LL.norm(X,covmat)
    }
    out
}
LL.concordia.comp.wetherill <- function(mu,x,mswd=FALSE,...){
    LL.concordia.comp(mu$x,x,wetherill=TRUE,mswd=mswd)
}
LL.concordia.comp.terawasserburg <- function(mu,x,mswd=FALSE,...){
    LL.concordia.comp(mu$x,x,wetherill=FALSE,mswd=mswd)
}

LL.concordia.age <- function(tt,x,mswd=FALSE,exterr=TRUE){
    if (hasClass(x,'wetherill'))
        y <- age_to_wetherill_ratios(tt)
    else if (hasClass(x,'terawasserburg'))
        y <- age_to_terawasserburg_ratios(tt)
    dx <- matrix(x$x-y$x,1,2)
    covmat <- x$cov
    if (exterr) covmat <- x$cov + y$cov
    if (mswd) out <- get.concordia.SS(dx,covmat)
    else out <- LL.norm(dx,covmat)
    out
}

get.concordia.SS <- function(x,covmat){
    x %*% solve(covmat) %*% t(x)
}
