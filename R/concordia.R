#' Concordia diagram
#'
#' Plot U-Pb data on Wetherill and Tera-Wasserburg concordia diagrams,
#' calculate concordia ages and compositions, evaluates the
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
#' @param sigdig number of significant digits for the concordia/discordia age
#' @importFrom grDevices rgb
#' @importFrom graphics polygon title points text
#' @importFrom stats pchisq
#' @examples
#' data(examples)
#' concordia(examples$UPb)
#' @export
concordia <- function(x,limits=NULL,alpha=0.05,wetherill=TRUE,show.numbers=FALSE,
                      ellipse.col=rgb(0,1,0,0.5),concordia.col='darksalmon',
                      exterr=TRUE,show.age=1,sigdig=2){
    concordia.line(x,limits,wetherill,concordia.col,alpha,exterr)
    if (show.age==3){
        fit <- discordia.age(x,wetherill)
        discordia.plot(fit,wetherill)
        title(discordia.title(fit,wetherill,sigdig=sigdig))
    }
    vars <- get.UPb.labels(wetherill)
    for (i in 1:nrow(x$x)){
        x0 <- x$x[i,vars[1]]
        y0 <- x$x[i,vars[2]]
        covmat <- get.covmat(x,i)[vars,vars]
        ell <- ellipse(x0,y0,covmat,alpha=alpha)
        polygon(ell,col=ellipse.col)
        points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) { text(x0,y0,i) }
    }
    if (show.age==2){
        fit <- concordia.age(x,wetherill,exterr)
        ell <- ellipse(fit$x[1],fit$x[2],fit$cov)
        polygon(ell,col='white')
        title(concordia.title(fit,sigdig=sigdig))
    }
}

# helper function for plot.concordia
concordia.line <- function(X,limits,wetherill,col,alpha=0.05,exterr=TRUE){
    lims <- get.concordia.limits(X,limits,wetherill)
    if (wetherill){
        x.lab <- expression(paste(""^"207","Pb/"^"235","U"))
        y.lab <- expression(paste(""^"206","Pb/"^"238","U"))
    } else {
        x.lab <- expression(paste(""^"238","U/"^"206","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"206","Pb"))
    }
    graphics::plot(c(lims$min.x,lims$max.x),c(lims$min.y,lims$max.y),
                   type='n',xlab=x.lab,ylab=y.lab)
    range.t <- lims$max.t-lims$min.t
    m <- max(0.8*lims$min.t,lims$min.t-range.t/20)
    M <- min(1.2*lims$max.t,lims$max.t+range.t/20)
    nn <- 30 # number of segments into which the concordia line is divided
    if (wetherill){ tt <- seq(from=m,to=M,length.out=nn) }
    else { tt <- exp(seq(from=log(m),to=log(M),length.out=nn)) }
    concordia <- list(x=NULL,y=NULL)
    for (i in 1:nn){ # build the concordia line
        UPbratios <- get.ratios.UPb(tt[i])
        if (wetherill){
            xc <- UPbratios$x['Pb207U235']
            yc <- UPbratios$x['Pb206U238']
        } else {
            xc <- UPbratios$x['U238Pb206']
            yc <- UPbratios$x['Pb207Pb206']
        }
        if (exterr){ # show decay constant uncertainty   
            if (wetherill){ covmat <- UPbratios$cov[c('Pb207U235','Pb206U238'),
                                                    c('Pb207U235','Pb206U238')] }
            else { covmat <- UPbratios$cov[c('U238Pb206','Pb207Pb206'),
                                           c('U238Pb206','Pb207Pb206')] }
            if (i > 1) oldell <- ell
            ell <- ellipse(xc,yc,covmat,alpha=alpha)
            if (i > 1){
                xycd <- rbind(oldell,ell)
                ii <- grDevices::chull(xycd)
                graphics::polygon(xycd[ii,],col=col,border=NA)
            }
        }
        concordia$x <- c(concordia$x,xc)
        concordia$y <- c(concordia$y,yc)
    }
    graphics::lines(concordia$x,concordia$y,col=col,lwd=2)
    # prepare and plot ticks
    ticks <- pretty(tt)
    for (i in 1:length(ticks)){
        UPbratios <- get.ratios.UPb(ticks[i])
        if (wetherill){
            xt <- UPbratios$x['Pb207U235']
            yt <- UPbratios$x['Pb206U238']
        } else {
            xt <- UPbratios$x['U238Pb206']
            yt <- UPbratios$x['Pb207Pb206']
        }
        if (exterr){ # show ticks as ellipse
            if (wetherill){ covmat <- UPbratios$cov[c('Pb207U235','Pb206U238'),
                                                    c('Pb207U235','Pb206U238')] }
            else { covmat <- UPbratios$cov[c('U238Pb206','Pb207Pb206'),
                                           c('U238Pb206','Pb207Pb206')] }
            ell <- ellipse(xt,yt,covmat,alpha=alpha)
            graphics::polygon(ell,col='white')
        } else {
            graphics::points(xt,yt,pch=21,bg='white')
        }
        pos <- 2
        if (exterr & (wetherill & diff(range(concordia$x))<0.05)
            | (!wetherill & diff(range(concordia$x))<2.5)){ pos <- NULL }
        graphics::text(xt,yt,as.character(ticks[i]),pos=pos)
    }
}

get.concordia.limits <- function(X,limits,wetherill){
    out <- list()
    nse <- 3 # number of standard errors used for buffer
    if (!is.null(limits) && wetherill){
        out$min.t <- limits[1]
        out$max.t <- limits[2]
        out$min.x <- get.ratios.UPb(out$min.t)$x['Pb207U235']
        out$max.x <- get.ratios.UPb(out$max.t)$x['Pb207U235']
        out$min.y <- get.ratios.UPb(out$min.t)$x['Pb206U238']
        out$max.y <- get.ratios.UPb(out$max.t)$x['Pb206U238']
    } else if (!is.null(limits) && !wetherill){
        if (limits[1] > 0){
            out$min.t <- limits[1]
        } else {
            out$max.x <- max(X$x[,'U238Pb206']+nse*X$x[,'errU238Pb206'])
            out$min.t <- get.Pb206U238age(1/out$max.x)[1]
        }
        out$max.t <- limits[2]
        out$min.x <- get.ratios.UPb(out$max.t)$x['U238Pb206']
        out$min.y <- get.ratios.UPb(out$max.t)$x['Pb207Pb206']
        out$max.x <- get.ratios.UPb(out$min.t)$x['U238Pb206']
        out$max.y <- get.ratios.UPb(out$min.t)$x['Pb207Pb206']
    } else if (is.null(limits) && wetherill) {
        out$min.x <- min(X$x[,'Pb207U235']-nse*X$x[,'errPb207U235'])
        out$max.x <- max(X$x[,'Pb207U235']+nse*X$x[,'errPb207U235'])
        out$min.y <- min(X$x[,'Pb206U238']-nse*X$x[,'errPb206U238'])
        out$max.y <- max(X$x[,'Pb206U238']+nse*X$x[,'errPb206U238'])
        out$min.t <- get.Pb206U238age(out$min.y)[1]
        out$max.t <- get.Pb207U235age(out$max.x)[1]
    } else if (is.null(limits) && !wetherill){
        out$min.x <- min(X$x[,'U238Pb206']-nse*X$x[,'errU238Pb206'])
        out$max.x <- max(X$x[,'U238Pb206']+nse*X$x[,'errU238Pb206'])
        out$min.y <- min(X$x[,'Pb207Pb206']-nse*X$x[,'errPb207Pb206'])
        out$max.y <- max(X$x[,'Pb207Pb206']+nse*X$x[,'errPb207Pb206'])
        out$min.t <- min(get.Pb206U238age(1/out$max.x)[1],get.Pb207Pb206age(out$min.y)[1])
        out$max.t <- max(get.Pb206U238age(1/out$min.x)[1],get.Pb207Pb206age(out$max.y)[1])
    }
    out
}

concordia.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age,fit$age.err,sigdig=sigdig)
    line1 <- substitute('concordia age ='~a%+-%b~'[Ma] (1'~sigma~')',
                        list(a=rounded.age$x, b=rounded.age$err))
    line2 <- substitute('MSWD (concordance) ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd$concordance,2),
                             b=signif(fit$p.value$concordance,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}

concordia.age <- function(x,...){ UseMethod("concordia.age",x) }
concordia.age.default <- function(x,wetherill=TRUE,exterr=TRUE,...){
    out <- x
    t.init <- initial.concordia.age.guess(out$x,wetherill)
    fit.age <- stats::optim(t.init, LL.concordia.age, x=out$x,
                            covmat=out$cov, wetherill=wetherill,
                            exterr=exterr, method="BFGS", hessian=TRUE)
    out$age <- fit.age$par
    out$age.err <- as.numeric(sqrt(solve(fit.age$hessian)))
    out
}
concordia.age.UPb <- function(x,wetherill=TRUE,exterr=TRUE,i=NA,...){
    if (is.na(i)){
        X <- UPb.preprocess(x,wetherill)
        CC <- concordia.comp(X,wetherill,exterr)
        out <- concordia.age.default(CC,wetherill,exterr,...)
        mswd <- mswd.concordia(X,out$x,out$cov,out$age,wetherill)
        out$mswd <- mswd$mswd
        out$p.value <- mswd$p.value
    } else {
        CC <- UPb.preprocess(x,wetherill,i)
        out <- concordia.age.default(CC,wetherill,exterr,...)
    }
    out
}

# x is a list of lists of U-Pb compositions and covariance
# matrices prepared by the UPb.preprocess function
concordia.comp <- function(x,wetherill=TRUE,exterr=TRUE){
    xy <- initialise.concordant.composition(x,wetherill)
    fit.comp <- stats::optim(xy, LL.concordia.comp, x=x, method="BFGS", hessian=TRUE)
    out <- list()
    out$x <- fit.comp$par
    out$cov <- solve(fit.comp$hessian)
    selection <- names(x[[1]]$x)
    names(out$x) <- selection
    colnames(out$cov) <- selection
    rownames(out$cov) <- selection
    out
}

# x = object of class UPb
# generates a list of lists containing U-Pb/Pb-Pb pairs and
# covariance matrices. Is used to calculate discordia lines
UPb.preprocess <- function(x,wetherill,i=NA){
    selection <- get.selection(x,wetherill)
    if (!is.na(i)){
        X <- x$x[i,selection]
        covmat <- get.covmat(x,i)[selection,selection]
        out <- list(x=X,cov=covmat)
    } else {
        out <- list()
        for (i in 1:nrow(x$x))
            out[[i]] <- UPb.preprocess(x,wetherill,i)
    }
    out
}

# x is a list containing single aliquots of U-Pb analyses
initialise.concordant.composition <- function(x,wetherill=TRUE){
    selection <- get.selection(x,wetherill)
    rs <- x[[1]]$x # running sum
    ns <- length(x) # number of samples
    for (i in 2:length(x)){
        rs <- rs + x[[i]]$x
    }
    rs/length(x)
}

initial.concordia.age.guess <- function(x,wetherill=TRUE){
    if (wetherill){
        Pb207U235age <- get.Pb207U235age(x['Pb207U235'])[1]
        Pb206U238age <- get.Pb206U238age(x['Pb206U238'])[1]
        out <- mean(c(Pb207U235age,Pb206U238age))
    } else {
        Pb206U238age <- get.Pb206U238age(1/x['U238Pb206'])[1]
        Pb207Pb206age <- get.Pb207Pb206age(x['Pb207Pb206'])[1]
        out <- mean(c(Pb206U238age,Pb207Pb206age))
    }
    out
}

mswd.concordia <- function(x,mu,covmat,age,wetherill=TRUE,exterr=TRUE){
    out <- list()
    SS.equivalence <- LL.concordia.comp(mu,x,TRUE)
    SS.concordance <- LL.concordia.age(age,mu,covmat,wetherill,mswd=TRUE,exterr=exterr)
    df.equivalence <- 2*length(x)-1
    df.concordance <- 1
    out$mswd <- list(equivalence = SS.equivalence/df.equivalence,
                     concordance = SS.concordance/df.concordance)
    out$p.value <- list(equivalence = 1-pchisq(SS.equivalence,df.equivalence),
                        concordance = 1-pchisq(SS.concordance,df.concordance))
    out
}

LL.concordia.comp <- function(mu,x,mswd=FALSE){
    out <- 0
    for (i in 1:length(x)){
        X <- matrix(x[[i]]$x-mu,1,2)
        covmat <- x[[i]]$cov
        if (mswd) out <- out + get.concordia.SS(X,covmat)
        else out <- out + LL.norm(X,covmat)
    }
    out
}

LL.concordia.age <- function(age,x,covmat,wetherill=TRUE,mswd=FALSE,exterr=TRUE){
    UPbratios <- get.ratios.UPb(age)
    selection <- get.UPb.labels(wetherill)
    X <- matrix(x[selection]-UPbratios$x[selection],1,2)
    COVMAT <- covmat[selection,selection]
    if (exterr) COVMAT <- COVMAT + UPbratios$cov[selection,selection]
    if (mswd) out <- get.concordia.SS(X,COVMAT)
    else out <- LL.norm(X,COVMAT)
    out
}

get.UPb.labels <- function(wetherill=TRUE){
    if (wetherill) selection <- c('Pb207U235','Pb206U238')
    else selection <- c('U238Pb206','Pb207Pb206')
    selection
}
