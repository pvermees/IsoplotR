#' Concordia diagram
#'
#' Wetherill and Tera-Wasserburg concordia diagrams
#'
#' @param x an object of class \code{UPb}
#' @param limits age limits of the concordia line
#' @param alpha confidence cutoff for the error ellipses
#' @param wetherill boolean flag (FALSE for Tera-Wasserburg)
#' @param show.numbers boolean flag (TRUE to show grain numbers)
#' @param ellipse.col background colour of the error ellipses
#' @param concordia.col colour of the concordia line
#' @param dcu show decay constant uncertainty?
#' @param show.age one of either
#'
#' \code{0}: don't show the age
#'
#' \code{1}: calculate the concordia age
#'
#' \code{2}: fit a discordia line
#'
#' @importFrom grDevices rgb
#' @importFrom graphics polygon title points text
#' @importFrom stats pchisq
#' @examples
#' data(UPb)
#' concordia.plot(UPb)
#' @export
concordia.plot <- function(x,limits=NULL,alpha=0.05,wetherill=TRUE,show.numbers=FALSE,
                           ellipse.col=rgb(0,1,0,0.5),concordia.col='darksalmon',
                           dcu=TRUE, show.age=0){
    concordia.line(x,limits,wetherill,concordia.col,alpha,dcu)
    if (show.age==2){
        fit <- discordia.age(x,wetherill)
        discordia.plot(fit,wetherill)
        title(discordia.title(fit,wetherill))
    }
    vars <- get.UPb.labels(wetherill)
    for (i in 1:nrow(x$x)){
        x0 <- x$x[i,vars[1]]
        y0 <- x$x[i,vars[2]]
        covmat <- get.covmat.UPb(x,i)[vars,vars]
        ell <- ellipse(x0,y0,covmat,alpha=alpha)
        polygon(ell,col=ellipse.col)
        points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) { text(x0,y0,i) }
    }
    if (show.age==1){
        fit <- concordia.age(x,wetherill,dcu)
        ell <- ellipse(fit$x[1],fit$x[2],fit$x.cov)
        polygon(ell,col='white')
        title(concordia.title(fit))
    }
}

# helper function for plot.concordia
concordia.line <- function(X,limits,wetherill,col,alpha=0.05,dcu=TRUE){
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
    nn <- 30
    if (wetherill){ tt <- seq(from=m,to=M,length.out=nn) }
    else { tt <- exp(seq(from=log(m),to=log(M),length.out=nn)) }
    concordia <- list(x=NULL,y=NULL)
    for (i in 1:nn){
        UPbratios <- get.ratios.UPb(tt[i])
        if (wetherill){
            xc <- UPbratios$x['Pb207U235']
            yc <- UPbratios$x['Pb206U238']
        } else {
            xc <- UPbratios$x['U238Pb206']
            yc <- UPbratios$x['Pb207Pb206']
        }
        if (dcu){ # show decay constant uncertainty   
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
        if (dcu){ # show ticks as ellipse
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
        if (dcu & (wetherill & diff(range(concordia$x))<0.05)
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
            out$min.t <- get.Pb206U238age(1/out$max.x)
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
        out$min.t <- get.Pb206U238age(out$min.y)
        out$max.t <- get.Pb207U235age(out$max.x)
    } else if (is.null(limits) && !wetherill){
        out$min.x <- min(X$x[,'U238Pb206']-nse*X$x[,'errU238Pb206'])
        out$max.x <- max(X$x[,'U238Pb206']+nse*X$x[,'errU238Pb206'])
        out$min.y <- min(X$x[,'Pb207Pb206']-nse*X$x[,'errPb207Pb206'])
        out$max.y <- max(X$x[,'Pb207Pb206']+nse*X$x[,'errPb207Pb206'])
        out$min.t <- min(get.Pb206U238age(1/out$max.x),get.Pb207Pb206age(out$min.y))
        out$max.t <- max(get.Pb206U238age(1/out$min.x),get.Pb207Pb206age(out$max.y))
    }
    out
}

discordia.title <- function(fit,wetherill){
    if (wetherill){
        lower.age <- roundit(fit$x[1],sqrt(fit$cov[1,1]))
        upper.age <- roundit(fit$x[2],sqrt(fit$cov[2,2]))
        line1 <- substitute('lower intercept ='~a%+-%b~'[Ma]',
                            list(a=lower.age$x, b=lower.age$err))
        line2 <- substitute('upper intercept ='~a%+-%b~'[Ma]',
                            list(a=upper.age$x, b=upper.age$err))
    } else {
        lower.age <- roundit(fit$x[1],sqrt(fit$cov[1,1]))
        intercept <- roundit(fit$x[2],sqrt(fit$cov[2,2]))
        line1 <- substitute('age ='~a%+-%b~'[Ma]',
                            list(a=lower.age$x, b=lower.age$err))
        line2 <- substitute('('^207*'Pb/'^206*'Pb)'[0]~'='~a%+-%b,
                              list(a=intercept$x, b=intercept$err))
    }
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}

concordia.title <- function(fit){
    rounded.age <- roundit(fit$age,fit$age.err)
    line1 <- substitute('concordia age ='~a%+-%b~'[Ma] (1'~sigma~')',
                        list(a=rounded.age$x, b=rounded.age$err))
    line2 <- substitute('MSWD (concordance) ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd$concordance,2),
                             b=signif(fit$p.value$concordance,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}
