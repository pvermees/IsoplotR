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
#' @importFrom grDevices rgb
#' @examples
#' data(UPb)
#' concordia.plot(UPb)
#' @export
concordia.plot <- function(x,limits=NULL,alpha=0.05,wetherill=TRUE,show.numbers=FALSE,
                           ellipse.col=rgb(0,1,0,0.5),concordia.col='darksalmon',
                           dcu=TRUE){
    concordia.line(x,limits,wetherill,concordia.col,alpha,dcu)
    if (wetherill){
        vars <- c('Pb207U235','Pb206U238')
    } else {
        vars <- c('Pb206U238','Pb207Pb206')
    }
    for (i in 1:nrow(x$x)){
        x0 <- x$x[i,vars[1]]
        y0 <- x$x[i,vars[2]]
        covmat <- get.covmat.UPb(x,i)[vars,vars]
        if (!wetherill){
            J <- matrix(c(-1/(x0*x0),0,0,1),nrow=2)
            covmat <- J %*% covmat %*% t(J)
            x0 <- 1/x0
        }
        ell <- get.ellipse(x0,y0,covmat,alpha=alpha)
        graphics::polygon(ell,col=ellipse.col)
        graphics::points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) { graphics::text(x0,y0,i) }
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
            xc <- UPbratios$x[3]
            yc <- UPbratios$x[2]
        } else {
            xc <- UPbratios$x[4]
            yc <- UPbratios$x[1]
        }
        if (dcu){ # show decay constant uncertainty   
            if (wetherill){ covmat <- UPbratios$cov[c(3,2),c(3,2)] }
            else { covmat <- UPbratios$cov[c(4,1),c(4,1)] }
            if (i > 1) oldell <- ell
            ell <- get.ellipse(xc,yc,covmat,alpha=alpha)
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
            xt <- UPbratios$x[3]
            yt <- UPbratios$x[2]
        } else {
            xt <- UPbratios$x[4]
            yt <- UPbratios$x[1]
        }
        if (dcu){ # show ticks as ellipse
            if (wetherill){ covmat <- UPbratios$cov[c(3,2),c(3,2)] }
            else { covmat <- UPbratios$cov[c(4,1),c(4,1)] }
            ell <- get.ellipse(xt,yt,covmat,alpha=alpha)
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
    if (!is.null(limits) && wetherill){
        out$min.t <- limits[1]
        out$max.t <- limits[2]
        out$min.x <- get.ratios.UPb(out$min.t)$x[3]
        out$max.x <- get.ratios.UPb(out$max.t)$x[3]
        out$min.y <- get.ratios.UPb(out$min.t)$x[2]
        out$max.y <- get.ratios.UPb(out$max.t)$x[2]
    } else if (!is.null(limits) && !wetherill){
        if (limits[1] > 0){
            out$min.t <- limits[1]
        } else {
            out$max.x <- 1/min(X$x[,'Pb206U238']-2*X$x[,'errPb206U238'])
            out$min.t <- get.Pb206U238age(1/out$max.x)
        }
        out$max.t <- limits[2]
        out$min.x <- get.ratios.UPb(out$max.t)$x[4]
        out$min.y <- get.ratios.UPb(out$max.t)$x[1]
        out$max.x <- get.ratios.UPb(out$min.t)$x[4]
        out$max.y <- get.ratios.UPb(out$min.t)$x[1]
    } else if (is.null(limits) && wetherill) {
        out$min.x <- min(X$x[,'Pb207U235']-2*X$x[,'errPb207U235'])
        out$max.x <- max(X$x[,'Pb207U235']+2*X$x[,'errPb207U235'])
        out$min.y <- min(X$x[,'Pb206U238']-2*X$x[,'errPb206U238'])
        out$max.y <- max(X$x[,'Pb206U238']+2*X$x[,'errPb206U238'])
        out$min.t <- get.Pb206U238age(out$min.y)
        out$max.t <- get.Pb207U235age(out$max.x)
    } else if (is.null(limits) && !wetherill){
        out$min.x <- 1/(max(X$x[,'Pb206U238']+2*X$x[,'errPb206U238']))
        out$max.x <- 1/min(X$x[,'Pb206U238']-2*X$x[,'errPb206U238'])
        out$min.y <- min(X$x[,'Pb207Pb206']-2*X$x[,'errPb207Pb206'])
        out$max.y <- max(X$x[,'Pb207Pb206']+2*X$x[,'errPb207Pb206'])
        out$min.t <- get.Pb206U238age(1/out$max.x)
        out$max.t <- get.Pb207Pb206age(out$max.y)
    }
    out
}
