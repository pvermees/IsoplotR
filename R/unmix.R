#' @title Unmix detrital age distributions
#'
#' @description
#' Estimates the proportions of specified source datasets in each sink
#' dataset by fitting a weighted combination of the sources' empirical
#' cumulative distribution functions.
#'
#' @details
#' The input is a named list of numeric vectors. The first
#' \code{nsources} vectors are treated as sources by default; all
#' remaining vectors are treated as sinks. The source proportions are
#' constrained to be between zero and one and to sum to one. With
#' \code{boot=TRUE}, the source observations are resampled to estimate
#' confidence limits for the proportions.
#'
#' @param x a named list of numeric vectors containing the source and
#'     sink datasets.
#' @param nsources the number of source datasets at the start of
#'     \code{x}, used when \code{source_names} is not supplied.
#' @param source_names names of the source datasets. Datasets in
#'     \code{x} not named here are treated as sinks.
#' @param plot logical flag indicating whether to plot the empirical
#'     cumulative distributions, fitted distributions, and estimated
#'     proportions.
#' @param boot logical flag indicating whether to bootstrap the source
#'     datasets to estimate confidence limits for the proportions.
#' @param nboot number of bootstrap replicates.
#' @param recurse logical flag indicating whether or not to
#'     recursively test solutions in which one or more source
#'     proportions are zero.
#' @param hide vector with indices of samples that should be removed
#'     from the plot.
#' @param ... additional graphical parameters passed to the plots.
#' @return If \code{boot=FALSE}, a matrix with one row per sink and
#'     columns for the estimated source proportions and the Cramer-von
#'     Mises fit statistic. If \code{boot=TRUE}, a list containing
#'     this matrix (\code{prop}) and matrices of lower (\code{ll}) and
#'     upper (\code{ul}) bootstrap confidence limits for the source
#'     proportions.
#' @examples
#' attach(examples)
#' fit <- unmix(DZ, sources=c('N1','N4','N14'), plot=FALSE)
#' print(fit)
#' @export
unmix <- function(x,
                  nsources=2,
                  source_names=names(x)[1:nsources],
                  plot=TRUE,
                  boot=(length(x)-length(source_names)==1),
                  nboot=500,
                  recurse=TRUE,
                  hide=NULL,...){
    # 1. data prep
    wrong_names <- !(source_names %in% names(x))
    if (any(wrong_names)){
        missing_names <- paste(source_names[wrong_names],collapse=", ")
        num_missing <- length(missing_names)
        do_es <- ifelse(num_missing>1," do "," does ")
        stop(missing_names,do_es,"not exist in this dataset.")
    }
    if (is.character(hide)){
        hide <- which(names(x)%in%hide)
    }
    x2calc <- clear(x,hide)
    sorted_data <- lapply(x2calc, sort)
    sample_names <- names(sorted_data)
    is_source <- sample_names %in% source_names
    sink_names <- sample_names[!is_source]
    num_sources <- length(source_names)
    if (num_sources<2){
        stop("You must specify at least two sources")
    }
    num_sinks <- length(sink_names)
    prop <- matrix(NA,nrow=num_sinks,ncol=num_sources+1)
    rownames(prop) <- sink_names
    colnames(prop) <- c(source_names,'CvM')
    # 2. get the mixing proportions
    XYWCvM <- list()
    prop <- matrix(NA,nrow=num_sinks,ncol=num_sources+1)
    rownames(prop) <- sink_names
    colnames(prop) <- c(source_names,'CvM')
    for (sink_name in sink_names){
        sink_data <- sorted_data[[sink_name]]
        length_sink <- length(sink_data)
        X <- matrix(0,nrow=length_sink,ncol=num_sources)
        colnames(X) <- source_names
        for (source_name in colnames(X)){
            source_data <- sorted_data[[source_name]]
            X[,source_name] <- stats::ecdf(source_data)(sink_data)
        }
        Y <- matrix(seq(from=1/(2*length_sink),
                        to=(2*length_sink-1)/(2*length_sink),
                        length.out=length_sink),
                    nrow=length_sink,ncol=1)
        XYWCvM[[sink_name]] <- list(X=X,Y=Y)
        fit <- XY2WCvM(X=X,Y=Y,recurse=recurse)
        XYWCvM[[sink_name]]$W <- prop[sink_name,source_names] <- fit$W
        XYWCvM[[sink_name]]$CvM <- prop[sink_name,'CvM'] <- fit$CvM
    }
    if (boot){
        W_replicates <- array(NA, dim=c(nboot, num_sinks, num_sources))
        dimnames(W_replicates) <- list(NULL, sink_names, source_names)
        for (b in 1:nboot) {
            x_boot <- x2calc
            for (source_name in source_names) {
                n_source <- length(x2calc[[source_name]])
                idx <- sample(1:n_source, replace=TRUE)
                x_boot[[source_name]] <- x2calc[[source_name]][idx]
            }
            fit <- unmix(x_boot, nsources=nsources,
                         source_names=source_names,
                         plot=FALSE, boot=FALSE,
                         recurse=recurse)
            W_replicates[b, , ] <- fit[1:num_sources]
        }
        alpha <- settings('alpha')
        ll <- apply(W_replicates, c(2,3), stats::quantile, probs=alpha/2)
        ul <- apply(W_replicates, c(2,3), stats::quantile, probs=1-alpha/2)
        out <- list(prop=prop,ll=ll,ul=ul)
    } else {
        out <- prop
    }
    if (plot){
        plot_unmix(XYWCvM=XYWCvM,
                   sorted_data=sorted_data,
                   source_names=source_names,
                   sink_names=sink_names,...)
    }
    out
}

XY2WCvM <- function(X,Y,recurse=TRUE){
    num_sources <- ncol(X)
    tXXinv <- solve(t(X)%*%X)
    WOLS <- tXXinv%*%t(X)%*%Y
    c1 <- matrix(1,nrow=num_sources,ncol=1)
    r1 <- matrix(1,nrow=1,ncol=num_sources)
    fact <- (1-r1%*%WOLS)/(r1%*%tXXinv%*%c1)
    W <- WOLS + fact[1,1] * tXXinv %*% c1
    CvM <- sum((X %*% W - Y)^2)
    bad <- any(W<0) || any(W>1)
    if (recurse && num_sources>1 && bad){
        CvM <- Inf
        for (i in 1:num_sources){
            fit <- XY2WCvM(X=X[,-i,drop=FALSE],Y=Y,recurse=recurse)
            if (!any(fit$W < 0) && !any(fit$W > 1) && fit$CvM < CvM){
                W <- rep(0,num_sources)
                W[-i] <- fit$W
                CvM <- fit$CvM
            }
        }
    }
    list(W=W,CvM=CvM)
}

empty_ecdf_plot <- function(data_range,...){
    plot(x=data_range,y=c(0,1),type='n',xaxt='n',
         yaxt='n',xlab=NULL,ylab=NULL,...)
}

plot_unmix <- function(XYWCvM,
                       sorted_data,
                       source_names,sink_names,
                       xlab='t (Ma)',...){
    num_sources <- length(source_names)
    num_sinks <- length(sink_names)
    cads <- c(3,3:(2+num_sinks),2+num_sinks)
    zigzag <- c(0,rep(3+num_sinks,num_sinks),0)
    top <- c(0,1,2,2,0)
    bottom <- cbind(0,cads,0,zigzag,0)
    layout_matrix <- rbind(0,top,bottom,0)
    w <- c(1,6,1,12,1)
    if (num_sinks>1){
        h <- c(3,4,2,2,rep(4,num_sinks-2),2,2,3)
    } else {
        h <- c(3,4,2,2,2,3)
    }
    graphics::layout(layout_matrix,widths=w,heights=h)
    colours <- grDevices::hcl.colors(num_sources)
    op <- graphics::par(mar=rep(0,4),mgp=c(1.5,0.75,0))
    graphics::plot.new()
    graphics::legend('topright',legend=source_names,
                     bty='n',lty=rep(1,num_sources),col=colours)
    graphics::legend('topleft',legend=c('observed','fitted'),
                     bty='n',lty=rep(1,2),col=c('red','blue'))
    data_range <- range(unlist(sorted_data))
    empty_ecdf_plot(data_range=data_range,...)
    for (i in 1:num_sources){
        source_name <- source_names[i]
        samp <- sorted_data[[source_name]]
        graphics::plot(stats::ecdf(samp),verticals=TRUE,pch=NA,
                       col=colours[i],add=TRUE,col.01line=NULL)
    }
    graphics::axis(side=3)
    graphics::mtext(text=xlab,side=3,line=1.75,cex=0.8)
    W <- matrix(NA,nrow=num_sinks,ncol=num_sources)
    rownames(W) <- sink_names
    colnames(W) <- source_names
    for (sink_name in sink_names){
        samp <- sorted_data[[sink_name]]
        Y <- XYWCvM[[sink_name]]$X %*% XYWCvM[[sink_name]]$W
        empty_ecdf_plot(data_range=data_range,...)
        graphics::plot(stats::ecdf(samp),verticals=TRUE,
                       pch=NA,main='',col='red',xaxt='n',yaxt='n',
                       xlab=NULL,ylab=NULL,add=TRUE,col.01line=NULL)
        graphics::lines(x=samp,y=c(Y),type='s',col='blue')
        CvM <- signif(XYWCvM[[sink_name]]$CvM,2)
        graphics::legend('bottomright',
                         legend=c(sink_name,paste0('CvM = ',CvM)),
                         bty='n',text.col=c('blue','black'))
        W[sink_name,] <- 100*XYWCvM[[sink_name]]$W
    }
    graphics::axis(side=1)
    graphics::mtext(text=xlab,side=1,line=1.75,cex=0.8)
    if (num_sinks>1){
        graphics::plot(x=c(0,100),y=c(1,num_sinks),type='n',
                       xaxs = "i", yaxs = "i",
                       xaxt='n',yaxt='n',
                       xlab=NULL,ylab=NULL)
        sum_of_weights <- cbind(0,t(apply(W,1,cumsum)))
        for (i in num_sources:1){
            graphics::polygon(x=c(0,sum_of_weights[,i+1],0),
                              y=c(num_sinks,num_sinks:1,1),
                              col=colours[i])
            xtext <- (sum_of_weights[,i+1]+sum_of_weights[,i])/2
            graphics::text(x=xtext,y=num_sinks:1,labels=round(W[,i]),
                           xpd=NA,col='grey40')
        }
        graphics::axis(side=1)
        graphics::mtext(text='%',side=1,line=1.75,cex=0.8)
    } else {
        graphics::plot(1,type="n",xlim=c(0,10),ylim=c(0,10),
                       xlab="",ylab="",axes=FALSE)
        graphics::legend('center',
                         legend=paste0(source_names,' = ',
                                       round(W[1,]),'%'),
                         bty='n',cex=1.2,xpd=NA)
    }
    graphics::par(op)
}
