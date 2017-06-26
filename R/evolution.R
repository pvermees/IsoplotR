evolution <- function(x,xlim=NA,ylim=NA,alpha=0.05,
                      show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),
                      lwd=1,line.col='darksalmon',exterr=TRUE,
                      show.age=1,sigdig=2,...){
    fit <- isochron.ThU(x,type=3,plot=FALSE,exterr=exterr)
    d <- data2tit.ThU(x,osmond=TRUE)
    id <- c('X','sX','Z','sZ','rXZ')
    x.lab <- expression(paste(""^"232","Th/"^"238","U"))
    y.lab <- expression(paste(""^"230","Th/"^"238","U"))
    scatterplot(d[,id],xlim=xlim,ylim=ylim,alpha=alpha,
                show.numbers=show.numbers,ellipse.col=ellipse.col,
                a=fit$a[1],b=fit$b[1],line.col=line.col,lwd=lwd)
    title(isochron.title(fit,sigdig=sigdig),xlab=x.lab,ylab=y.lab)
}
