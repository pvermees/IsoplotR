roundit <- function(age,err){
    out <- list()
    out$err <- signif(err,2)
    nd <- log10(trunc(abs(age)/err))+2
    out$x <- signif(age,nd)
    out
}

# set minimum and maximum values of a dataset
getmM <- function(x,from=NA,to=NA,log=FALSE){
    if (is.na(from)) { from <- min(x); getm <- TRUE }
    else { getm <- FALSE }
    if (is.na(to)) { to <- max(x); getM <- TRUE }
    else { getM <- FALSE }
    if (getm) {
        if (log) { from <- from/2 }
        else {
            if (2*from-to<0) {from <- 0}
            else {from <- from-(to-from)/10}
        }
    }
    if (getM) {
        if (log) { to <- 2*to }
        else { to <- to+(to-from)/10 }
    }
    list(m=from,M=to)
}

emptyplot <- function(){
    graphics::plot(c(0,1),c(0,1),type='n',axes=FALSE,xlab="",ylab="")
}
