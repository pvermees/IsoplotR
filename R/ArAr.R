get.selection.ArAr <- function(x,inverse=TRUE,...){
    if (inverse) selection <- c('Ar39Ar40','Ar36Ar40')
    else selection <- c('Ar39Ar36','Ar40Ar36')
    selection
}

get.covmat.ArAr <- function(x,i,...){
    out <- matrix(rep(0,16),nrow=4)
    rownames(out) <- c('Ar39Ar40','Ar36Ar40','Ar39Ar36','Ar40Ar36')
    colnames(out) <- rownames(out)
    if (x$format == 1){
        out['Ar39Ar40','Ar39Ar40'] <- x$x[i,'errAr39Ar40']^2
        out['Ar36Ar40','Ar36Ar40'] <- x$x[i,'errAr36Ar40']^2
        out['Ar39Ar36','Ar39Ar36'] <- x$x[i,'errAr39Ar36']^2
        out['Ar40Ar36','Ar40Ar36'] <- x$x[i,'errAr40Ar36']^2

        out['Ar39Ar40','Ar36Ar40'] <- get.cov.xzyz(
            x$x[i,'Ar39Ar40'],x$x[i,'errAr39Ar40'],
            x$x[i,'Ar36Ar40'],x$x[i,'errAr36Ar40'],x$x[i,'errAr39Ar36']
        )
        out['Ar39Ar40','Ar39Ar36'] <- get.cov.zxzy(
            x$x[i,'Ar39Ar40'],x$x[i,'errAr39Ar40'],
            x$x[i,'Ar39Ar36'],x$x[i,'errAr39Ar36'],x$x[i,'errAr40Ar36']
        )
        out['Ar39Ar40','Ar40Ar36'] <- get.cov.xzzy(
            x$x[i,'Ar39Ar40'],x$x[i,'errAr39Ar40'],
            x$x[i,'Ar40Ar36'],x$x[i,'errAr40Ar36'],x$x[i,'errAr39Ar36']
        )
        out['Ar36Ar40','Ar39Ar36'] <- get.cov.xzzy(
            x$x[i,'Ar39Ar36'],x$x[i,'errAr39Ar36'],
            x$x[i,'Ar36Ar40'],x$x[i,'errAr36Ar40'],x$x[i,'errAr39Ar40']
        )
        out['Ar36Ar40','Ar40Ar36'] <- -1
        out['Ar39Ar36','Ar40Ar36'] <- get.cov.xzyz(
            x$x[i,'Ar39Ar36'],x$x[i,'errAr39Ar36'],
            x$x[i,'Ar40Ar36'],x$x[i,'errAr40Ar36'],x$x[i,'errAr39Ar40']
        )
        out['Ar36Ar40','Ar39Ar40'] <- out['Ar39Ar40','Ar36Ar40']
        out['Ar39Ar36','Ar39Ar40'] <- out['Ar39Ar40','Ar39Ar36']
        out['Ar40Ar36','Ar39Ar40'] <- out['Ar39Ar40','Ar40Ar36']
        out['Ar39Ar36','Ar36Ar40'] <- out['Ar36Ar40','Ar39Ar36']
        out['Ar40Ar36','Ar36Ar40'] <- out['Ar36Ar40','Ar40Ar36']
        out['Ar40Ar36','Ar39Ar36'] <- out['Ar39Ar36','Ar40Ar36']
    }
    out
}
