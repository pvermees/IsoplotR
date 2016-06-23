get.covmat.ArAr <- function(x,i){
    covmat <- matrix(rep(0,16),nrow=4)
    rownames(covmat) <- c('Ar39Ar40','Ar36Ar40','Ar39Ar36','Ar40Ar36')
    colnames(covmat) <- rownames(covmat)
    if (x$format == 1){
        covmat['Ar39Ar40','Ar39Ar40'] <- x$x[i,'errAr39Ar40']^2
        covmat['Ar36Ar40','Ar36Ar40'] <- x$x[i,'errAr36Ar40']^2
        covmat['Ar39Ar36','Ar39Ar36'] <- x$x[i,'errAr39Ar36']^2
        covmat['Ar40Ar36','Ar40Ar36'] <- x$x[i,'errAr40Ar36']^2

        covmat['Ar39Ar40','Ar36Ar40'] <- get.covmat.xzyz(
            x$x[i,'Ar39Ar40'],x$x[i,'errAr39Ar40'],
            x$x[i,'Ar36Ar40'],x$x[i,'errAr36Ar40'],x$x[i,'errAr39Ar36']
        )
        covmat['Ar39Ar40','Ar39Ar36'] <- get.covmat.zxzy(
            x$x[i,'Ar39Ar40'],x$x[i,'errAr39Ar40'],
            x$x[i,'Ar39Ar36'],x$x[i,'errAr39Ar36'],x$x[i,'errAr40Ar36']
        )
        covmat['Ar39Ar40','Ar40Ar36'] <- get.covmat.xzzy(
            x$x[i,'Ar39Ar40'],x$x[i,'errAr39Ar40'],
            x$x[i,'Ar40Ar36'],x$x[i,'errAr40Ar36'],x$x[i,'errAr39Ar36']
        )
        covmat['Ar36Ar40','Ar39Ar36'] <- get.covmat.xzzy(
            x$x[i,'Ar39Ar36'],x$x[i,'errAr39Ar36'],
            x$x[i,'Ar36Ar40'],x$x[i,'errAr36Ar40'],x$x[i,'errAr39Ar40']
        )
        covmat['Ar36Ar40','Ar40Ar36'] <- -1
        covmat['Ar39Ar36','Ar40Ar36'] <- get.covmat.xzyz(
            x$x[i,'Ar39Ar36'],x$x[i,'errAr39Ar36'],
            x$x[i,'Ar40Ar36'],x$x[i,'errAr40Ar36'],x$x[i,'errAr39Ar40']
        )
        covmat['Ar36Ar40','Ar39Ar40'] <- covmat['Ar39Ar40','Ar36Ar40']
        covmat['Ar39Ar36','Ar39Ar40'] <- covmat['Ar39Ar40','Ar39Ar36']
        covmat['Ar40Ar36','Ar39Ar40'] <- covmat['Ar39Ar40','Ar40Ar36']
        covmat['Ar39Ar36','Ar36Ar40'] <- covmat['Ar36Ar40','Ar39Ar36']
        covmat['Ar40Ar36','Ar36Ar40'] <- covmat['Ar36Ar40','Ar40Ar36']
        covmat['Ar40Ar36','Ar39Ar36'] <- covmat['Ar39Ar36','Ar40Ar36']
    }
    covmat
}
