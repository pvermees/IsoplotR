PbPb_age <- function(x,exterr=FALSE,i=NULL,common.Pb=0,
                     omit4c=NULL,projerr=FALSE){
    y <- data2york(x,inverse=TRUE)
    if (common.Pb == 0){
        PbPb <- y[,c('Y','sY')]
    } else {
        if (common.Pb == 1){
            b <- settings('iratio','Pb207Pb204')
        } else if (common.Pb == 2) {
            b <- regression(y,model=1,omit=omit4c)$b
        } else {
            stop('illegal common.Pb option')
        }
        ns <- nrow(y)
        PbPb <- matrix(0,ns,2)
        PbPb[,1] <- y[,'Y'] - b[1]*y[,'X']
        J1 <- -b[1]
        J2 <- rep(1,ns)
        J3 <- -y[,'X']
        E11 <- y[,'sX']^2
        E22 <- y[,'sY']^2
        E12 <- y[,'sX']*y[,'sY']*y[,'rXY']
        E33 <- b[2]^2
        if (projerr) PbPb[,2] <- sqrt(errorprop1x3(J1,J2,J3,E11,E22,E33,E12))
        else PbPb[,2] <- sqrt(errorprop1x2(J1,J2,E11,E22,E12))
    }
    PbPb2t(PbPb,exterr=exterr,i=i)
}

PbPb2t <- function(PbPb,exterr=FALSE,i=NULL){
    ns <- nrow(PbPb)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    for (j in 1:ns){
        out[j,] <- get_Pb207Pb206_age(PbPb[j,1],PbPb[j,2],exterr=exterr)
    }
    if (!is.null(i)) out <- out[i,]
    out
}
