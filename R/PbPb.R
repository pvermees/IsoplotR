PbPb.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,common.Pb=0,omit4c=NULL){
    if (common.Pb == 0){
        y <- data2york(x,inverse=TRUE)
        PbPb <- y[,c('Y','sY')]
    } else {
        y <- data2york(x,inverse=TRUE)
        if (common.Pb == 1){
            b <- settings('iratio','Pb207Pb204')[1]
        } else if (common.Pb == 2) {
            b <- regression(y,model=1,omit=omit4c)$b[1]
        } else {
            stop('illegal common.Pb option')
        }
        ns <- nrow(y)
        PbPb <- matrix(0,ns,2)
        PbPb[,1] <- y[,'Y'] - b*y[,'X']
        J1 <- -b
        J2 <- rep(1,ns)
        E11 <- y[,'sX']^2
        E22 <- y[,'sY']^2
        E12 <- y[,'sX']*y[,'sY']*y[,'rXY']
        PbPb[,2] <- errorprop1x2(J1,J2,E11,E22,E12)
    }
    PbPb2t(PbPb,exterr=exterr,sigdig=sigdig,i=i)
}

PbPb2t <- function(PbPb,exterr=FALSE,sigdig=NA,i=NA){
    ns <- nrow(PbPb)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    for (j in 1:ns){
        tt <- get.Pb207Pb206.age(PbPb[j,1],PbPb[j,2],exterr=exterr)
        out[j,] <- roundit(tt[1],tt[2],sigdig=sigdig)
    }
    if (!is.na(i)) out <- out[i,]
    out
}
