PbPb.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,common.Pb=0){
    if (common.Pb == 0){
        y <- data2york(x,inverse=TRUE)
        PbPb <- y[,c('Y','sY')]
    } else if (common.Pb == 1){
        y <- data2york(x,inverse=FALSE)
        r64 <- y[,'X'] - settings('iratio','Pb206Pb204')[1]
        r74 <- y[,'Y'] - settings('iratio','Pb207Pb204')[1]
        PbPb <- quotient(r74,y[,'sY'],r64,y[,'sX'],y[,'rXY'])
        # alternative implementation: project along inverse isochron:
        # y <- data2york(x,inverse=TRUE)
        # i74 <- settings('iratio','Pb207Pb204')[1]
        # i64 <- settings('iratio','Pb206Pb204')[1]
        # PbPb <- get.76(y,a=i74/i64,b=i74)
    } else if (common.Pb == 2){
        y <- data2york(x,inverse=TRUE)
        fit <- regression(y,model=1)
        PbPb <- get.76(y,a=fit$a[1],b=fit$b[1])
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

get.76 <- function(y,a,b){
    ns <- nrow(y)
    out <- matrix(0,ns,2)
    out[,1] <- 2*a + b*y[,'X'] - y[,'Y']
    J1 <- b
    J2 <- rep(-1,ns)
    E11 <- y[,'sX']^2
    E22 <- y[,'sY']^2
    E12 <- y[,'sX']*y[,'sY']*y[,'rXY']
    out[,2] <- errorprop1x2(J1,J2,E11,E22,E12)
    out
}
