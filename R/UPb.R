wetherill <- function(x,i=1,format){
    if (missing(format)){
        X <- x$x
        format <- x$format
    } else {
        X <- x
    }
    out <- list()
    if (format < 4) labels <- c('Pb207U235','Pb206U238')
    else if (format < 7) labels <- c('Pb207U235','Pb206U238','Pb204U238')
    else if (format < 9) labels <- c('Pb207U235','Pb206U238','Pb208Th232','Th232U238')
    else if (format == 9) labels <- c('Pb206U238','Pb204U238')
    else if (format == 10) labels <- c('Pb207U235','Pb204U235')
    else if (format == 11) labels <- c('Pb206U238','Pb208Th232','Th232U238')
    else if (format == 12) labels <- c('Pb207U235','Pb208Th232','Th232U238')
    else if (format == 85) labels <- c('Pb207U235','Pb206U238','Pb208U238')
    else if (format == 119) labels <- c('Pb206U238','Pb208U238')
    else if (format == 1210) labels <- c('Pb207U235','Pb208U235')
    else stop('Invalid U-Pb format for wetherill() function.')
    if (format %in% c(1,3)){
        out$x <- X[i,labels]
        out$cov <- cor2cov2(X[i,'errPb207U235'],X[i,'errPb206U238'],X[i,'rXY'])
    } else if (format == 2){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*X[i,'Pb207Pb206']/X[i,'U238Pb206']
        Pb206U238 <- 1/X[i,'U238Pb206']
        out$x <- c(Pb207U235,Pb206U238)
        J <- matrix(0,2,2)
        J[1,1] <- -Pb207U235/X[i,'U238Pb206']
        J[1,2] <- U238U235/X[i,'U238Pb206']
        J[2,1] <- -1/X[i,'U238Pb206']^2
        E <- cor2cov2(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],X[i,'rXY'])
        out$cov <- J %*% E %*% t(J)
    } else if (format==4){
        out$x <- X[i,labels]
        out$cov <- cor2cov3(X[i,'errPb207U235'],X[i,'errPb206U238'],
                            X[i,'errPb204U238'],X[i,'rXY'],
                            X[i,'rXZ'],X[i,'rYZ'])
    } else if (format == 5){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*X[i,'Pb207Pb206']/X[i,'U238Pb206']
        Pb206U238 <- 1/X[i,'U238Pb206']
        Pb204U238 <- X[i,'Pb204Pb206']/X[i,'U238Pb206']
        out$x <- c(Pb207U235,Pb206U238,Pb204U238)
        J <- matrix(0,3,3)
        J[1,1] <- -Pb207U235/X[i,'U238Pb206']
        J[1,2] <- U238U235/X[i,'U238Pb206']
        J[2,1] <- -Pb206U238/X[i,'U238Pb206']
        J[3,1] <- -Pb204U238/X[i,'U238Pb206']
        J[3,3] <- 1/X[i,'U238Pb206']
        E <- cor2cov3(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],
                      X[i,'errPb204Pb206'],X[i,'rXY'],
                      X[i,'rXZ'],X[i,'rYZ'])
        out$cov <- J %*% E %*% t(J)        
    } else if (format == 6){
        out$x <- X[i,labels]
        out$cov <- matrix(0,3,3)
        diag(out$cov) <-
            X[i,c('errPb207U235','errPb206U238','errPb204U238')]^2
        out$cov[1,2] <-
            get_cov_75_68(X[i,'Pb207U235'],X[i,'errPb207U235'],
                          X[i,'Pb206U238'],X[i,'errPb206U238'],
                          X[i,'Pb207Pb206'],X[i,'errPb207Pb206'])
        out$cov[1,3] <-
            get_cov_75_48(X[i,'Pb207U235'],X[i,'errPb207U235'],
                          X[i,'Pb204U238'],X[i,'errPb204U238'],
                          X[i,'Pb204Pb207'],X[i,'errPb204Pb207'])
        out$cov[2,3] <-
            get_cov_68_48(X[i,'Pb206U238'],X[i,'errPb206U238'],
                          X[i,'Pb204U238'],X[i,'errPb204U238'],
                          X[i,'Pb204Pb206'],X[i,'errPb204Pb206'])
        out$cov[2,1] <- out$cov[1,2]
        out$cov[3,1] <- out$cov[1,3]
        out$cov[3,2] <- out$cov[2,3]
    } else if (format == 7){
        out$x <- X[i,labels]
        out$cov <- cor2cov4(X[i,'errPb207U235'],X[i,'errPb206U238'],
                            X[i,'errPb208Th232'],X[i,'errTh232U238'],
                            X[i,'rXY'],X[i,'rXZ'],X[i,'rXW'],
                            X[i,'rYZ'],X[i,'rYW'],X[i,'rZW'])
    } else if (format == 8){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*X[i,'Pb207Pb206']/X[i,'U238Pb206']
        Pb206U238 <- 1/X[i,'U238Pb206']
        Pb208Th232 <- X[i,'Pb208Pb206']/(X[i,'U238Pb206']*X[i,'Th232U238'])
        Th232U238 <- X[i,'Th232U238']
        out$x <- c(Pb207U235,Pb206U238,Pb208Th232,Th232U238)
        J <- matrix(0,4,4)
        J[1,1] <- -Pb207U235/X[i,'U238Pb206']
        J[1,2] <- U238U235/X[i,'U238Pb206']
        J[2,1] <- -Pb206U238/X[i,'U238Pb206']
        J[3,1] <- -Pb208Th232/X[i,'U238Pb206']
        J[3,3] <- 1/(X[i,'U238Pb206']*X[i,'Th232U238'])
        J[3,4] <- -Pb208Th232/X[i,'Th232U238']
        J[4,4] <- 1
        E <- cor2cov4(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],
                      X[i,'errPb208Pb206'],X[i,'errTh232U238'],
                      X[i,'rXY'],X[i,'rXZ'],X[i,'rXW'],
                      X[i,'rYZ'],X[i,'rYW'],X[i,'rZW'])
        out$cov <- J %*% E %*% t(J)
    } else if (format == 9){
        Pb206U238 <- 1/X[i,'U238Pb206']
        Pb204U238 <- X[i,'Pb204Pb206']/X[i,'U238Pb206']
        out$x <- c(Pb206U238,Pb204U238)
        J <- matrix(0,2,2)
        J[1,1] <- -Pb206U238/X[i,'U238Pb206']
        J[2,1] <- -Pb204U238/X[i,'U238Pb206']
        J[2,2] <- 1/X[i,'U238Pb206']
        E <- cor2cov2(X[i,'errU238Pb206'],X[i,'errPb204Pb206'],X[i,'rXY'])
        out$cov <- J %*% E %*% t(J)
    } else if (format == 10){
        Pb207U235 <- 1/X[i,'U235Pb207']
        Pb204U235 <- X[i,'Pb204Pb207']/X[i,'U235Pb207']
        out$x <- c(Pb207U235,Pb204U235)
        J <- matrix(0,2,2)
        J[1,1] <- -Pb207U235/X[i,'U235Pb207']
        J[2,1] <- -Pb204U235/X[i,'U235Pb207']
        J[2,2] <- 1/X[i,'U235Pb207']
        E <- cor2cov2(X[i,'errU235Pb207'],X[i,'errPb204Pb207'],X[i,'rXY'])
        out$cov <- J %*% E %*% t(J)
    } else if (format == 11){
        Pb206U238 <- 1/X[i,'U238Pb206']
        Pb208Th232 <- X[i,'Pb208Pb206']/(X[i,'U238Pb206']*X[i,'Th232U238'])
        Th232U238 <- X[i,'Th232U238']
        out$x <- c(Pb206U238,Pb208Th232,Th232U238)
        J <- matrix(0,3,3)
        J[1,1] <- -Pb206U238/X[i,'U238Pb206']
        J[2,1] <- -Pb208Th232/X[i,'U238Pb206']
        J[2,2] <- 1/(X[i,'U238Pb206']*X[i,'Th232U238'])
        J[2,3] <- -Pb208Th232/X[i,'Th232U238']
        J[3,3] <- 1
        E <- cor2cov3(X[i,'errU238Pb206'],X[i,'errPb208Pb206'],X[i,'errTh232U238'],
                      X[i,'rXY'],X[i,'rXZ'],X[i,'rYZ'])
        out$cov <- J %*% E %*% t(J)
    } else if (format == 12){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- 1/X[i,'U235Pb207']
        Pb208Th232 <- X[i,'Pb208Pb207']/(X[i,'U235Pb207']*U238U235*X[i,'Th232U238'])
        Th232U238 <- X[i,'Th232U238']
        out$x <- c(Pb207U235,Pb208Th232,Th232U238)
        J <- matrix(0,3,3)
        J[1,1] <- -Pb207U235/X[i,'U235Pb207']
        J[2,1] <- -Pb208Th232/X[i,'U235Pb207']
        J[2,2] <- 1/(X[i,'U235Pb207']*U238U235*X[i,'Th232U238'])
        J[2,3] <- -Pb208Th232/X[i,'Th232U238']
        J[3,3] <- 1
        E <- cor2cov3(X[i,'errU235Pb207'],X[i,'errPb208Pb207'],X[i,'errTh232U238'],
                      X[i,'rXY'],X[i,'rXZ'],X[i,'rYZ'])
        out$cov <- J %*% E %*% t(J)
    } else if (format == 85){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*X[i,'Pb207Pb206']/X[i,'U238Pb206']
        Pb206U238 <- 1/X[i,'U238Pb206']
        Pb208U238 <- X[i,'Pb208Pb206']/X[i,'U238Pb206']
        out$x <- c(Pb207U235,Pb206U238,Pb208U238)
        J <- matrix(0,3,3)
        J[1,1] <- -Pb207U235/X[i,'U238Pb206']
        J[1,2] <- U238U235/X[i,'U238Pb206']
        J[2,1] <- -Pb206U238/X[i,'U238Pb206']
        J[3,1] <- -Pb208U238/X[i,'U238Pb206']
        J[3,3] <- 1/X[i,'U238Pb206']
        E <- cor2cov3(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],
                      X[i,'errPb208Pb206'],X[i,'rXY'],X[i,'rXZ'],X[i,'rYZ'])
        out$cov <- J %*% E %*% t(J)
    } else if (format == 119){
        Pb206U238 <- 1/X[i,'U238Pb206']
        Pb208U238 <- X[i,'Pb208Pb206']/X[i,'U238Pb206']
        out$x <- c(Pb206U238,Pb208U238)
        J <- matrix(0,2,2)
        J[1,1] <- -Pb206U238/X[i,'U238Pb206']
        J[2,1] <- -Pb208U238/X[i,'U238Pb206']
        J[2,2] <- 1/X[i,'U238Pb206']
        E <- cor2cov2(X[i,'errU238Pb206'],X[i,'errPb208Pb206'],X[i,'rXY'])
        out$cov <- J %*% E %*% t(J)
    } else if (format == 1210){
        Pb207U235 <- 1/X[i,'U235Pb207']
        Pb208U235 <- X[i,'Pb208Pb207']/X[i,'U235Pb207']
        out$x <- c(Pb207U235,Pb208U235)
        J <- matrix(0,2,2)
        J[1,1] <- -Pb207U235/X[i,'U235Pb207']
        J[2,1] <- -Pb208U235/X[i,'U235Pb207']
        J[2,2] <- 1/X[i,'U235Pb207']
        E <- cor2cov2(X[i,'errU235Pb207'],X[i,'errPb208Pb207'],X[i,'rXY'])
        out$cov <- J %*% E %*% t(J)
    }
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    class(out) <- "wetherill"
    out
}
tera_wasserburg <- function(x,i=1,format){
    if (missing(format)){
        X <- x$x
        format <- x$format
    } else {
        X <- x
    }
    out <- list()
    if (format < 4) labels <- c('U238Pb206','Pb207Pb206')
    else if (format < 7) labels <- c('U238Pb206','Pb207Pb206','Pb204Pb206')
    else if (format < 9) labels <- c('U238Pb206','Pb207Pb206','Pb208Pb206','Th232U238')
    else if (format == 9) labels <- c('U238Pb206','Pb204Pb206')
    else if (format == 10) labels <- c('U235Pb207','Pb204Pb207')
    else if (format == 11) labels <- c('U238Pb206','Pb208Pb206','Th232U238')
    else if (format == 12) labels <- c('U235Pb207','Pb208Pb207','Th232U238')
    else if (format == 85) labels <- c('U238Pb206','Pb207Pb206','Pb208Pb206')
    else if (format == 119) labels <- c('U238Pb206','Pb208Pb206')
    else if (format == 1210) labels <- c('U235Pb207','Pb208Pb207')
    else stop('tera_wasserburg() is not available for this U-Pb format')
    if (format==1){
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/X[i,'Pb206U238']
        Pb207Pb206 <- X[i,'Pb207U235']/(X[i,'Pb206U238']*U238U235)
        J <- matrix(0,2,2)
        J[1,2] <- -1/X[i,'Pb206U238']^2
        J[2,1] <- 1/(X[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/X[i,'Pb206U238']
        E <- cor2cov2(X[i,'errPb207U235'],X[i,'errPb206U238'],X[i,'rXY'])
        out$x <- c(U238Pb206,Pb207Pb206)
        out$cov <- J %*% E %*% t(J)
    } else if (format == 2){
        out$x <- X[i,labels]
        out$cov <- cor2cov2(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],X[i,'rXY'])
    } else if (format == 3){
        U238Pb206 <- 1/X[i,'Pb206U238']
        out$x <- c(U238Pb206,X[i,'Pb207Pb206'])
        J <- matrix(0,2,2)
        E <- cor2cov2(X[i,'errPb206U238'],X[i,'errPb207Pb206'],X[i,'rYZ'])
        J[1,1] <- -U238Pb206^2
        J[2,2] <- 1
        out$cov <- J %*% E %*% t(J)
    } else if (format == 4){
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/X[i,'Pb206U238']
        Pb207Pb206 <- X[i,'Pb207U235']/(X[i,'Pb206U238']*U238U235)
        Pb204Pb206 <- X[i,'Pb204U238']/X[i,'Pb206U238']
        E <- cor2cov3(X[i,'errPb207U235'],X[i,'errPb206U238'],
                      X[i,'errPb204U238'],X[i,'rXY'],
                      X[i,'rXZ'],X[i,'rYZ'])
        J <- matrix(0,3,3)
        J[1,2] <- -U238Pb206/X[i,'Pb206U238']
        J[2,1] <- 1/(X[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/X[i,'Pb206U238']
        J[3,2] <- -Pb204Pb206/X[i,'Pb206U238']
        J[3,3] <- 1/X[i,'Pb206U238']
        out$x <- c(U238Pb206,Pb207Pb206,Pb204Pb206)
        out$cov <- J %*% E %*% t(J)
    } else if (format == 5){
        out$x <- X[i,labels]
        out$cov <- cor2cov3(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],
                            X[i,'errPb204Pb206'],X[i,'rXY'],X[i,'rXZ'],X[i,'rYZ'])
    } else if (format == 6){
        U238Pb206 <- 1/X[i,'Pb206U238']
        Pb207Pb206 <- X[i,'Pb207Pb206']
        Pb204Pb206 <- X[i,'Pb204Pb206']
        out$x <- c(U238Pb206,Pb207Pb206,Pb204Pb206)
        out$cov <- matrix(0,3,3)
        out$cov[1,1] <- (U238Pb206*X[i,'errPb206U238']/X[i,'Pb206U238'])^2
        out$cov[2,2] <- X[i,'errPb207Pb206']^2
        out$cov[3,3] <- X[i,'errPb204Pb206']^2
        out$cov[1,2] <- get_cov_76_86(X[i,'Pb207Pb206'],X[i,'errPb207Pb206'],
                                      X[i,'Pb206U238'],X[i,'errPb206U238'],
                                      X[i,'Pb207U235'],X[i,'errPb207U235'])
        out$cov[1,3] <- get_cov_46_86(X[i,'Pb204Pb206'],X[i,'errPb204Pb206'],
                                      X[i,'Pb206U238'],X[i,'errPb206U238'],
                                      X[i,'Pb204U238'],X[i,'errPb204U238'])
        out$cov[2,3] <- get_cov_46_76(X[i,'Pb204Pb206'],X[i,'errPb204Pb206'],
                                      X[i,'Pb207Pb206'],X[i,'errPb207Pb206'],
                                      X[i,'Pb204Pb207'],X[i,'errPb204Pb207'])
        out$cov[2,1] <- out$cov[1,2]
        out$cov[3,1] <- out$cov[1,3]
        out$cov[3,2] <- out$cov[2,3]
    } else if (format == 7){
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/X[i,'Pb206U238']
        Pb207Pb206 <- X[i,'Pb207U235']/(X[i,'Pb206U238']*U238U235)
        Pb208Pb206 <- X[i,'Pb208Th232']*X[i,'Th232U238']/X[i,'Pb206U238']
        Th232U238 <- X[i,'Th232U238']
        E <- cor2cov4(X[i,'errPb207U235'],X[i,'errPb206U238'],
                      X[i,'errPb208Th232'],X[i,'errTh232U238'],
                      X[i,'rXY'],X[i,'rXZ'],X[i,'rXW'],
                      X[i,'rYZ'],X[i,'rYW'],X[i,'rZW'])
        J <- matrix(0,4,4)
        J[1,2] <- -U238Pb206/X[i,'Pb206U238']
        J[2,1] <- 1/(X[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/X[i,'Pb206U238']
        J[3,2] <- -Pb208Pb206/X[i,'Pb206U238']
        J[3,3] <- X[i,'Th232U238']/X[i,'Pb206U238']
        J[3,4] <- X[i,'Pb208Th232']/X[i,'Pb206U238']
        J[4,4] <- 1
        out$x <- c(U238Pb206,Pb207Pb206,Pb208Pb206,Th232U238)
        out$cov <- J %*% E %*% t(J)
    } else if (format == 8){
        out$x <- X[i,labels]
        out$cov <- cor2cov4(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],
                            X[i,'errPb208Pb206'],X[i,'errTh232U238'],
                            X[i,'rXY'],X[i,'rXZ'],X[i,'rXW'],
                            X[i,'rYZ'],X[i,'rYW'],X[i,'rZW'])
    } else if (format == 9){
        out$x <- X[i,labels]
        out$cov <- cor2cov2(X[i,'errU238Pb206'],X[i,'errPb204Pb206'],X[i,'rXY'])
    } else if (format == 10){
        out$x <- X[i,labels]
        out$cov <- cor2cov2(X[i,'errU235Pb207'],X[i,'errPb204Pb207'],X[i,'rXY'])
    } else if (format == 11){
        out$x <- X[i,labels]
        out$cov <- cor2cov3(X[i,'errU238Pb206'],X[i,'errPb208Pb206'],
                            X[i,'errTh232U238'],X[i,'rXY'],X[i,'rXZ'],X[i,'rYZ'])
    } else if (format == 12){
        out$x <- X[i,labels]
        out$cov <- cor2cov3(X[i,'errU235Pb207'],X[i,'errPb208Pb207'],
                            X[i,'errTh232U238'],X[i,'rXY'],X[i,'rXZ'],X[i,'rYZ'])
    } else if (format == 85){
        out$x <- X[i,labels]
        out$cov <- cor2cov3(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],
                            X[i,'errPb208Pb206'],X[i,'rXY'],X[i,'rXZ'],X[i,'rYZ'])
    } else if (format == 119){
        out$x <- X[i,labels]
        out$cov <- cor2cov2(X[i,'errU238Pb206'],X[i,'errPb208Pb206'],X[i,'rXY'])
    } else if (format == 1210){
        out$x <- X[i,labels]
        out$cov <- cor2cov2(X[i,'errU235Pb207'],X[i,'errPb208Pb207'],X[i,'rXY'])
    }
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    class(out) <- "terawasserburg"
    out
}
get_UPb_isochron_ratios_20x <- function(x,i=NULL){
    if (x$format%in%c(4,5,6)){
        labels <- c('U238Pb206','Pb204Pb206','U235Pb207','Pb204Pb207')
    } else if (x$format==9){
        labels <- c('U238Pb206','Pb204Pb206')
    } else if (x$format==10){
        labels <- c('U235Pb207','Pb204Pb207')
    } else if (x$format==85){
        labels <- c('U238Pb206','Pb208Pb206','U235Pb207','Pb208Pb207')
    } else if (x$format==119){
        labels <- c('U238Pb206','Pb208Pb206')
    } else if (x$format==1210){
        labels <- c('U235Pb207','Pb208Pb207')
    } else {
        stop('Invalid U-Pb format for get_UPb_isochron_ratios_20x.')
    }
    if (is.null(i)){
        ns <- length(x)
        out <- matrix(0,ns,length(labels))
        for (j in 1:ns){
            out[j,] <- get_UPb_isochron_ratios_20x(x,i=j)$x
        }
        colnames(out) <- labels
        return(out)
    }
    out <- list()
    if (x$format%in%c(9,10,119,1210)){
        out$x <- x$x[i,c(1,3)]
        out$cov <- cor2cov2(sX=x$x[i,2],sY=x$x[i,4],rXY=x$x[i,'rXY'])
    } else { # formats 4, 5, and 6
        Pbx6label <- ifelse(x$format<11,'Pb204Pb206','Pb208Pb206')
        U <- iratio('U238U235')[1]
        tw <- tera_wasserburg(x,i) # 38/06, 07/06 and 0x/06
        U8Pb6 <- tw$x['U238Pb206']
        Pbx6 <- tw$x[Pbx6label]
        U5Pb7 <- tw$x['U238Pb206']/(U*tw$x['Pb207Pb206'])
        Pbx7 <- tw$x[Pbx6label]/tw$x['Pb207Pb206']
        J <- matrix(0,4,3)
        J[1,1] <- 1
        J[2,3] <- 1
        J[3,1] <- 1/(U*tw$x['Pb207Pb206'])
        J[3,2] <- -U5Pb7/tw$x['Pb207Pb206']
        J[4,2] <- -Pbx7/tw$x['Pb207Pb206']
        J[4,3] <- 1/tw$x['Pb207Pb206']
        out$x <- c(U8Pb6,Pbx6,U5Pb7,Pbx7)
        out$cov <- J %*% tw$cov %*% t(J)
    }
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    out
}
get_UPb_isochron_ratios_208 <- function(x,i=NULL,tt=0){
    if (x$format%in%c(7,8)){
    labels <- c('U238Pb206','Pb208cPb206','U235Pb207',
                'Pb208cPb207','Th232U238','Th232Pb208',
                'Pb206cPb208','Pb207cPb208')
    } else if (x$format==11){
        labels <- c('U238Pb206','Pb208cPb206',
                    'Th232Pb208','Pb206cPb208')
    } else if (x$format==12){
        labels <- c('U235Pb207','Pb208cPb207',
                    'Th232Pb208','Pb207cPb208')
    } else {
        stop('Invalid format for get_UPb_isochron_ratios_208')
    }
    McL <- mclean(tt,d=x$d[i])
    if (is.null(i)){
        ns <- length(x)
        out <- matrix(0,ns,length(labels))
        for (j in 1:ns){
            out[j,] <- get_UPb_isochron_ratios_208(x,i=j,tt=tt)$x
        }
        colnames(out) <- labels
        return(out)
    } else if (x$format%in%c(7,8)){
        l2 <- settings('lambda','Th232')[1]
        U85 <- iratio('U238U235')[1]
        tw <- tera_wasserburg(x,i) # 38/06, 07/06, 08/06, 32/38
        U8Pb6 <- tw$x['U238Pb206']
        Pb8c6 <- tw$x['Pb208Pb206'] -
            tw$x['Th232U238']*tw$x['U238Pb206']*McL$Pb208Th232
        U5Pb7 <- tw$x['U238Pb206']/(U85*tw$x['Pb207Pb206'])
        Pb8c7 <- Pb8c6/tw$x['Pb207Pb206']
        Th2Pb8 <- tw$x['Th232U238']*tw$x['U238Pb206']/tw$x['Pb208Pb206']
        Pb6c8 <- (1 - McL$Pb206U238*tw$x['U238Pb206'])/tw$x['Pb208Pb206']
        Pb7c8 <- tw$x['Pb207Pb206']/tw$x['Pb208Pb206'] -
            McL$Pb207U235*tw$x['U238Pb206']/(U85*tw$x['Pb208Pb206'])
        J <- matrix(0,8,4)
        J[1,1] <- 1
        J[2,1] <- -tw$x['Th232U238']*(exp(l2*tt)-1)
        J[2,3] <- 1
        J[2,4] <- -tw$x['U238Pb206']*(exp(l2*tt)-1)
        J[3,1] <- 1/(U85*tw$x['Pb207Pb206'])
        J[3,2] <- -U5Pb7/tw$x['Pb207Pb206']
        J[4,1] <- J[2,1]/tw$x['Pb207Pb206']
        J[4,2] <- -Pb8c7/tw$x['Pb207Pb206']
        J[4,3] <- J[2,3]/tw$x['Pb207Pb206']
        J[4,4] <- J[2,4]/tw$x['Pb207Pb206']
        J[5,4] <- 1
        J[6,1] <- tw$x['Th232U238']/tw$x['Pb208Pb206']
        J[6,3] <- -Th2Pb8/tw$x['Pb208Pb206']
        J[6,4] <- tw$x['U238Pb206']/tw$x['Pb208Pb206']
        J[7,1] <- -McL$Pb206U238/tw$x['Pb208Pb206']
        J[7,3] <- -Pb6c8/tw$x['Pb208Pb206']
        J[8,1] <- -McL$Pb207U235/(U85*tw$x['Pb208Pb206'])
        J[8,2] <- 1/tw$x['Pb208Pb206']
        J[8,3] <- -Pb7c8/tw$x['Pb208Pb206']
        out <- list()
        out$x <- c(U8Pb6,Pb8c6,U5Pb7,Pb8c7,
                   tw$x['Th232U238'],Th2Pb8,Pb6c8,Pb7c8)
        E <- tw$cov
    } else if (x$format==11){
        l2 <- settings('lambda','Th232')[1]
        U8Pb6 <- x$x[i,'U238Pb206']
        Pb8c6 <- x$x[i,'Pb208Pb206'] -
            x$x[i,'Th232U238']*x$x[i,'U238Pb206']*(exp(l2*tt)-1)
        Th2Pb8 <- x$x[i,'Th232U238']*x$x[i,'U238Pb206']/x$x[i,'Pb208Pb206']
        Pb6c8 <- (1 - McL$Pb206U238*x$x[i,'U238Pb206'])/x$x[i,'Pb208Pb206']
        J <- matrix(0,4,3)
        J[1,1] <- J[2,2] <- 1
        J[2,1] <- -x$x[i,'Th232U238']*(exp(l2*tt)-1)
        J[2,3] <- -x$x[i,'U238Pb206']*(exp(l2*tt)-1)
        J[3,1] <- x$x[i,'Th232U238']/x$x[i,'Pb208Pb206']
        J[3,2] <- -Th2Pb8/x$x[i,'Pb208Pb206']
        J[3,3] <- x$x[i,'U238Pb206']/x$x[i,'Pb208Pb206']
        J[4,1] <- -McL$Pb206U238/x$x[i,'Pb208Pb206']
        J[4,2] <- -Pb6c8/x$x[i,'Pb208Pb206']
        out <- list()
        out$x <- c(U8Pb6,Pb8c6,Th2Pb8,Pb6c8)
        E <- cor2cov3(sX=x$x[i,'errU238Pb206'],
                      sY=x$x[i,'errPb208Pb206'],
                      sZ=x$x[i,'errTh232U238'],
                      rXY=x$x[i,'rXY'],
                      rXZ=x$x[i,'rXZ'],
                      rYZ=x$x[i,'rYZ'])
    } else if (x$format==12){
        U85 <- iratio('U238U235')[1]
        l2 <- settings('lambda','Th232')[1]
        U5Pb7 <- x$x[i,'U235Pb207']
        Pb8c7 <- x$x[i,'Pb208Pb207'] -
            x$x[i,'Th232U238']*U85*x$x[i,'U235Pb207']*(exp(l2*tt)-1)
        Th2Pb8 <- x$x[i,'Th232U238']*U85*x$x[i,'U235Pb207']/x$x[i,'Pb208Pb207']
        Pb7c8 <- (1 - McL$Pb207U235*x$x[i,'U235Pb207'])/x$x[i,'Pb208Pb207']
        J <- matrix(0,4,3)
        J[1,1] <- J[2,2] <- 1
        J[2,1] <- -x$x[i,'Th232U238']*U85*(exp(l2*tt)-1)
        J[2,3] <- -x$x[i,'U235Pb207']*U85*(exp(l2*tt)-1)
        J[3,1] <- x$x[i,'Th232U238']*U85/x$x[i,'Pb208Pb207']
        J[3,2] <- -Th2Pb8/x$x[i,'Pb208Pb207']
        J[3,3] <- U85*x$x[i,'U235Pb207']/x$x[i,'Pb208Pb207']
        J[4,1] <- -McL$Pb207U235/x$x[i,'Pb208Pb207']
        J[4,2] <- -Pb7c8/x$x[i,'Pb208Pb207']
        out <- list()
        out$x <- c(U5Pb7,Pb8c7,Th2Pb8,Pb7c8)
        E <- cor2cov3(sX=x$x[i,'errU235Pb207'],
                      sY=x$x[i,'errPb208Pb207'],
                      sZ=x$x[i,'errTh232U238'],
                      rXY=x$x[i,'rXY'],
                      rXZ=x$x[i,'rXZ'],
                      rYZ=x$x[i,'rYZ'])
    } else {
        stop('Invalid U-Pb format for get_UPb_isochron_ratios_208')
    }
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    out
}

w2tw <- function(w,format){
    if (format < 4){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206','rXY')
    } else if (format < 7){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206',
                    'Pb206Pb206','errPb204Pb206',
                    'rXY','rXZ','rYZ')
    } else if (format < 9){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206',
                    'Pb208Pb206','errPb208Pb206',
                    'Th232U238','errTh232U238',
                    'rXY','rXZ','rXW',
                    'rYZ','rYW','rZW')
    } else if (format==9){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb204Pb206','errPb204Pb206','rXY')
    } else if (format==10){
        cnames <- c('U235Pb207','errU235Pb207',
                    'Pb204Pb207','errPb204Pb207','rXY')
    } else if (format==11){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb208Pb206','errPb208Pb206',
                    'Th232U238','errTh232U238',
                    'rXY','rXZ','rYZ')
    } else if (format==12){
        cnames <- c('U235Pb207','errU235Pb207',
                    'Pb208Pb207','errPb208Pb207',
                    'Th232U238','errTh232U238',
                    'rXY','rXZ','rYZ')
    } else if (format==85){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206',
                    'Pb208Pb206','errPb208Pb206',
                    'rXY','rXZ','rYZ')
    } else if (format==119){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb208Pb206','errPb208Pb206','rXY')
    } else if (format==1210){
        cnames <- c('U235Pb207','errU235Pb207',
                    'Pb208Pb207','errPb208Pb207','rXY')
    } else {
        stop('Invalid input format.')
    }
    ns <- nrow(w)
    out <- w*0
    colnames(out) <- cnames
    for (i in 1:ns){
        tw <- tera_wasserburg(x=w,i=i,format=format)
        out[i,] <- wtw_helper(x=tw$x,covmat=tw$cov,cnames=cnames)
    }
    out
}

tw2w <- function(tw,format){
    if (format < 4){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238','rXY')
    } else if (format < 7){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238',
                    'Pb204U238','errPb204U238',
                    'rXY','rXZ','rYZ')
    } else if (format < 9){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238',
                    'Pb208Th232','errPb208Th232',
                    'Th232U238','errTh232U238',
                    'rXY','rXZ','rXW',
                    'rYZ','rYW','rZW')
    } else if (format==9){
        cnames <- c('Pb206U238','errPb206U238',
                    'Pb204U238','errPb204U238','rXY')
    } else if (format==10){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb204U235','errPb204U235','rXY')
    } else if (format==11){
        cnames <- c('Pb206U238','errPb206U238',
                    'Pb208Th232','errPb208Th232',
                    'Th232U238','errTh232U238',
                    'rXY','rXZ','rYZ')
    } else if (format==12){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb208Th232','errPb208Th232',
                    'Th232U238','errTh232U238',
                    'rXY','rXZ','rYZ')
    } else if (format==85){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238',
                    'Pb208U238','errPb208U238',
                    'rXY','rXZ','rYZ')
    } else if (format==119){
        cnames <- c('Pb206U238','errPb206U238',
                    'Pb208U238','errPb208U238','rXY')
    } else if (format==1210){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb208U235','errPb208U235','rXY')
    } else {
        stop('Invalid input format.')
    }
    ns <- nrow(tw)
    out <- tw*0
    colnames(out) <- cnames
    for (i in 1:ns){
        w <- wetherill(x=tw,i=i,format=format)
        out[i,] <- wtw_helper(x=w$x,covmat=w$cov,cnames=cnames)
    }
    out
}

wtw_helper <- function(x,covmat,cnames){
    nc <- length(cnames)
    if (any(is.na(x))){
        out <- rep(NA,nc)
        names(out) <- cnames
    } else {
        out <- rep(0,nc)
        names(out) <- cnames
        err <- sqrt(diag(covmat))
        cormat <- 0*covmat
        pos <- which(diag(covmat)>0)
        cormat[pos,pos] <- stats::cov2cor(covmat[pos,pos])
        out[c(1,3)] <- x[1:2]
        out[c(2,4)] <- err[1:2]
        out['rXY'] <- cormat[1,2]
        if (nc>5){
            out[5] <- x[3]
            out[6] <- err[3]
            out['rXZ'] <- cormat[1,3]
            out['rYZ'] <- cormat[2,3]
        }
        if (nc>9){
            out[7] <- x[4]
            out[8] <- err[4]
            out['rXW'] <- cormat[1,4]
            out['rYW'] <- cormat[2,4]
            out['rZW'] <- cormat[3,4]
        }
    }
    out
}

age_to_wetherill_ratios <- function(tt,st=0,exterr=FALSE,d=diseq()){
    out <- list()
    labels <- c('Pb207U235','Pb206U238')
    l8 <- settings('lambda','U238')[1]
    l5 <- settings('lambda','U235')[1]
    McL <- mclean(tt=tt,d=d,exterr=exterr)
    out$x <- c(McL$Pb207U235,McL$Pb206U238)
    E <- matrix(0,3,3)
    diag(E) <- c(st,lambda('U235')[2],lambda('U238')[2])^2
    J <- matrix(0,2,3)
    J[1,1] <- McL$dPb207U235dt
    J[2,1] <- McL$dPb206U238dt
    if (exterr){
        J[1,2] <- McL$dPb207U235dl35
        J[2,3] <- McL$dPb206U238dl38
    }
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}
age_to_terawasserburg_ratios <- function(tt,st=0,exterr=FALSE,d=diseq()){
    out <- list()
    labels <- c('U238Pb206','Pb207Pb206')
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    tt <- check_zero_UPb(tt)
    McL <- mclean(tt=tt,d=d,exterr=exterr)
    U238Pb206 <- 1/McL$Pb206U238
    Pb207Pb206 <- McL$Pb207U235/(U*McL$Pb206U238)
    d75dt <- McL$dPb207U235dt
    d68dt <- McL$dPb206U238dt
    d75dl35 <- McL$dPb207U235dl35
    d68dl38 <- McL$dPb206U238dl38
    out$x <- c(U238Pb206,Pb207Pb206)
    E <- matrix(0,4,4)
    diag(E) <- c(st,lambda('U235')[2],lambda('U238')[2],iratio('U238U235')[2])^2
    J <- matrix(0,2,4)
    J[1,1] <- -d68dt/McL$Pb206U238^2
    J[2,1] <- (d75dt*McL$Pb206U238-McL$Pb207U235*d68dt)/(U*McL$Pb206U238^2)
    if (exterr){
        J[1,3] <- -d68dl38/McL$Pb206U238^2
        J[2,2] <- d75dl35/(U*McL$Pb206U238)
        J[2,3] <- -Pb207Pb206*d68dl38/McL$Pb206U238
        J[2,4] <- -Pb207Pb206/U
    }
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}
age_to_cottle_ratios <- function(tt,st=0,exterr=FALSE,d=diseq(),option=1){
    out <- list()
    l2 <- settings('lambda','Th232')[1]
    McL <- mclean(tt=tt,d=d,exterr=exterr)
    Pb8Th2 <- exp(l2*tt)-1
    E <- matrix(0,3,3)
    J <- matrix(0,2,3)
    J[2,1] <- l2*exp(l2*tt)
    if (exterr) J[2,3] <- tt*exp(l2*tt)
    if (option==2){
        labels <- c('Pb207U235','Pb208Th232')
        l5 <- settings('lambda','U235')[1]
        Pb7U5 <- McL$Pb207U235
        out$x <- c(Pb7U5,Pb8Th2)
        diag(E) <- c(st,lambda('U235')[2],lambda('Th232')[2])^2
        J[1,1] <- McL$dPb207U235dt
        if (exterr) J[1,2] <- McL$dPb207U235dl35
    } else {
        labels <- c('Pb206U238','Pb208Th232')
        l8 <- settings('lambda','U238')[1]
        Pb6U8 <- McL$Pb206U238
        out$x <- c(Pb6U8,Pb8Th2)
        diag(E) <- c(st,lambda('U238')[2],lambda('Th232')[2])^2
        J[1,1] <- McL$dPb206U238dt
        if (exterr) J[1,2] <- McL$dPb206U238dl38
    }
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

age_to_Pb207U235_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('75','s[75]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb207U235_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l5 <- lambda('U235')[1]
        McL <- mclean(tt=tt,d=d,exterr=exterr)
        R <- McL$Pb207U235
        J <- McL$dPb207U235dt
        if (exterr){
            J <- c(J,McL$dPb207U235dl35)
            sl5 <- lambda('U235')[2]
            E <- diag(c(st,sl5))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            R.err <- abs(J*st)
        }
        out <- cbind(R,R.err)
        colnames(out) <- c('75','s[75]')
    }
    out
}
age_to_U235Pb207_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('57','s[57]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_U235Pb207_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l5 <- lambda('U235')[1]
        McL <- mclean(tt=tt,d=d,exterr=exterr)
        R <- 1/McL$Pb207U235
        if (exterr){
            J <- -c(McL$dPb207U235dt,McL$dPb207U235l35)/McL$Pb207U235^2
            sl5 <- lambda('U235')[2]
            E <- diag(c(st,sl5))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            J <- -McL$dPb207U235dt/McL$Pb207U235^2
            R.err <- abs(J*st)
        }
        out <- cbind(R,R.err)
        colnames(out) <- c('57','s[57]')
    }
    out
}
age_to_Pb206U238_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('68','s[68]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb206U238_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l8 <- lambda('U238')[1]
        McL <- mclean(tt=tt,d=d,exterr=exterr)
        R <- McL$Pb206U238
        J <- McL$dPb206U238dt
        if (exterr){
            J <- c(J,McL$dPb206U238dl38)
            sl8 <- lambda('U238')[2]
            E <- diag(c(st,sl8))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            R.err <- abs(J*st)
        }
        out <- cbind(R,R.err)
        colnames(out) <- c('68','s[68]')
    }
    out
}
age_to_U238Pb206_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('86','s[86]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_U238Pb206_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l8 <- lambda('U238')[1]
        tt <- check_zero_UPb(tt)
        McL <- mclean(tt=tt,d=d,exterr=exterr)
        R <- 1/McL$Pb206U238
        if (exterr){
            J <- -c(McL$dPb206U238dt,McL$dPb206U238l38)/McL$Pb206U238^2
            sl8 <- lambda('U238')[2]
            E <- diag(c(st,sl8))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            J <- -McL$dPb206U238dt/McL$Pb206U238^2
            R.err <- abs(J*st)
        }
        R.err <- abs(J*st)
        out <- cbind(R,R.err)
        colnames(out) <- c('86','s[86]')
    }
    out
}
age_to_Pb207Pb206_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('76','s[76]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb207Pb206_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l5 <- lambda('U235')[1]
        l8 <- lambda('U238')[1]
        tt <- check_zero_UPb(tt)
        U <- iratio('U238U235')[1]
        McL <- mclean(tt=tt,d=d,exterr=exterr)
        R <- (1/U)*McL$Pb207U235/McL$Pb206U238
        J <- (1/U)*(McL$dPb207U235dt*McL$Pb206U238 -
                    McL$Pb207U235*McL$dPb206U238dt)/McL$Pb206U238^2
        if (exterr){
            sl5 <- lambda('U235')[2]
            sl8 <- lambda('U238')[2]
            d76dl35 <- (1/U)*McL$dPb207U235dl35/McL$Pb206U238
            d76dl38 <- -(1/U)*McL$Pb207U235*McL$dPb206U238dl38/McL$Pb206U238^2
            J <- c(J,d76dl35,d76dl38)
            E <- diag(c(st,sl5,sl8))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            R.err <- abs(J*st)
        }
        out <- cbind(R,R.err)
        colnames(out) <- c('76','s[76]')
    }
    out
}
age_to_Pb208Th232_ratio <- function(tt,st=0,exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('82','s[82]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb208Th232_ratio(tt[i],st=sti)
        }
    } else {
        l2 <- lambda('Th232')[1]
        R <- exp(l2*tt)-1
        J <- l2*exp(l2*tt)
        if (exterr){
            sl2 <- lambda('Th232')[2]
            J <- c(J,tt*exp(l2*tt))
            E <- diag(c(st,sl2))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            R.err <- abs(J*st)
        }
        out <- cbind(R,R.err)
        colnames(out) <- c('82','s[82]')
    }
    out    
}

check_zero_UPb <- function(tt){
    smallnum <- 2*.Machine$double.neg.eps/lambda('U238')[1]
    if (length(tt)>1){
        pos <- which(tt>0)
        out <- tt
        out[!pos] <- smallnum
    } else {
        out <- max(tt,smallnum)
    }
    out
}

get_Pb207U235_ratios <- function(x,exterr=FALSE){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb207U235','errPb207U235')
    colnames(out) <- labels
    if (x$format %in% c(1,3,4,6,7)){
        out <- subset(x$x,select=labels)
    } else if (x$format %in% c(2,5,8,85)){
        R <- iratio('U238U235')[1]
        sR <- iratio('U238U235')[2]
        X <- x$x[,'U238Pb206']
        sX <- x$x[,'errU238Pb206']
        Y <- x$x[,'Pb207Pb206']
        sY <- x$x[,'errPb207Pb206']
        covXY <- x$x[,'rXY']*sX*sY
        out[,'Pb207U235'] <- R*Y/X
        relerr2 <- (sX/X)^2 -2*covXY/(X*Y) + (sY/Y)^2
        if (exterr) relerr2 <- relerr2 + (sR/R)^2
        out[,'errPb207U235'] <- sqrt(relerr2)*out[,'Pb207U235']
    } else if (x$format %in% c(10,12,1210)){
        out[,'Pb207U235'] <- 1/x$x[,'U235Pb207']
        out[,'errPb207U235'] <-
            out[,'Pb207U235']*x$x[,'errU235Pb207']/x$x[,'U235Pb207']
    } else {
        stop('Invalid U-Pb format for get_Pb207U235_ratios')
    }
    out
}
get_Pb206U238_ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb206U238','errPb206U238')
    colnames(out) <- labels
    if (x$format %in% c(1,3,4,6,7)){
        out <- subset(x$x,select=labels)
    } else if (x$format %in% c(2,5,8,9,11,85,119)){
        out[,'Pb206U238'] <- 1/x$x[,'U238Pb206']
        out[,'errPb206U238'] <- out[,'Pb206U238']*
            x$x[,'errU238Pb206']/x$x[,'U238Pb206']
    } else {
        stop('Invalid U-Pb format for get_Pb206U238_ratios')
    }
    out
}
get_U238Pb206_ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('U238Pb206','errU238Pb206')
    colnames(out) <- labels
    if (x$format %in% c(1,3,4,6,7)){
        out[,'U238Pb206'] <- 1/x$x[,'Pb206U238']
        out[,'errU238Pb206'] <- out[,'U238Pb206']*
            x$x[,'errPb206U238']/x$x[,'Pb206U238']
    } else if (x$format %in% c(2,5,8,9,11,85,119)){
        out <- subset(x$x,select=labels)
    } else {
        stop('Invalid U-Pb format for get_U238Pb206_ratios')
    }
    out
}
get_Pb207Pb206_ratios <- function(x,exterr=FALSE){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb207Pb206','errPb207Pb206')
    colnames(out) <- labels
    if (x$format %in% c(1,4,7)){
        R <- iratio('U238U235')[1]
        sR <- iratio('U238U235')[2]
        X <- x$x[,'Pb207U235']
        sX <- x$x[,'errPb207U235']
        Y <- x$x[,'Pb206U238']
        sY <- x$x[,'errPb206U238']
        rXY <- x$x[,'rXY']
        covXY <- rXY*sX*sY
        out[,'Pb207Pb206'] <- X/(R*Y)
        relerr2 <- (sX/X)^2 - 2*covXY/(X*Y) + (sY/Y)^2
        if (exterr) relerr2 <- relerr2 + (sR/R)^2
        out[,'errPb207Pb206'] <- sqrt(relerr2)*out[,'Pb207Pb206']
    } else if (x$format %in% c(2,3,5,6,8,85)){
        out <- subset(x$x,select=labels)
    } else {
        stop('Invalid U-Pb format for get_Pb207Pb206_ratios')
    }
    out
}
get_Pb208Th232_ratios <- function(x){
    labels <- c('Pb208Th232','errPb208Th232')
    if (x$format == 7){
        out <- x$x[,labels,drop=FALSE]
    } else if (x$format%in%c(8,11)){
        ns <- length(x)
        out <- matrix(0,ns,2)
        colnames(out) <- labels
        out[,'Pb208Th232'] <-
            x$x[,'Pb208Pb206']/(x$x[,'U238Pb206']*x$x[,'Th232U238'])
        J1 <- -out[,'Pb208Th232']/x$x[,'U238Pb206']
        J2 <- 1/(x$x[,'U238Pb206']*x$x[,'Th232U238'])
        J3 <- -out[,'Pb208Th232']/x$x[,'Th232U238']
        E11 <- x$x[,'errU238Pb206']^2
        E22 <- x$x[,'errPb208Pb206']^2
        E33 <- x$x[,'errTh232U238']^2
        if (x$format == 8){
            E12 <- x$x[,'rXZ']*x$x[,'errU238Pb206']*x$x[,'errPb208Pb206']
            E13 <- x$x[,'rXW']*x$x[,'errU238Pb206']*x$x[,'errTh232U238']
            E23 <- x$x[,'rZW']*x$x[,'errPb208Pb206']*x$x[,'errTh232U238']
        } else {
            E12 <- x$x[,'rXY']*x$x[,'errU238Pb206']*x$x[,'errPb208Pb206']
            E13 <- x$x[,'rXZ']*x$x[,'errU238Pb206']*x$x[,'errTh232U238']
            E23 <- x$x[,'rYZ']*x$x[,'errPb208Pb206']*x$x[,'errTh232U238']
        }
        out[,'errPb208Th232'] <-
            sqrt(errorprop1x3(J1,J2,J3,E11,E22,E33,E12,E13,E23))
    } else if (x$format == 12){
        ns <- length(x)
        out <- matrix(0,ns,2)
        colnames(out) <- labels
        U85 <- iratio('U238U235')[1]
        out[,'Pb208Th232'] <-
            x$x[,'Pb208Pb207']/(x$x[,'U235Pb207']*U85*x$x[,'Th232U238'])
        J1 <- -out[,'Pb208Th232']/x$x[,'U235Pb207']
        J2 <- 1/(x$x[,'U235Pb207']*U85*x$x[,'Th232U238'])
        J3 <- -out[,'Pb208Th232']/x$x[,'Th232U238']
        E11 <- x$x[,'errU235Pb207']^2
        E22 <- x$x[,'errPb208Pb207']^2
        E33 <- x$x[,'errTh232U238']^2
        E12 <- x$x[,'rXY']*x$x[,'errU235Pb207']*x$x[,'errPb208Pb207']
        E13 <- x$x[,'rXZ']*x$x[,'errU235Pb207']*x$x[,'errTh232U238']
        E23 <- x$x[,'rYZ']*x$x[,'errPb208Pb207']*x$x[,'errTh232U238']
    } else {
        stop('Wrong input format: no Pb208 or Th232 present in this dataset.')
    }
    out
}
get_Pb208Pb206_ratios <- function(x){
    labels <- c('Pb208Pb206','errPb208Pb206')
    if (x$format == 7){
        ns <- length(x)
        out <- matrix(0,ns,2)
        colnames(out) <- labels
        out[,'Pb208Pb206'] <- x$x[,'Pb208Th232']*x$x[,'Th232U238']/x$x[,'Pb206U238']
        J1 <- -out[,'Pb208Pb206']/x$x[,'Pb206U238']
        J2 <- x$x[,'Th232U238']/x$x[,'Pb206U238']
        J3 <- x$x[,'Pb208Th232']/x$x[,'Pb206U238']
        E11 <- x$x[,'errPb206U238']^2
        E22 <- x$x[,'errPb208Th232']^2
        E33 <- x$x[,'Th232U238']^2
        E12 <- x$x[,'rXZ']*x$x[,'Pb206U238']*x$x[,'Pb208Th232']
        E13 <- x$x[,'rXW']*x$x[,'Pb206U238']*x$x[,'Th232U238']
        E23 <- x$x[,'rZW']*x$x[,'Pb208Th232']*x$x[,'Th232U238']
        out[,'errPb208Th232'] <-
            sqrt(errorprop1x3(J1,J2,J3,E11,E22,E33,E12,E13,E23))
    } else if (x$format %in% c(8,11,85,119)){
        out <- x$x[,labels]
    } else {
        stop('Wrong input format: no Pb208 present in this dataset.')
    }
    out
}

#' Compute Pb207/U235 age(s)
#' @param x either an object of class \code{UPb} or a Pb207/U235-ratio
#' @noRd
get_Pb207U235_age <- function(x,...){
    UseMethod("get_Pb207U235_age",x)
}
#' @param sx the standard error of \code{x}
#' @param exterr propagate the decay constant uncertainty?
#' @param d an object of class \code{diseq}
#' @noRd
get_Pb207U235_age.default <- function(x,sx=0,exterr=FALSE,d=diseq(),...){
    ns <- length(x)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('t75','s[t75]')
        for (i in 1:ns){
            if (length(sx) < ns) sxi <- sx[1]
            else sxi <- sx[i]
            out[i,] <- get_Pb207U235_age(x[i],sxi,exterr=exterr,d=d[i],...)
        }
    } else {
        if (is.na(x)) return(c('t75'=NA,'s[t75]'=NA))
        l5 <- lambda('U235')[1]
        t.init <- ifelse(x>-1,log(1+x)/l5,0)
        E <- matrix(0,3,3)
        J <- matrix(0,1,3)
        E[1,1] <- sx^2
        E[2:3,2:3] <- getEl('U235')
        if (d$equilibrium | d$PaU$option<1){
            t.75 <- t.init
            J[1,1] <- 1/(l5*(1+x))                       # dt/dx
            if (exterr & x>-1) J[1,2] <- log(1+x)/l5^2   # dt/dl5
        } else { # apply a disequilibrium correction
            t.75 <- stats::optimise(diseq_75_misfit,x=x,d=d,
                                    interval=t.init*c(.5,2))$minimum
            McL <- mclean(tt=t.75,d=d,exterr=exterr)    
            J[1,1] <- 1/McL$dPb207U235dt                 # dt/dx
            J[1,2] <- -McL$dPb207U235dl35/McL$dPb207U235dt # dt/dl35
            J[1,3] <- -McL$dPb207U235dl31/McL$dPb207U235dt # dt/dl31
        }
        st.75 <- sqrt(J %*% E %*% t(J))
        out <- c(t.75,st.75)
        names(out) <- c('t75','s[t75]')
    }
    out
}
#' @param i an aliquot number. If \code{NULL}, returns all ages
#' @noRd
get_Pb207U235_age.UPb <- function(x,i=NULL,exterr=FALSE,...){
    r75 <- get_Pb207U235_ratios(x)
    if (is.null(i)){
        ns <- length(x)
        out <- matrix(0,ns,2)
        for (j in 1:ns){
            out[j,] <- get_Pb207U235_age.UPb(x,i=j,exterr=exterr,...)
        }
    } else {
        out <- get_Pb207U235_age(r75[i,'Pb207U235'],r75[i,'errPb207U235'],
                                 exterr=exterr,d=x$d[i],...)
    }
    out
}
#' @noRd
get_Pb207U235_age.wetherill <- function(x,exterr=FALSE,...){
    i <- 'Pb207U235'
    r75 <- x$x[i]
    sr75 <- sqrt(x$cov[i,i])
    get_Pb207U235_age(r75,sr75,exterr=exterr,d=x$d,...)
}

#' Compute Pb206/U238 age(s)
#' @param x either an object of class \code{UPb} or a Pb206/U238-ratio
#' @noRd
get_Pb206U238_age <- function(x,...){
    UseMethod("get_Pb206U238_age",x)
}
#' @param sx the standard error of \code{x}
#' @param exterr propagate the decay constant uncertainty?
#' @param d an object of class \code{diseq}
#' @param bayes if \code{TRUE}, computes the posterior distribution
#'     for the disequilibrium correction
#' @param plot if \code{TRUE}, plots the aforementioned posterior
#'     distribution
#' @noRd
get_Pb206U238_age.default <- function(x,sx=0,exterr=FALSE,d=diseq(),
                                      bayes=FALSE,plot=FALSE,...){
    ns <- length(x)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('t68','s[t68]')
        for (i in 1:ns){
            if (length(sx) < ns) sxi <- sx[1]
            else sxi <- sx[i]
            out[i,] <- get_Pb206U238_age(x[i],sxi,exterr=exterr,d=d[i])
        }
    } else {
        if (is.na(x)) return(c('t68'=NA,'s[t68]'=NA))
        E <- matrix(0,5,5)
        E[1,1] <- sx^2
        E[2:5,2:5] <- getEl('U238')
        J <- matrix(0,1,5)
        l8 <- lambda('U238')[1]
        if (x>-1) t.init <- log(1+x)/l8 else t.init <- 0
        if (d$equilibrium){
            t.68 <- t.init
            J[1,1] <- 1/(l8*(1+x))                       # dt/dx
            if (exterr & x>-1) J[1,2] <- log(1+x)/l8^2   # dt/dl38
        } else { # apply a disequilibrium correction
            if (measured_disequilibrium(d)) tlim <- c(0,meas_diseq_maxt(d))
            else tlim <- c(0,4600)
            t.68 <- tryCatch({
                stats::optimise(diseq_68_misfit,interval=tlim,x=x,d=d)$minimum
            }, error = function(error_condition) {
                warning("Failed to solve the Pb206U238 disequilibrium equation.")
                t.init
            })
            McL <- mclean(tt=t.68,d=d,exterr=exterr)
            J[1,1] <- 1/McL$dPb206U238dt
            if (exterr){
                J[1,2] <- -McL$dPb206U238dl38/McL$dPb206U238dt # dt/dl38
                J[1,3] <- -McL$dPb206U238dl34/McL$dPb206U238dt # dt/dl34
                J[1,4] <- -McL$dPb206U238dl30/McL$dPb206U238dt # dt/dl30
                J[1,5] <- -McL$dPb206U238dl26/McL$dPb206U238dt # dt/dl26
            }
        }
        st.68 <- sqrt(J %*% E %*% t(J))
        out <- c(t.68,st.68)
        names(out) <- c('t68','s[t68]')
    }
    out
}
#' @param i an aliquote number. If \code{NULL}, returns all ages
#' @noRd
get_Pb206U238_age.UPb <- function(x,i=NULL,exterr=FALSE,...){
    r68 <- get_Pb206U238_ratios(x)
    if (is.null(i)){
        ns <- length(x)
        out <- matrix(0,ns,2)
        for (j in 1:ns){
            out[j,] <- get_Pb206U238_age.UPb(x,i=j,exterr=exterr,...)
        }
    } else {
        out <- get_Pb206U238_age(r68[i,'Pb206U238'],
                                 r68[i,'errPb206U238'],
                                 exterr=exterr,d=x$d[i],...)
    }
    out
}
#' @noRd
get_Pb206U238_age.wetherill <- function(x,exterr=FALSE,...){
    i <- 'Pb206U238'
    r68 <- x$x[i]
    sr68 <- sqrt(x$cov[i,i])
    get_Pb206U238_age(r68,sr68,exterr=exterr,d=x$d,...)
}
#' @noRd
get_Pb206U238_age.terawasserburg <- function(x,exterr=FALSE,...){
    i <- 'U238Pb206'
    r86 <- x$x[i]
    r68 <- 1/r86
    sr68 <- sqrt(x$cov[i,i])/r86
    get_Pb206U238_age(r68,sr68,exterr=exterr,d=x$d,...)
}

twslope <- function(tt=0,d=diseq()){
    McL <- mclean(tt=tt,d=d)
    McL$dPb207U235dt/McL$dPb206U238dt
}

#' Compute Pb207/Pb206 age(s)
#' @param x either an object of class \code{UPb} or a Pb207/Pb206-ratio
#' @noRd
get_Pb207Pb206_age <- function(x,...){
    UseMethod("get_Pb207Pb206_age",x)
}
#' @param sx the standard error of \code{x}
#' @param exterr propagate the decay constant uncertainty?
#' @param d an object of class \code{diseq}
#' @param t.68 a Pb206/U238-age to initialise the calculation
#' @noRd
get_Pb207Pb206_age.default <- function(x,sx=0,exterr=FALSE,d=diseq(),t.68=NULL,...){
    ns <- length(x)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('t76','s[t76]')
        for (i in 1:ns){
            if (length(sx) < ns) sxi <- sx[1]
            else sxi <- sx[i]
            out[i,] <- get_Pb207Pb206_age(x[i],sxi,exterr=exterr,d=d[i],t.68=t.68)
        }
    } else {
        if (is.na(x)) return(c('t68'=NA,'s[t68]'=NA))
        interval <- c(1/10000,10000)
        if (!d$equilibrium & !is.null(t.68)){
            midpoint <- stats::optimise(twslope,d=d,interval=interval)$minimum
            if (t.68<midpoint){
                interval[2] <- midpoint
            } else {
                interval[1] <- midpoint
            }
        }
        if (measured_disequilibrium(d)) interval[2] <- meas_diseq_maxt(d)
        t.76 <- stats::optimise(get.76.misfit,x=x,d=d,interval=interval)$minimum
        McL <- mclean(tt=t.76,d=d,exterr=exterr)
        E <- matrix(0,9,9)
        E[1,1] <- sx^2
        E[2,2] <- iratio('U238U235')[2]^2
        E[3:9,3:9] <- getEl()
        J <- matrix(0,1,9)
        J[1,1] <- -1/McL$dPb207Pb206dt                # dt/dx
        J[1,2] <- McL$dPb207Pb206dU/McL$dPb207Pb206dt   # dt/dU
        J[1,3] <- McL$dPb207Pb206dl38/McL$dPb207Pb206dt # dt/dl38
        J[1,4] <- McL$dPb207Pb206dl35/McL$dPb207Pb206dt # dt/dl35
        J[1,5] <- McL$dPb207Pb206dl34/McL$dPb207Pb206dt # dt/dl34
        J[1,7] <- McL$dPb207Pb206dl31/McL$dPb207Pb206dt # dt/dl31
        J[1,8] <- McL$dPb207Pb206dl30/McL$dPb207Pb206dt # dt/dl30
        J[1,9] <- McL$dPb207Pb206dl26/McL$dPb207Pb206dt # dt/dl26
        st.76 <- sqrt( J %*% E %*% t(J) )
        out <- c(t.76,st.76)
        names(out) <- c('t76','s[t76]')
    }
    out
}
#' @param i and aliquot number
#' @noRd
get_Pb207Pb206_age.UPb <- function(x,i=1,exterr=FALSE,...){
    r76 <- get_Pb207Pb206_ratios(x)
    get_Pb207Pb206_age(r76[i,'Pb207Pb206'],r76[i,'errPb207Pb206'],
                       exterr=exterr,d=x$d[i],...)
}
#' @noRd
get_Pb207Pb206_age.wetherill <- function(x,exterr=FALSE,...){
    U <- iratio('U238U235')[1]
    r76 <- x$x['Pb207U235']/(U*x$x['Pb206U238'])
    J <- matrix(0,1,2)
    E <- x$cov
    J[1,1] <- 1/(U*x$x['Pb206U238'])                   # d76d75
    J[1,2] <- -x$x['Pb207U235']/(U*x$x['Pb206U238']^2) # d76d68
    sr76 <- J %*% E %*% t(J)
    get_Pb207Pb206_age(r76,sr76,exterr=exterr,d=x$d)
}
#' @noRd
get_Pb207Pb206_age.terawasserburg <- function(x,exterr=FALSE,...){
    r76 <- x$x['Pb207Pb206']
    sr76 <- x$cov['Pb207Pb206','Pb207Pb206']
    get_Pb207Pb206_age(r76,sr76,exterr=exterr,d=x$d)
}

#' Compute Pb208/Th232 age(s)
#' @param x either an object of class \code{UPb} or a Pb208/Th232-ratio
#' @noRd
get_Pb208Th232_age <- function(x,...){
    UseMethod("get_Pb208Th232_age",x)
}
#' @param sx the standard error of \code{x}
#' @param exterr propagate the decay constant uncertainty?
#' @noRd
get_Pb208Th232_age.default <- function(x,sx=0,exterr=FALSE,...){
    ns <- length(x)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('t82','s[t82]')
        for (i in 1:ns){
            if (length(sx) < ns) sxi <- sx[1]
            else sxi <- sx[i]
            out[i,] <- get_Pb208Th232_age(x[i],sxi,exterr=exterr)
        }
    } else {
        if (is.na(x)) return(c('t82'=NA,'s[t82]'=NA))
        l2 <- lambda('Th232')[1]
        sl2 <- lambda('Th232')[2]
        if (x>-1) t.init <- log(1+x)/l2 else t.init <- 0
        J <- matrix(0,1,2)
        t.82 <- t.init
        J[1,1] <- 1/(l2*(1+x))                       # dt/dx
        if (exterr & x>-1) J[1,2] <- log(1+x)/l2^2   # dt/dl2
        E <- matrix(0,2,2)
        E[1,1] <- sx^2
        E[2,2] <- sl2^2
        st.82 <- sqrt(J %*% E %*% t(J))
        out <- c(t.82,st.82)
        names(out) <- c('t82','s[t82]')
    }
    out
}
#' @param i an aliquot number
#' @noRd
get_Pb208Th232_age.UPb <- function(x,i=NULL,exterr=FALSE,...){
    r82 <- get_Pb208Th232_ratios(x)
    if (is.null(i)){
        ns <- length(x)
        out <- matrix(0,ns,2)
        for (j in 1:ns){
            out[j,] <- get_Pb208Th232_age.UPb(x,i=j,exterr=exterr,...)
        }
    } else {
        out <- get_Pb208Th232_age(r82[i,'Pb208Th232'],r82[i,'errPb208Th232'],
                                 exterr=exterr,...)
    }
    out
}

# x is an object of class \code{UPb}
# returns a matrix of 7/5, 6/8, 7/6
# and concordia_ages and their uncertainties.
UPb_age <- function(x,exterr=FALSE,i=NULL,conc=TRUE,omit4c=NULL,
                    discordance=discfilter(),common.Pb=0,...){
    if (discordance$option==0 | discordance$before) xd <- x
    else xd <- Pb0corr(x,option=common.Pb,omit4c=omit4c)
    if (common.Pb==0) X <- x
    else X <- Pb0corr(x,option=common.Pb,omit4c=omit4c)
    if (!is.null(i)){
        out <- UPb_age_helper(x=x,X=X,xd=xd,i=i,exterr=exterr,
                              conc=conc,discordance=discordance)
    } else {
        nn <- length(x)
        out <- NULL
        for (i in 1:nn){
            ti <- UPb_age_helper(x=x,X=X,xd=xd,i=i,exterr=exterr,
                                 conc=conc,discordance=discordance)
            out <- rbind(out,ti)
        }
    }
    out
}

# x = raw data
# X = common Pb corrected data (if common.Pb>0)
# xd = data to be used for concordia_age calculation 
#      (raw if before==TRUE, common Pb corrected if before==FALSE)
UPb_age_helper <- function(x,X,xd,i=1,exterr=FALSE,
                           conc=TRUE,discordance=discfilter(),...){
    Xi <- subset(X,subset=((1:length(X))%in%i))
    do68 <- do75 <- do76 <- do82 <- FALSE 
    if (x$format<7 || x$format==85){
        do68 <- do75 <- do76 <- TRUE
        labels <- c('t.75','s[t.75]','t.68','s[t.68]','t.76','s[t.76]')
    } else if (x$format<9){
        do68 <- do75 <- do76 <- do82 <- TRUE
        labels <- c('t.75','s[t.75]','t.68','s[t.68]',
                    't.76','s[t.76]','t.82','s[t.82]')
    } else if (x$format%in%c(9,119)){
        labels <- c('t.68','s[t.68]')
        do68 <- TRUE
    } else if (x$format%in%c(10,1210)){
        labels <- c('t.75','s[t.75]')
        do75 <- TRUE
    } else if (x$format==11){
        do68 <- do82 <- TRUE
        labels <- c('t.68','s[t.68]','t.82','s[t.82]')
    } else if (x$format==12){
        do75 <- do82 <- TRUE
        labels <- c('t.75','s[t.75]','t.82','s[t.82]')
    } else {
        stop('Invalid U-Pb format')
    }
    t.75 <- t.68 <- t.76 <- t.82 <- t.conc <- dif <- pval <- NULL
    if (do75){
        t.75 <- get_Pb207U235_age(Xi,exterr=exterr)
    }
    if (do68){
        t.68 <- get_Pb206U238_age(Xi,exterr=exterr)
    }
    if (do76){
        t.76 <- get_Pb207Pb206_age(Xi,exterr=exterr,
                                   t.68=subset(t.68,select=1))
    }
    if (do82) {
        t.82 <- get_Pb208Th232_age(Xi,exterr=exterr)
    }
    if (x$format<9 || x$format==85){
        if (conc){
            labels <- c(labels,'t.conc','s[t.conc]')
            t.conc <- concordia_age(x=Xi,i=1,exterr=exterr)$age
        }
        if (discordance$option>0){
            xdi <- subset(xd,subset=((1:length(xd))%in%i))
        }
        if (discordance$option%in%c(1,'t',2,'r',3,'sk',4,'a',5,'c')){
            labels <- c(labels,'disc')
            xi <- subset(x,subset=((1:length(x))%in%i))
            dif <- discordance(x=xi,X=xdi,option=discordance$option)
        } else if (discordance$option%in%c(6,'p')){
            labels <- c(labels,'p[conc]')
            SS.concordance <-
                LL_concordia_age(pars=t.conc[1],
                                 cc=wetherill(xdi,i=1),
                                 mswd=TRUE,exterr=exterr,d=xdi$d)
            pval <- 1-stats::pchisq(SS.concordance,1)
        }
    }
    out <- c(t.75,t.68,t.76,t.82,t.conc,dif,pval)
    names(out) <- labels
    out
}

getEl <- function(parent){
    out <- matrix(0,7,7)
    out[1,1] <- lambda('U238')[2]^2
    out[2,2] <- lambda('U235')[2]^2
    out[3,3] <- (lambda('U234')[2]*1000)^2
    out[4,4] <- lambda('Th232')[2]^2
    out[5,5] <- (lambda('Pa231')[2]*1000)^2
    out[6,6] <- (lambda('Th230')[2]*1000)^2
    out[7,7] <- (lambda('Ra226')[2]*1000)^2
    rcnames <- c('U238','U235','U234','Th232','Pa231','Th230','Ra226')
    rownames(out) <- colnames(out) <- rcnames
    if (missing(parent)){
        # do nothing
    } else if (identical(parent,'U238')){
        rcnames <- c('U238','U234','Th230','Ra226')
    } else if (identical(parent,'U235')){
        rcnames <- c('U235','Pa231')
    } else if (identical(parent,'Th232')){
        rcnames <- c('Th232')
    } else {
        # do nothing
    }
    out[rcnames,rcnames]
}
