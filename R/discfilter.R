discfilter <- function(option='c',before=TRUE,cutoff){
    out <- list()
    out$option <- option
    out$before <- before
    if (missing(cutoff)){
        if (option%in%c(0,'t')) cutoff <- c(-5,50)
        else if (option%in%c(1,'r')) cutoff <- c(-5,15)
        else if (option%in%c(2,'sk')) cutoff <- c(-0.001,0.01)
        else if (option%in%c(3,'a')) cutoff <- c(-1,5)
        else if (option%in%c(4,'c')) cutoff <- c(-1,5)
        else stop('Invalid discordance filter option.')
    }
    out$cutoff <- cutoff
    class(out) <- 'discfilter'
    out
}

filter.UPb.ages <- function(x,type=5,cutoff.76=1100,exterr=FALSE,
                            cutoff.disc=discfilter(),common.Pb=0){
    if (is.null(cutoff.disc)){
        is.concordant <- rep(TRUE,length(x))
    } else {
        conc <- (type==5) | (cutoff.disc$option=='c')
        if (cutoff.disc$before | common.Pb==0){
            X <- x
            tt <- UPb.age(X,exterr=exterr,conc=conc)
        } else {
            X <- Pb0corr(x,option=common.Pb)
            tt <- UPb.age(X,exterr=exterr,conc=conc,common.Pb=common.Pb)
        }
        if (cutoff.disc$option%in%c(0,'t')){
            dif <- tt[,'t.76']-tt[,'t.68']
        } else if (cutoff.disc$option%in%c(1,'r')){
            dif <- (1-tt[,'t.68']/tt[,'t.76'])*100
        } else if (cutoff.disc$option%in%c(2,'sk')){
            U8Pb6.raw <- get.U238Pb206.ratios(x)[,'U238Pb206']
            U8Pb6.corr <- get.U238Pb206.ratios(X)[,'U238Pb206']
            dif <- (1-U8Pb6.raw/U8Pb6.corr)*100
        } else if (cutoff.disc$option%in%c(3,'a')){
            X1 <- age_to_U238Pb206_ratio(tt[,'t.76'])[,1]
            Y1 <- age_to_Pb207Pb206_ratio(tt[,'t.68'])[,1]
            U8Pb6 <- get.U238Pb206.ratios(X)[,'U238Pb206']
            Pb76 <- get.Pb207Pb206.ratios(X)[,'Pb207Pb206']
            DX <- log(U8Pb6) - log(X1)
            DY <- log(Pb76) - log(Y1)
            dif <- 100*DX*sin(atan(DY/DX))
        } else if (cutoff.disc$option%in%c(4,'c')){
            x1 <- age_to_U238Pb206_ratio(tt[,'t.conc'])[,1]
            U8Pb6 <- get.U238Pb206.ratios(X)[,'U238Pb206']
            dx <- log(x1) - log(U8Pb6)
            x2 <- age_to_Pb207Pb206_ratio(tt[,'t.conc'])[,1]
            Pb76 <- get.Pb207Pb206.ratios(X)[,'Pb207Pb206']
            dy <- log(x2) - log(Pb76)
            dif <- 100*sqrt(dx^2+dy^2)
        } else {
            stop('Invalid discordance filter option.')
        }
        is.concordant <- (dif>cutoff.disc$cutoff[1]) & (dif<cutoff.disc$cutoff[2])
    }
    if (!any(is.concordant)){
        stop(paste0('There are no concordant grains in this sample.',
                    'Try adjusting the discordance limits OR ',
                    'apply a common-Pb correction OR ',
                    '(if you have already applied a common-Pb correction), ',
                    'apply the discordance filter before the ',
                    'common-Pb correction.'))
    }
    if (!cutoff.disc$before & common.Pb!=0){
        tt <- UPb.age(X,exterr=exterr,conc=conc,common.Pb=common.Pb)
    }
    out <- matrix(NA,length(x),2)
    if (type==1){
        out[is.concordant,] <- tt[is.concordant,c('t.75','s[t.75]'),drop=FALSE]
    } else if (type==2){
        out[is.concordant,] <- tt[is.concordant,c('t.68','s[t.68]'),drop=FALSE]
    } else if (type==3){
        out[is.concordant,] <- tt[is.concordant,c('t.76','s[t.76]'),drop=FALSE]
    } else if (type==4){
        do.76 <- (tt[,'t.68']>cutoff.76)
        i.76 <- as.vector(which(do.76 & is.concordant))
        i.68 <- as.vector(which(!do.76 & is.concordant))
        out[i.76,] <- tt[i.76,c('t.76','s[t.76]'),drop=FALSE]
        out[i.68,] <- tt[i.68,c('t.68','s[t.68]'),drop=FALSE]
    } else if (type==5){
        out[is.concordant,] <- tt[is.concordant,c('t.conc','s[t.conc]'),drop=FALSE]
    } else if (type==6){
        out[is.concordant,] <- tt[is.concordant,c('t.82','s[t.82]'),drop=FALSE]
    }
    colnames(out) <- c('t','s[t]')
    out
}
