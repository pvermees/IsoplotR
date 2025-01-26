#' @title Set up a discordance filter
#' 
#' @description Define a discordance cutoff to filter U--Pb data.
#' 
#' @details The most reliable U--Pb age constraints are obtained from
#'     (zircon) grains whose \eqn{^{206}}Pb/\eqn{^{238}}U and
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb ages are statistically
#'     indistinguishable from each other. U--Pb compositions that
#'     fulfil this requirements are called `concordant'. Those that
#'     violate it are called `discordant'. The discordance of the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U and \eqn{^{207}}Pb/\eqn{^{206}}Pb
#'     systems can be defined in five different ways. By setting a
#'     cutoff for any of these criteria, U--Pb data can be filtered
#'     for data quality.
#'
#' @param option one of five options:
#'
#' \code{0}: do not apply a discordance filter
#' 
#' \code{1} or \code{'t'}: the absolute age difference (Ma) between
#' the \eqn{^{206}}Pb/\eqn{^{238}}U and \eqn{^{207}}Pb/\eqn{^{206}}Pb
#' ages.
#' 
#' \code{2} or \code{'r'}: the relative age difference (\%) between the
#' \eqn{^{206}}Pb/\eqn{^{238}}U and \eqn{^{207}}Pb/\eqn{^{206}}Pb ages.
#'
#' \code{3} or \code{'sk'}: percentage of common Pb measured along a
#' mixing line connecting the measured composition and the
#' Stacey-Kramers mantle composition in Tera-Wasserburg space.
#'
#' \code{4} or \code{'a'}: logratio distance (\%) measured along a
#' perpendicular line connecting Tera-Wasserburg concordia and the
#' measured composition.
#'
#' \code{5} or \code{'c'}: logratio distance (\%) measured along a
#' line connecting the measured composition and the corresponding
#' single grain concordia_age composition.
#'
#' Further details in Vermeesch (2021).
#'
#' @param before logical flag indicating whether the discordance
#'     filter should be applied before (\code{TRUE}) or after
#'     (\code{FALSE}) the common-Pb correction.
#'
#' @param cutoff a two-element vector with the minimum (negative) and
#'     maximum (positive) allowed discordance. Default values vary
#'     between the different options. To view them, enter
#'     \code{discfilter(option)} at the \code{R} command line.
#' 
#' @return a list with the input parameters. Default values for
#'     \code{cutoff} are
#'
#' \code{c(-48,140)} if \code{option=='t'};
#'
#' \code{c(-5,15)} if \code{option=='r'};
#' 
#' \code{c(-0.36,0.96)} if \code{option=='sk'};
#'
#' \code{c(-1.6,4.7)} if \code{option=='a'}; and
#' 
#' \code{c(-2,5.8)} if \code{option=='c'}.
#' 
#' @seealso \code{\link{cad}}, \code{\link{kde}},
#'     \code{\link{radialplot}}
#' @references Vermeesch (2021) ``On the treatment of discordant data
#'     in detrital zircon U--Pb geochronology'', Geochronology.
#' @examples
#' dscf <- discfilter(option='c',before=TRUE,cutoff=c(-1,1))
#' weightedmean(x=examples$UPb,exterr=FALSE,sigdig=2,
#'              cutoff.disc=dscf,common.Pb=3)
#' 
#' @export
discfilter <- function(option=0,before=TRUE,cutoff){
    out <- list()
    out$option <- option
    out$before <- before
    if (missing(cutoff)){
        if (option%in%c(1,'t')) cutoff <- c(-48,140)
        else if (option%in%c(2,'r')) cutoff <- c(-5,15)
        else if (option%in%c(3,'sk')) cutoff <- c(-0.36,0.96)
        else if (option%in%c(4,'a')) cutoff <- c(-1.6,4.7)
        else if (option%in%c(5,'c')) cutoff <- c(-2,5.8)
        else cutoff <- c(-Inf,Inf)
    }
    out$cutoff <- cutoff
    class(out) <- 'discfilter'
    out
}

filter_UPb_ages <- function(x,type=5,cutoff.76=1100,exterr=FALSE,
                            cutoff.disc=discfilter(),common.Pb=0,omit4c=NULL){
    if (x$format>8){ # override type if necessary:
        if (x$format==9) type <- 2
        else if (x$format==10) type <- 1
        else if (x$format==11 && type!=6) type <- 2
        else if (x$format==12 && type!=6) type <- 1
        else if (x$format==85 && type==6) type <- 2
        else if (x$format==119) type <- 2
        else if (x$format==1210) type <- 1
        else stop('Invalid UPb format')
    }
    tt <- UPb_age(x,exterr=exterr,conc=(type==5),omit4c=omit4c,
                  common.Pb=common.Pb,discordance=cutoff.disc)
    if (cutoff.disc$option==0){
        is.concordant <- rep(TRUE,length(x))
    } else {
        dcol <- which(colnames(tt)%in%c('disc','p[conc]'))
        is.concordant <- (tt[,dcol]>cutoff.disc$cutoff[1]) &
            (tt[,dcol]<cutoff.disc$cutoff[2])
    }
    if (!any(is.concordant)){
        stop(paste0('There are no concordant grains in this sample.',
                    'Try adjusting the discordance limits OR ',
                    'apply a common-Pb correction OR ',
                    '(if you have already applied a common-Pb correction), ',
                    'apply the discordance filter before the ',
                    'common-Pb correction.'))
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

# x: raw data, xd: common Pb corrected data (or not)
discordance <- function(x,xd=x,i=NULL,option=4){
    if (is.null(i)){
        ns <- length(x)
        dif <- rep(0,ns)
        for (j in 1:ns){
            dif[j] <- discordance(x,xd=x,i=j,option=option)
        }
    } else {
        xi <- subset(x,subset=((1:length(x))%in%i))
        xdi <- subset(xd,subset=((1:length(xd))%in%i))
        t.68 <- get_Pb206U238_age(xdi)[1]
        t.76 <- get_Pb207Pb206_age(xdi,t.68=t.68)[1]
        if (option%in%c(5,'c')){
            t.conc <- concordia_age(x=xd,i=i)$age[1]
        }
        if (option%in%c(1,'t')){
            dif <- t.76-t.68
        } else if (option%in%c(2,'r')){
            dif <- (1-t.68/t.76)*100
        } else if (option%in%c(3,'sk')){
            x.corr <- Pb0corr(xi,option=3)
            U8Pb6.raw <- get_U238Pb206_ratios(xi)[,'U238Pb206']
            U8Pb6.corr <- get_U238Pb206_ratios(x.corr)[,'U238Pb206']
            dif <- (1-U8Pb6.raw/U8Pb6.corr)*100
        } else if (option%in%c(4,'a')){
            U8Pb6 <- get_U238Pb206_ratios(xdi)[,'U238Pb206']
            Pb76 <- get_Pb207Pb206_ratios(xdi)[,'Pb207Pb206']
            r86.76 <- age_to_U238Pb206_ratio(t.76)[,1]
            r76.68 <- age_to_Pb207Pb206_ratio(t.68)[,1]
            DX <- (log(U8Pb6) - log(r86.76))/sqrt(2)
            DY <- (log(Pb76) - log(r76.68))/sqrt(2/3)
            dif <- 100*DX*sin(atan(DY/DX))
        } else if (option%in%c(5,'c')){
            U8Pb6 <- get_U238Pb206_ratios(xdi)[,'U238Pb206']
            Pb76 <- get_Pb207Pb206_ratios(xdi)[,'Pb207Pb206']
            c86 <- age_to_U238Pb206_ratio(t.conc)[,1]
            c76 <- age_to_Pb207Pb206_ratio(t.conc)[,1]
            dx <- (log(U8Pb6) - log(c86))/sqrt(2)
            dy <- (log((Pb76^2)/U8Pb6) - log((c76^2)/c86))/sqrt(6)
            dif <- 100*sign(t.76-t.68)*sqrt(dx^2+dy^2)
        } else {
            dif <- 0
        }
    }
    dif
}

is.discordant <- function(x,xd=x,cutoff.disc=discfilter()){
    disc <- discordance(x=x,xd=xd,option=cutoff.disc$option)
    out <- which(disc < cutoff.disc$cutoff[1] | disc > cutoff.disc$cutoff[2])
    out
}
