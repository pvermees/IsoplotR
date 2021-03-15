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
#' \code{2} or \code{'r'}: the relative age difference (%) between the
#' \eqn{^{206}}Pb/\eqn{^{238}}U and \eqn{^{207}}Pb/\eqn{^{206}}Pb ages.
#'
#' \code{3} or \code{'sk'}: percentage of common Pb measured along a
#' mixing line connecting the measured composition and the
#' Stacey-Kramers mantle composition in Tera-Wasserburg space.
#'
#' \code{4} or \code{'a'}: logratio distance (%) measured along a
#' perpendicular line connecting Tera-Wasserburg concordia and the
#' measured composition.
#'
#' \code{5} or \code{'c'}: logratio distance (%) measured along a line
#' connecting the measured composition and the corresponding single
#' grain concordia age composition.
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
#' \code{c(-50,140)} if \code{option=='t'};
#'
#' \code{c(-5,15)} if \code{option=='r'};
#' 
#' \code{c(-0.3,1)} if \code{option=='sk'};
#'
#' \code{c(-2,6)} if \code{option=='a'}; and
#' 
#' \code{c(-2,7)} if \code{option=='c'}.
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
        if (option%in%c(1,'t')) cutoff <- c(-50,140)
        else if (option%in%c(2,'r')) cutoff <- c(-5,15)
        else if (option%in%c(3,'sk')) cutoff <- c(-0.3,1)
        else if (option%in%c(4,'a')) cutoff <- c(-2,6)
        else if (option%in%c(5,'c')) cutoff <- c(-2,7)
        else cutoff <- c(-Inf,Inf)
    }
    out$cutoff <- cutoff
    class(out) <- 'discfilter'
    out
}

filter.UPb.ages <- function(x,type=5,cutoff.76=1100,exterr=FALSE,
                            cutoff.disc=discfilter(),common.Pb=0,omit4c=NULL){
    tt <- UPb.age(x,exterr=exterr,conc=(type==5),omit4c=omit4c,
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

# x: raw data, X: common Pb corrected data (or not)
discordance <- function(x,X,tt=NULL,option=4){
    t.68 <- get.Pb206U238.age(X)[1]
    t.76 <- get.Pb207Pb206.age(X,t.68=t.68)[1]
    if (option%in%c(5,'c')){
        t.conc <- concordia.age(x=X,i=1)$age[1]
    }
    if (option%in%c(1,'t')){
        dif <- t.76-t.68
    } else if (option%in%c(2,'r')){
        dif <- (1-t.68/t.76)*100
    } else if (option%in%c(3,'sk')){
        x.corr <- Pb0corr(x,option=3)
        U8Pb6.raw <- get.U238Pb206.ratios(x)[,'U238Pb206']
        U8Pb6.corr <- get.U238Pb206.ratios(x.corr)[,'U238Pb206']
        dif <- (1-U8Pb6.raw/U8Pb6.corr)*100
    } else if (option%in%c(4,'a')){
        x76 <- age_to_U238Pb206_ratio(t.76)[,1]
        y68 <- age_to_Pb207Pb206_ratio(t.68)[,1]
        U8Pb6 <- get.U238Pb206.ratios(X)[,'U238Pb206']
        Pb76 <- get.Pb207Pb206.ratios(X)[,'Pb207Pb206']
        DX <- log(U8Pb6) - log(x76)
        DY <- log(Pb76) - log(y68)
        dif <- 100*DX*sin(atan(DY/DX))
    } else if (option%in%c(5,'c')){
        xc <- age_to_U238Pb206_ratio(t.conc)[,1]
        U8Pb6 <- get.U238Pb206.ratios(X)[,'U238Pb206']
        dx <- log(xc) - log(U8Pb6)
        yc <- age_to_Pb207Pb206_ratio(t.conc)[,1]
        Pb76 <- get.Pb207Pb206.ratios(X)[,'Pb207Pb206']
        dy <- log(yc) - log(Pb76)
        dif <- 100*sign(t.76-t.68)*sqrt(dx^2+dy^2)
    } else {
        stop('Invalid discordance filter option.')
    }
    dif
}
