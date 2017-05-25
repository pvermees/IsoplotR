#' Calculate isotopic ages
#'
#' Calculates  ages and propagates their analytical
#' uncertainties.
#'
#' @param x can be:
#' \itemize{
#' \item a scalar containing an isotopic ratio,
#'
#' \item a two element vector containing an isotopic ratio and its standard
#' error, or the spontaneous and induced track densities \code{Ns} and
#' \code{Ni} (if \code{method='fissiontracks'}),
#'
#' \item a four element vector containing \code{Ar40Ar39},
#' \code{s[Ar40Ar39]}, \code{J}, \code{s[J]},
#'
#' \item a six element vector containing \code{U}, \code{s[U]}, \code{Th},
#' \code{s[Th]}, \code{He} and \code{s[He]},
#'
#' \item an eight element vector containing \code{U}, \code{s[U]},
#' \code{Th}, \code{s[Th]}, \code{He}, \code{s[He]}, \code{Sm} and
#' \code{s[Sm]}
#'
#' \item a six element vector containing \code{Rb}, \code{s[Rb]},
#' \code{Sr}, \code{s[Sr]}, \code{Sr87Sr86}, and \code{s[Sr87Sr86]}
#'
#' \item a six element vector containing \code{Re}, \code{s[Re]},
#' \code{Os}, \code{s[Os]}, \code{Os187Os188}, and \code{s[Os187Os188]}
#'
#' \item a six element vector containing \code{Sm}, \code{s[Sm]},
#' \code{Nd}, \code{s[Nd]}, \code{Nd143Nd144}, and \code{s[Nd144Nd143]}
#'
#' \item a six element vector containing \code{Lu}, \code{s[Lu]},
#' \code{Hf}, \code{s[Hf]}, \code{Hf176Hf177}, and \code{s[Hf176Hf177]}
#' }
#'
#' OR
#'
#' \itemize{
#' \item an object of class \code{UPb}, \code{ArAr},
#' \code{RbSr}, \code{SmNd}, \code{ReOs}, \code{LuHf} \code{UThHe} or
#' \code{fissiontracks}.
#' }
#'
#' @param method one of either \code{'U238-Pb206'}, \code{'U235-Pb207'},
#'     \code{'Pb207-Pb206'}, \code{'Ar-Ar'}, \code{'Re-Os'},
#'     \code{'Sm-Nd'}, \code{'Rb-Sr'}, \code{'Lu-Hf'}, \code{'U-Th-He'} or
#'     \code{'fissiontracks'}
#' 
#' @param exterr propagate the external (decay constant and
#'     calibration factor) uncertainties?
#' 
#' @param i (optional) index of a particular aliquot
#' 
#' @param ... additional arguments
#' 
#' @rdname age
#' @export
age <- function(x,...){ UseMethod("age",x) }
#' @rdname age
#' @export
age.default <- function(x,method='U238-Pb206',exterr=TRUE,J=c(NA,NA),
                        zeta=c(NA,NA),rhoD=c(NA,NA),...){
    if (length(x)==1) X <- c(x,0)
    else X <- x[1:2]
    if (identical(method,'U235-Pb207')){
        out <- get.Pb207U235.age(X[1],X[2],exterr)
    } else if (identical(method,'U238-Pb206')){
        out <- get.Pb206U238.age(X[1],X[2],exterr)
    } else if (identical(method,'Pb206-Pb207')){
        out <- get.Pb207Pb206.age(X[1],X[2],exterr)
    } else if (identical(method,'Ar-Ar')){
        out <- get.ArAr.age(X[1],X[2],X[3],X[4],exterr)
    } else if (identical(method,'Re-Os')){
        out <- get.ReOs.age(X[1],X[2],exterr)
    } else if (identical(method,'Rb-Sr')){
        out <- get.RbSr.age(X[1],X[2],exterr)
    } else if (identical(method,'Sm-Nd')){
        out <- get.SmNd.age(X[1],X[2],exterr)
    } else if (identical(method,'Lu-Hf')){
        out <- get.LuHf.age(X[1],X[2],exterr)
    } else if (identical(method,'U-Th-He')){
        if (length(x)==6)
            out <- get.UThHe.age(X[1],X[2],X[3],X[4],X[5],X[6])
        else if (length(x)==8)
            out <- get.UThHe.age(X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8])
    } else if (identical(method,'fissiontracks')){
        out <- get.EDM.age(X[1],X[2],zeta,rhoD)
    } else {
        out <- x
    }
    out
}
#' @param type scalar flag indicating whether each U-Pb analysis
#'     should be considered separately (\code{type=1}), a
#'     concordia age should be calculated from all U-Pb analyses
#'     together (\code{type=2}), or a discordia line should be
#'     fit through all the U-Pb analyses (\code{type=3}).
#' 
#' @param wetherill logical flag to indicate whether the data should
#'     be evaluated in Wetherill (\code{TRUE}) or Tera-Wasserburg
#'     (\code{FALSE}) space.  This option is only used when
#'     \code{type=2}
#' 
#' @param sigdig number of significant digits for the uncertainty
#'     estimate (only used if \code{type=1}, \code{isochron=FALSE}
#'     or \code{central=FALSE}).
#'
#' @return 
#'
#' \enumerate{
#'
#' \item if \code{x} is a scalar or a vector, returns the age using
#' the geochronometer given by \code{method} and its standard error.
#'
#' \item if \code{x} has class \code{UPb} and \code{type=1},
#' returns a table with the following columns: \code{t.75},
#' \code{err[t.75]}, \code{t.68}, \code{err[t.68]}, \code{t.76},
#' \code{err[t.76]}, \code{t.conc}, \code{err[t.conc]}, containing the
#' \eqn{^{207}}Pb/\eqn{^{235}}U-age and standard error, the
#' \eqn{^{206}}Pb/\eqn{^{238}}U-age and standard error, the
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb-age and standard error, and the
#' concordia age and standard error, respectively.
#'  
#' \item if \code{x} has class \code{UPb} and \code{type=2},
#' returns a list with the following items:
#'
#' \describe{
#' \item{x}{ a named vector with the (weighted mean) U-Pb composition }
#' 
#' \item{cov}{ the covariance matrix of the (mean) U-Pb composition }
#' 
#' \item{age}{ the concordia age (in Ma) }
#' 
#' \item{age.err}{ the standard error of the concordia age }
#' 
#' \item{mswd}{ a list with two items (\code{equivalence} and
#' \code{concordance}) containing the MSWD (Mean of the Squared
#' Weighted Deviates, a.k.a the reduced Chi-squared statistic outside
#' of geochronology) of isotopic equivalence and age concordance,
#' respectively. }
#' 
#' \item{p.value}{ a list with two items (\code{equivalence} and
#' \code{concordance}) containing the p-value of the Chi-square test
#' for isotopic equivalence and age concordance, respectively. }
#' }
#' 
#' \item if \code{x} has class \code{UPb} and \code{type=3},
#' returns a list with the following items:
#'
#' \describe{
#' \item{x}{ a two element vector with the upper and lower intercept
#' ages (if \code{wetherill=TRUE}) or the lower intercept age and
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb intercept (for Tera-Wasserburg) }
#' 
#' \item{cov}{ the covariance matrix of the elements in \code{x} }
#' }
#'
#' \item if \code{x} has class \code{ArAr}, \code{RbSr}, \code{SmNd},
#' \code{ReOs}, \code{LuHf} and \code{isochron=FALSE}, returns a table of Ar-Ar,
#' Rb-Sr, Sm-Nd, Re-Os or Lu-Hf ages and standard errors.
#'
#' \item if \code{x} has class \code{ArAr}, \code{RbSr}, \code{SmNd},
#' \code{ReOs} or \code{LuHf} and \code{isochron=TRUE}, returns a list with the
#' following items:
#'
#' \describe{
#'
#' \item{a}{ the intercept of the straight line fit and its standard
#' error. }
#' 
#' \item{b}{ the slope of the fit and its standard error. }
#' 
#' \item{y0}{ the atmospheric \eqn{^{40}}Ar/\eqn{^{36}}Ar or initial
#' \eqn{^{187}}Os/\eqn{^{188}}Os, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#' \eqn{^{143}}Nd/\eqn{^{144}}Nd or \eqn{^{176}}Hf/\eqn{^{177}}Hf
#' ratio and its standard error. }
#' 
#' \item{age}{ the \eqn{^{40}}Ar/\eqn{^{39}}Ar,
#' \eqn{^{187}}Os/\eqn{^{187}}Re, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#' \eqn{^{143}}Nd/\eqn{^{144}}Nd or \eqn{^{176}}Hf/\eqn{^{177}}Hf age
#' and its standard error. }
#' 
#' }
#' 
#' \item if \code{x} has class \code{UThHe} and \code{central=FALSE},
#' returns a table of U-Th-He ages and standard errors.
#' 
#' \item if \code{x} has class \code{UThHe} and \code{central=TRUE},
#' returns a list with the following items:
#'
#' \describe{
#'
#' \item{uvw}{ a three-element list with the weighted mean log[U/He],
#' log[Th/He] and log[Sm/He] compositions. }
#'
#' \item{covmat}{ a 3x3 covariance matrix for \code{uvw}}
#'
#' \item{mswd}{ the reduced Chi-square value for the
#' log[U/He]-log[Th/He] compositions. }
#'
#' \item{p.value}{ the p-value of concordance between the
#' log[U/He]-log[Th/He] compositions. }
#'
#' \item{age}{ two-element vector with the central age and its
#' standard error. }
#'
#' }
#'
#' \item if \code{x} has class \code{fissiontracks} and
#' \code{central=FALSE}, returns a table of fission track ages and
#' standard errors.
#' 
#' \item if \code{x} has class \code{fissiontracks} and
#' \code{central=TRUE}, returns a list with the following items:
#'
#' \describe{
#'
#' \item{mswd}{ the reduced Chi-square value for the fission track
#' ages. }
#'
#' \item{p.value}{ the p-value of concordance between the fission
#' track ages. }
#'
#' \item{age}{ a two-element vector with the central age and its
#' standard error. }
#'
#' \item{disp}{ the (over)dispersion of the single grain ages beyond
#' the formal analytical uncertainties. }
#' }
#' }
#' 
#' @examples
#' data(examples)
#' print(age(examples$UPb))
#' print(age(examples$UPb,type=1))
#' print(age(examples$UPb,type=2))
#' @rdname age
#' @export
age.UPb <- function(x,type=1,wetherill=TRUE,
                    exterr=TRUE,i=NA,sigdig=NA,...){
    if (type==1)
        out <- UPb.age(x,exterr=exterr,i=i,sigdig=sigdig,...)
    else if (type==2)
        out <- concordia.age(x,wetherill=TRUE,exterr=TRUE,...)
    else if (type==3)
        out <- discordia.age(x,wetherill=TRUE,exterr=TRUE,...)
    out
}
#' @param J two-element vector with the J-factor and its standard
#'     error.
#' @param isochron logical flag indicating whether each Ar-Ar analysis
#'     should be considered separately (\code{isochron=FALSE}) or an
#'     isochron age should be calculated from all Ar-Ar analyses
#'     together (\code{isochron=TRUE}).
#' @param i2i
#'     `isochron to intercept': calculates the initial (aka `inherited',
#'     `excess', or `common') \eqn{^{40}}Ar/\eqn{^{36}}Ar,
#'     \eqn{^{187}}Os/\eqn{^{188}}Os or \eqn{^{176}}Hf/\eqn{^{177}}Hf
#'     ratio from an isochron fit. Setting \code{i2i} to \code{FALSE}
#'     uses the default values stored in \code{settings('iratio',...)}
#' @rdname age
#' @export
age.ArAr <- function(x,isochron=FALSE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE,...)
    else out <- ArAr.age(x,exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
    out
}
#' @param central logical flag indicating whether each U-Th-He analysis
#'     should be considered separately (\code{central=FALSE}) or a
#'     central age should be calculated from all U-Th-He analyses
#'     together (\code{central=TRUE}).
#' @rdname age
#' @export
age.UThHe <- function(x,central=FALSE,i=NA,sigdig=NA,...){
    if (central) out <- central(x)
    else out <- UThHe.age(x,i=i,sigdig=sigdig)
    out
}
#' @param zeta two-element vector with the zeta-factor and its standard
#'     error.
#' @param rhoD two-element vector with the track density of the
#'     dosimeter glass and its standard error.
#' @rdname age
#' @export
age.fissiontracks <- function(x,central=FALSE,i=NA,sigdig=NA,exterr=TRUE,...){
    if (central) out <- central(x)
    else out <- fissiontrack.age(x,i=i,sigdig=sigdig,exterr=exterr)
    out
}
#' @rdname age
#' @export
age.ReOs <- function(x,isochron=TRUE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    age.PD(x,'Re187',isochron=isochron,i2i=i2i,exterr=exterr,i=i,sigdig=sigdig,...)
}
#' @rdname age
#' @export
age.SmNd <- function(x,isochron=TRUE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    age.PD(x,'Sm147',isochron=isochron,i2i=i2i,exterr=exterr,i=i,sigdig=sigdig,...)
}
#' @rdname age
#' @export
age.RbSr <- function(x,isochron=TRUE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    age.PD(x,'Rb87',isochron=isochron,i2i=i2i,exterr=exterr,i=i,sigdig=sigdig,...)
}
#' @rdname age
#' @export
age.LuHf <- function(x,isochron=TRUE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    age.PD(x,'Lu176',isochron=isochron,i2i=i2i,exterr=exterr,i=i,sigdig=sigdig,...)
}
age.PD <- function(x,nuclide,isochron=TRUE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE,sigdig=sigdig)
    else out <- PD.age(x,nuclide,exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
    out
}
# tt and st are the age and error (scalars produced by peakfit or weightedmean)
# calculated without taking into account the external errors
add.exterr <- function(x,tt,st,cutoff.76=1100,type=4){
    out <- c(0,0)
    if (hasClass(x,'UPb')){
        if (type==1){
            R <- age.to.Pb207U235.ratio(tt,st)
            out <- get.Pb207U235.age(R[,1],R[,2],exterr=TRUE)
        } else if (type==2 | (type==4 & (tt<cutoff.76)) | (type==5)){
            R <- age.to.Pb206U238.ratio(tt,st)
            out <- get.Pb206U238.age(R[,1],R[,2],exterr=TRUE)
        } else if (type==3 | (type==4 & (tt>=cutoff.76))){
            R <- age.to.Pb207Pb206.ratio(tt,st)
            out <- get.Pb207Pb206.age(R[,1],R[,2],exterr=TRUE)
        }
    } else if (hasClass(x,'ArAr')){
        R <- get.ArAr.ratio(tt,st,x$J[1],0,exterr=FALSE)
        out <- get.ArAr.age(R[1],R[2],x$J[1],x$J[2],exterr=TRUE)
    } else if (hasClass(x,'ReOs')){
        R <- get.ReOs.ratio(tt,st,exterr=FALSE)
        out <- get.ReOs.age(R[1],R[2],exterr=TRUE)
    } else if (hasClass(x,'SmNd')){
        R <- get.SmNd.ratio(tt,st,exterr=FALSE)
        out <- get.SmNd.age(R[1],R[2],exterr=TRUE)
    } else if (hasClass(x,'RbSr')){
        R <- get.RbSr.ratio(tt,st,exterr=FALSE)
        out <- get.RbSr.age(R[1],R[2],exterr=TRUE)
    } else if (hasClass(x,'LuHf')){
        R <- get.LuHf.ratio(tt,st,exterr=FALSE)
        out <- get.LuHf.age(R[1],R[2],exterr=TRUE)
    } else if (hasClass(x,'fissiontracks')){
        out[2] <- tt * sqrt( (x$zeta[2]/x$zeta[1])^2 + (st/tt)^2 )
    }
    out
}
