#' Calculate isotopic ages
#'
#' Calculates U-Pb, Pb-Pb, Ar-Ar, K-Ca, Re-Os, Sm-Nd, Rb-Sr, Lu-Hf, U-Th-He,
#' Th-U and fission track ages and propagates their analytical
#' uncertainties. Includes options for single grain, isochron and
#' concordia ages.
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
#' \item a two element vector containing \code{K40Ca40} and
#' \code{s[K40Ca40]},
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
#'
#' \item a five element vector containing \code{0/8}, \code{s[0/8]},
#' \code{4/8}, \code{s[4/8]} and \code{cov[0/8,4/8]}
#'
#' }
#'
#' OR
#'
#' \itemize{ \item an object of class \code{UPb}, \code{PbPb},
#' \code{ArAr}, \code{KCa}, \code{ThU}, \code{RbSr}, \code{SmNd},
#' \code{ReOs}, \code{LuHf}, \code{UThHe} or \code{fissiontracks}.  }
#'
#' @param method one of either \code{'U238-Pb206'},
#'     \code{'U235-Pb207'}, \code{'Pb207-Pb206'}, \code{'Ar-Ar'},
#'     \code{'K-Ca'}, \code{'Th-U'}, \code{'Re-Os'}, \code{'Sm-Nd'},
#'     \code{'Rb-Sr'}, \code{'Lu-Hf'}, \code{'U-Th-He'} or
#'     \code{'fissiontracks'}
#' 
#' @param exterr propagate the external (decay constant and
#'     calibration factor) uncertainties?
#' 
#' @param i (optional) index of a particular aliquot
#' 
#' @param U48 the initial \eqn{^{234}}U/\eqn{^{238}}U-activity ratio
#'     (only used if \code{method='U238-Pb206'}, \code{'U238-Pb206'}
#'     or \code{'Pb207-Pb206'}).
#' @param Th0U8 the initial \eqn{^{230}}Th/\eqn{^{238}}U-activity
#'     ratio (only used if \code{method='U238-Pb206'},
#'     \code{'U238-Pb206'} or \code{'Pb207-Pb206'}).
#' @param Ra6U8 the initial \eqn{^{226}}Ra/\eqn{^{238}}U-activity
#'     ratio (only used if \code{method='U238-Pb206'},
#'     \code{'U238-Pb206'} or \code{'Pb207-Pb206'}).
#' @param Pa1U5 the initial \eqn{^{231}}Pa/\eqn{^{235}}U-activity
#'     ratio (only used if \code{method='U238-Pb206'},
#'     \code{'U238-Pb206'} or \code{'Pb207-Pb206'}).
#'
#' @param ... additional arguments
#'
#' @rdname age
#' @export
age <- function(x,...){ UseMethod("age",x) }
#' @rdname age
#' @export
age.default <- function(x,method='U238-Pb206',exterr=TRUE,J=c(NA,NA),
                        zeta=c(NA,NA),rhoD=c(NA,NA),U48=1,Th0U8=1,Ra6U8=1,Pa1U5=1,...){
    d <- diseq(U48=U48,Th0U8=Th0U8,Ra6U8=Ra6U8,Pa1U5=Pa1U5)
    if (length(x)==1) x <- c(x,0)
    if (identical(method,'U235-Pb207')){
        out <- get.Pb207U235.age(x=x[1],sx=x[2],exterr=exterr,d=d)
    } else if (identical(method,'U238-Pb206')){
        out <- get.Pb206U238.age(x=x[1],sx=x[2],exterr=exterr,d=d)
    } else if (identical(method,'Pb207-Pb206')){
        out <- get.Pb207Pb206.age(x=x[1],sx=x[2],exterr,d=d)
    } else if (identical(method,'Ar-Ar')){
        out <- get.ArAr.age(Ar40Ar39=x[1],sAr40Ar39=x[2],
                            J=x[3],sJ=x[4],exterr=exterr)
    } else if (identical(method,'K-Ca')){
        out <- get.KCa.age(K40Ca40=x[1],sK40Ca40=x[2],exterr=exterr)
    } else if (identical(method,'Re-Os')){
        out <- get.ReOs.age(Os187Re187=x[1],sOs187Re187=x[2],exterr=exterr)
    } else if (identical(method,'Rb-Sr')){
        out <- get.RbSr.age(Rb87Sr86=x[1],sRb87Sr86=x[2],exterr)
    } else if (identical(method,'Sm-Nd')){
        out <- get.SmNd.age(Nd143Sm147=x[1],sNd143Sm147=x[2],exterr)
    } else if (identical(method,'Lu-Hf')){
        out <- get.LuHf.age(Hf176Lu176=x[1],sHf176Lu176=x[2],exterr)
    } else if (identical(method,'Th-U')){
        out <- get.ThU.age(Th230U238=x[1],sTh230U238=x[2],U234U238=x[3],
                           sU234U238=x[4],cov4808=x[5],exterr=exterr)
    } else if (identical(method,'U-Th-He') && length(x)==6){
        out <- get.UThHe.age(U=x[1],sU=x[2],Th=x[3],
                             sTh=x[4],He=x[5],sHe=x[6])
    } else if (identical(method,'U-Th-He') && length(x)==8){
        out <- get.UThHe.age(U=x[1],sU=x[2],Th=x[3],sTh=x[4],
                             He=x[5],sHe=x[6],Sm=x[7],sSm=x[8])
    } else if (identical(method,'fissiontracks')){
        out <- get.EDM.age(Ns=x[1],Ni=x[2],zeta=zeta,rhoD=rhoD)
    } else {
        out <- NA
    }
    out
}
#' @param type scalar flag indicating whether
#'
#' \code{1}: each U-Pb analysis should be considered separately,
#'
#' \code{2}: all the measurements should be combined to calculate a
#' concordia age,
#'
#' \code{3}: a discordia line should be fitted through all the U-Pb
#'     analyses using the maximum likelihood algorithm of Ludwig
#'     (1998), which assumes that the scatter of the data is solely
#'     due to the analytical uncertainties.
#'
#' \code{4}: a discordia line should be fitted ignoring the analytical
#' uncertainties.
#'
#' \code{5}: a discordia line should be fitted using a modified
#' maximum likelihood algorithm that accounts for overdispersion by
#' adding a geological (co)variance term.
#'
#' @param wetherill logical flag to indicate whether the data should
#'     be evaluated in Wetherill (\code{TRUE}) or Tera-Wasserburg
#'     (\code{FALSE}) space.  This option is only used when
#'     \code{type=2}
#' @param sigdig number of significant digits for the uncertainty
#'     estimate (only used if \code{type=1}, \code{isochron=FALSE} and
#'     \code{central=FALSE}).
#' @param common.Pb apply a common lead correction using one of three
#'     methods:
#'
#' \code{1}: use the Stacey-Kramer two-stage model to infer the initial
#' Pb-composition
#'
#' \code{2}: use the isochron intercept as the initial Pb-composition
#'
#' \code{3}: use the Pb-composition stored in
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')}
#'
#' @param show.p Show the p-value for concordance for each aliquot to
#'     the output table. Note: it would be unwise to use the p-value
#'     value as a concordance filter. Doing so would 'punish' high
#'     precision measurements, which are more likely to fail the
#'     Chi-square test than low precision measurements. The latter
#'     would therefore be 'rewarded' by such a criterion.
#'
#' @return
#' \enumerate{
#'
#' \item if \code{x} is a scalar or a vector, returns the age using
#' the geochronometer given by \code{method} and its standard error.
#'
#' \item if \code{x} has class \code{UPb} and \code{type=1}, returns a
#' table with the following columns: \code{t.75}, \code{err[t.75]},
#' \code{t.68}, \code{err[t.68]}, \code{t.76}, \code{err[t.76]},
#' \code{t.conc}, \code{err[t.conc]}, \code{err[p.conc]}, containing
#' the \eqn{^{207}}Pb/\eqn{^{235}}U-age and standard error, the
#' \eqn{^{206}}Pb/\eqn{^{238}}U-age and standard error, the
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb-age and standard error, the single
#' grain concordia age and standard error, and the p-value for
#' concordance, respectively.
#'
#' \item if \code{x} has class \code{UPb} and \code{type=2, 3, 4} or
#' \code{5}, returns the output of the \code{\link{concordia}}
#' function.
#'
#' \item if \code{x} has class \code{PbPb}, \code{ArAr}, \code{KCa},
#' \code{RbSr}, \code{SmNd}, \code{ReOs}, \code{LuHf}, \code{ThU} or
#' \code{UThHe} and \code{isochron=FALSE}, returns a table of Pb-Pb,
#' Ar-Ar, K-Ca, Rb-Sr, Sm-Nd, Re-Os, Lu-Hf, Th-U or U-Th-He ages and
#' their standard errors.
#'
#' \item if \code{x} has class \code{ThU} and \code{isochron=FALSE},
#' returns a 5-column table with the Th-U ages, their standard errors,
#' the initial \eqn{^{234}}U/\eqn{^{238}}U-ratios, their standard errors,
#' and the correlation coefficient between the ages and the initial
#' ratios.
#'
#' \item if \code{x} has class \code{PbPb}, \code{ArAr}, \code{KCa},
#' \code{RbSr}, \code{SmNd}, \code{ReOs}, \code{LuHf}, \code{UThHe} or
#' \code{ThU} and \code{isochron=TRUE}, returns the output of the
#' \code{\link{isochron}} function.
#'
#' \item if \code{x} has class \code{fissiontracks} and
#' \code{central=FALSE}, returns a table of fission track ages and
#' standard errors.
#'
#' \item if \code{x} has class \code{fissiontracks} or \code{UThHe}
#' and \code{central=TRUE}, returns the output of the
#' \code{\link{central}} function.
#'
#' }
#' @seealso \code{\link{concordia}}, \code{\link{isochron}},
#'     \code{\link{central}}
#' @examples
#' data(examples)
#' tUPb <- age(examples$UPb,type=1)
#' tconc <- age(examples$UPb,type=2)
#' tdisc <- age(examples$UPb,type=3)
#' tArAr <- age(examples$ArAr)
#' tiso <- age(examples$ArAr,isochron=TRUE,i2i=TRUE)
#' tcentral <- age(examples$FT1,central=TRUE)
#' @rdname age
#' @export
age.UPb <- function(x,type=1,wetherill=TRUE,exterr=TRUE,i=NA,
                    sigdig=NA,common.Pb=0,show.p=FALSE,...){
    if (common.Pb %in% c(1,2,3))
        X <- Pb0corr(x,option=common.Pb)
    else
        X <- x
    if (type==1){
        out <- UPb.age(X,exterr=exterr,i=i,sigdig=sigdig,show.p=show.p,...)
    } else if (type==2){
        out <- concordia.age(X,wetherill=wetherill,exterr=exterr)
    } else if (type %in% c(3,4,5)){
        out <- concordia.intersection.ludwig(x,wetherill=wetherill,
                                             exterr=exterr,model=type-2)
    }
    out
}
#' @rdname age
#' @export
age.PbPb <- function(x,isochron=TRUE,common.Pb=1,
                     exterr=TRUE,i=NA,sigdig=NA,...){
    if (isochron)
        out <- isochron(x,plot=FALSE,exterr=exterr,sigdig=sigdig,...)
    else
        out <- PbPb.age(x,exterr=exterr,i=i,sigdig=sigdig,common.Pb=common.Pb)
    out
}
#' @param J two-element vector with the J-factor and its standard
#'     error.
#' 
#' @param isochron logical flag indicating whether each Ar-Ar analysis
#'     should be considered separately (\code{isochron=FALSE}) or an
#'     isochron age should be calculated from all Ar-Ar analyses
#'     together (\code{isochron=TRUE}).
#' 
#' @param i2i `isochron to intercept': calculates the initial (aka
#'     `inherited', `excess', or `common')
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar, \eqn{^{40}}Ca/\eqn{^{44}}Ca,
#'     \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#'     \eqn{^{143}}Nd/\eqn{^{144}}Nd, \eqn{^{187}}Os/\eqn{^{188}}Os or
#'     \eqn{^{176}}Hf/\eqn{^{177}}Hf ratio from an isochron
#'     fit. Setting \code{i2i} to \code{FALSE} uses the default values
#'     stored in \code{settings('iratio',...)}. When applied to data
#'     of class \code{ThU}, setting \code{i2i} to \code{TRUE} applies
#'     a detrital Th-correction.
#'
#' @rdname age
#' @export
age.ArAr <- function(x,isochron=FALSE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE,exterr=exterr,sigdig=sigdig,...)
    else out <- ArAr.age(x,exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
    out
}
#' @rdname age
#' @export
age.KCa <- function(x,isochron=FALSE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE,exterr=exterr,sigdig=sigdig,...)
    else out <- KCa.age(x,exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
    out
}
#' @param central logical flag indicating whether each analysis should
#'     be considered separately (\code{central=FALSE}) or a central
#'     age should be calculated from all analyses together
#'     (\code{central=TRUE}).
#' @rdname age
#' @export
age.UThHe <- function(x,isochron=FALSE,central=FALSE,i=NA,sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE,sigdig=sigdig,...)
    else if (central) out <- central(x)
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
#' @param detritus detrital \eqn{^{230}}Th correction (only applicable
#'     when \code{x$format = 1} or \code{2}).
#'
#' \code{0}: no correction
#'
#' \code{1}: project the data along an isochron fit
#'
#' \code{2}: correct the data using an assumed initial
#' \eqn{^{230}}Th/\eqn{^{232}}Th-ratio for the detritus.
#'
#' \code{3}: correct the data using the measured present day
#' \eqn{^{230}}Th/\eqn{^{238}}U, \eqn{^{232}}Th/\eqn{^{238}}U and
#' \eqn{^{234}}U/\eqn{^{238}}U-ratios in the detritus.
#' 
#' @rdname age
#' @export
age.ThU <- function(x,isochron=FALSE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,detritus=0,...){
    if (isochron) out <- isochron(x,plot=FALSE,exterr=exterr,sigdig=sigdig,...)
    else out <- ThU.age(x,exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,detritus=detritus,...)
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
age.LuHf <- function(x,isochron=TRUE,i2i=TRUE,exterr=TRUE,i=NA,
                     sigdig=NA,...){
    age.PD(x,'Lu176',isochron=isochron,i2i=i2i,exterr=exterr,i=i,
           sigdig=sigdig,...)
}
age.PD <- function(x,nuclide,isochron=TRUE,i2i=TRUE,exterr=TRUE,i=NA,
                   sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE,sigdig=sigdig)
    else out <- PD.age(x,nuclide,exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
    out
}
# tt and st are the age and error (scalars produced by peakfit or weightedmean)
# calculated without taking into account the external errors
add.exterr <- function(x,tt,st,cutoff.76=1100,type=4){
    out <- c(tt,st)
    if (hasClass(x,'UPb')){
        if (type==1){
            R <- age_to_Pb207U235_ratio(tt,st)
            out <- get.Pb207U235.age(R[,1],R[,2],exterr=TRUE)
        } else if (type==2 | (type==4 & (tt<cutoff.76)) | (type==5)){
            R <- age_to_Pb206U238_ratio(tt,st)
            out <- get.Pb206U238.age(R[,1],R[,2],exterr=TRUE)
        } else if (type==3 | (type==4 & (tt>=cutoff.76))){
            R <- age_to_Pb207Pb206_ratio(tt,st)
            out <- get.Pb207Pb206.age(R[,1],R[,2],exterr=TRUE)
        }
    } else if (hasClass(x,'PbPb')){
        R <- age_to_Pb207Pb206_ratio(tt,st)
        out <- get.Pb207Pb206.age(R[,1],R[,2],exterr=TRUE)
    } else if (hasClass(x,'ArAr')){
        R <- get.ArAr.ratio(tt,st,x$J[1],0,exterr=FALSE)
        out <- get.ArAr.age(R[1],R[2],x$J[1],x$J[2],exterr=TRUE)
    } else if (hasClass(x,'KCa')){
        R <- get.KCa.ratio(tt,st,exterr=FALSE)
        out <- get.KCa.age(R[1],R[2],exterr=TRUE)
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

get.ages <- function(x,type=4,cutoff.76=1100,cutoff.disc=c(-15,5),
                     i2i=FALSE,common.Pb=0,detritus=0){
    if (hasClass(x,'UPb')){
        out <- filter.UPb.ages(x,type=type,cutoff.76=cutoff.76,
                               cutoff.disc=cutoff.disc,
                               exterr=FALSE,common.Pb=common.Pb)
    } else if (hasClass(x,'PbPb')){
        out <- PbPb.age(x,exterr=FALSE,common.Pb=common.Pb)
    } else if (hasClass(x,'ArAr')){
        out <- ArAr.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'KCa')){
        out <- KCa.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'UThHe')){
        out <- UThHe.age(x)
    } else if (hasClass(x,'ReOs')){
        out <- ReOs.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'SmNd')){
        out <- SmNd.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'RbSr')){
        out <- RbSr.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'LuHf')){
        out <- LuHf.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'fissiontracks')){
        out <- fissiontrack.age(x,exterr=FALSE)
    } else if (hasClass(x,'ThU')){
        out <- ThU.age(x,exterr=FALSE,i2i=i2i,detritus=detritus)
    }
    out
}
