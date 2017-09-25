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
#'
#' \item a five element vector containing \code{0/8}, \code{s[0/8]},
#' \code{4/8}, \code{s[4/8]} and \code{cov[0/8,4/8]}
#'
#' }
#'
#' OR
#'
#' \itemize{
#' \item an object of class \code{UPb}, \code{PbPb}, \code{ArAr}, \code{ThU},
#' \code{RbSr}, \code{SmNd}, \code{ReOs}, \code{LuHf}, \code{UThHe} or
#' \code{fissiontracks}.
#' }
#'
#' @param method one of either \code{'U238-Pb206'}, \code{'U235-Pb207'},
#'     \code{'Pb207-Pb206'}, \code{'Ar-Ar'}, \code{'Th-U'}, \code{'Re-Os'},
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
    if (length(x)==1) x <- c(x,0)
    if (identical(method,'U235-Pb207')){
        out <- get.Pb207U235.age(x[1],x[2],exterr)
    } else if (identical(method,'U238-Pb206')){
        out <- get.Pb206U238.age(x[1],x[2],exterr)
    } else if (identical(method,'Pb206-Pb207')){
        out <- get.Pb207Pb206.age(x[1],x[2],exterr)
    } else if (identical(method,'Ar-Ar')){
        out <- get.ArAr.age(x[1],x[2],x[3],x[4],exterr)
    } else if (identical(method,'Re-Os')){
        out <- get.ReOs.age(x[1],x[2],exterr)
    } else if (identical(method,'Rb-Sr')){
        out <- get.RbSr.age(x[1],x[2],exterr)
    } else if (identical(method,'Sm-Nd')){
        out <- get.SmNd.age(x[1],x[2],exterr)
    } else if (identical(method,'Lu-Hf')){
        out <- get.LuHf.age(x[1],x[2],exterr)
    } else if (identical(method,'Th-U')){
        out <- get.ThU.age(x[1],x[2],x[3],x[4],x[5],exterr)
    } else if (identical(method,'U-Th-He')){
        if (length(x)==6)
            out <- get.UThHe.age(x[1],x[2],x[3],x[4],x[5],x[6])
        else if (length(x)==8)
            out <- get.UThHe.age(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8])
    } else if (identical(method,'fissiontracks')){
        out <- get.EDM.age(x[1],x[2],zeta,rhoD)
    } else {
        out <- NA
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
#' @param common.Pb apply a common lead correction using one of three
#'     methods:
#'
#' \code{1}: use the isochron intercept as the initial Pb-composition
#'
#' \code{2}: use the Stacey-Kramer two-stage model to infer the initial
#' Pb-composition
#'
#' \code{3}: use the Pb-composition stored in
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')}
#'
#' @return
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
#'
#' \item{x}{ a named vector with the (weighted mean) U-Pb composition }
#' 
#' \item{cov}{ the covariance matrix of the (mean) U-Pb composition }
#' 
#' \item{age}{a 4-element vector with:
#'
#' \code{t}: the concordia age (in Ma)
#'
#' \code{s[t]}: the estimated uncertainty of \code{t}
#'
#' \code{ci[t]}: the 95\% confidence interval of \code{t} for the
#' appropriate degrees of freedom
#'
#' \code{disp[t]}: the 95\% confidence interval for \code{t} augmented
#' by \eqn{\sqrt{MSWD}} to account for overdispersed datasets.}
#' 
#' \item{mswd}{ a vector with three items (\code{equivalence},
#' \code{concordance} and \code{combined}) containing the MSWD (Mean
#' of the Squared Weighted Deviates, a.k.a the reduced Chi-squared
#' statistic outside of geochronology) of isotopic equivalence, age
#' concordance and combined goodness of fit, respectively. }
#' 
#' \item{p.value}{ a vector with three items (\code{equivalence},
#' \code{concordance} and \code{combined}) containing the p-value of
#' the Chi-square test for isotopic equivalence, age concordance and
#' combined goodness of fit, respectively. }
#' 
#' \item{df}{ the number of degrees of freedom used for the
#' \code{mswd} calculation.  These values are useful when expanding
#' the analytical uncertainties when \code{mswd>1}.}
#'
#' }
#' 
#' \item if \code{x} has class \code{UPb} and \code{type=3},
#' returns a list with the following items:
#'
#' \describe{
#'
#' \item{x}{ a two element vector with the upper and lower intercept
#' ages (if \code{wetherill=TRUE}) or the lower intercept age and
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb intercept (for Tera-Wasserburg) }
#' 
#' \item{cov}{ the covariance matrix of the elements in \code{x} }
#'
#' \item{err}{ a \code{3 x 2} matrix with the following rows:
#'
#' \code{s}: the estimated standard deviation for \code{x}
#' 
#' \code{ci}: the 95\% confidence interval of \code{x} for the
#' appropriate degrees of freedom
#'
#' \code{disp[t]}: the 95\% confidence interval for \code{x} augmented
#' by \eqn{\sqrt{MSWD}} to account for overdispersed datasets.}
#'
#' \item{df}{ the degrees of freedom of the concordia fit (concordance
#' + equivalence)}
#'
#' \item{p.value}{ p-value of a Chi-square test for age homogeneity }
#'
#' \item{mswd}{ mean square of the weighted deviates -- a
#' goodness-of-fit measure. \code{mswd > 1} indicates overdispersion
#' w.r.t the analytical uncertainties.}
#'
#' }
#'
#' \item if \code{x} has class \code{PbPb}, \code{ArAr}, \code{RbSr},
#' \code{SmNd}, \code{ReOs}, \code{LuHf} and \code{isochron=FALSE},
#' returns a table of Pb-Pb, Ar-Ar, Rb-Sr, Sm-Nd, Re-Os or Lu-Hf and
#' standard errors.
#' 
#' \item if \code{x} has class \code{PbPb}, \code{ArAr}, \code{RbSr},
#' \code{SmNd}, \code{ReOs} or \code{LuHf} and \code{isochron=TRUE},
#' returns a list with the following items:
#'
#' \describe{
#'
#' \item{a}{the intercept of the straight line fit and its standard
#' error.}
#' 
#' \item{b}{the slope of the fit and its standard error.}
#'
#' \item{cov.ab}{the covariance of the slope and intercept}
#' 
#' \item{mswd}{the mean square of the residuals (a.k.a
#'     `reduced Chi-square') statistic}
#'
#' \item{p.value}{the p-value of a Chi-square test for linearity}
#'
#' \item{y0}{a 4-element vector containing:
#'
#' \code{t}: the atmospheric \eqn{^{40}}Ar/\eqn{^{36}}Ar or initial
#' \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{187}}Os/\eqn{^{188}}Os,
#' \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd or
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf ratio.
#'
#' \code{s[t]}: the estimated standard deviation of \code{t}
#'
#' \code{ci[t]}: the 95\% confidence interval of \code{t} for the
#' appropriate degrees of freedom
#'
#' \code{disp[t]}: the 95\% confidence interval for \code{t} augmented
#' by \eqn{\sqrt{MSWD}} to account for overdispersed datasets.}
#' 
#' \item{age}{a 4-element vector containing:
#'
#' \code{y}: the \eqn{^{207}}Pb/\eqn{^{206}}Pb,
#' \eqn{^{40}}Ar/\eqn{^{39}}Ar, \eqn{^{187}}Os/\eqn{^{187}}Re,
#' \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd or
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf age and its standard error.
#'
#' \code{s[y]}: the estimated standard deviation of \code{y}
#'
#' \code{ci[y]}: the 95\% confidence interval of \code{y} for the
#' appropriate degrees of freedom
#'
#' \code{disp[y]}: the 95\% confidence interval for \code{y} augmented
#' by \eqn{\sqrt{MSWD}} to account for overdispersed datasets.}
#'
#' \item{tfact}{ the t-value for \code{df} degrees of freedom, which
#' is used for the construction of \code{ci[t]} and \code{ci[y]} }
#'
#' \item{df}{ the degrees of freedom for the isochron fit}
#'
#' \item{model}{ copied from the input parameters }
#'
#' }
#'
#' \item if \code{x} has class \code{ThU} and \code{isochron=FALSE},
#' returns a 5-column table with the Th-U ages, their standard errors,
#' the initial \eqn{^{234}U/^{238}U}-ratios, their standard errors,
#' and the correlation coefficient between the ages and the initial
#' ratios.
#'
#' \item if \code{x} has class \code{ThU} and \code{isochron=TRUE},
#' returns the output of an `Osmond Type-II' isochron, i.e.:
#'
#' \describe{
#'
#' \item{par}{the best fitting \eqn{^{234}U/^{238}U} intercept,
#' \eqn{^{234}U/^{232}Th} slope, \eqn{^{230}Th/^{238}U} intercept and
#' \eqn{^{230}Th/^{232}Th} slope.}
#'
#' \item{cov}{the covariance matrix of \code{par}.}
#'
#' \item{a}{the \eqn{^{234}U/^{238}U} intercept (i.e. the detrital
#' Th-corrected value) and its standard error.}
#' 
#' \item{b}{the \eqn{^{234}U/^{232}Th} slope and its standard error.}
#'
#' \item{cov.ab}{the covariance of \code{a} and \code{b}.}
#' 
#' \item{mswd}{the mean square of the residuals (a.k.a
#'     `reduced Chi-square') statistic.}
#'
#' \item{p.value}{the p-value of a Chi-square test for linearity.}
#'
#' \item{y0}{a 4-element vector containing:
#'
#' \code{y}: the initial \eqn{^{234}U/^{238}U}-ratio
#'
#' \code{s[y]}: the estimated standard deviation of \code{y}
#'
#' \code{ci[y]}: the 95\% confidence interval of \code{y} for the
#' appropriate degrees of freedom
#'
#' \code{disp[y]}: the 95\% confidence interval for \code{y} augmented
#' by \eqn{\sqrt{MSWD}} to account for overdispersed datasets.  }
#'
#' \item{age}{a 4-element vector containing:
#'
#' \code{t}: the Th-U isochron age
#'
#' \code{s[t]}: the estimated standard deviation of \code{t}
#'
#' \code{ci[t]}: the 95\% confidence interval of \code{t} for the
#' appropriate degrees of freedom
#'
#' \code{disp[t]}: the 95\% confidence interval for \code{t} augmented
#' by \eqn{\sqrt{MSWD}} to account for overdispersed datasets.  }
#'
#' \item{df}{the degrees of freedom for the isochron fit.}
#'
#' \item{tfact}{ the t-value for \code{df} degrees of freedom }
#'
#' }
#'
#' \item if \code{x} has class \code{UThHe} and \code{central=TRUE},
#' returns a list with the following items:
#'
#' \describe{
#'
#' \item{uvw}{ a three-element list with the weighted mean log[U/He],
#' log[Th/He] and log[Sm/He] compositions. }
#'
#' \item{covmat}{ a \code{3 x 3} covariance matrix for \code{uvw}}
#'
#' \item{mswd}{ the reduced Chi-square value for the
#' log[U/He]-log[Th/He] compositions. }
#'
#' \item{p.value}{ the p-value of concordance between the
#' log[U/He]-log[Th/He] compositions. }
#'
#' \item{age}{ a 4-element vector containing:
#'
#' \code{t}: the central age age, i.e. the U-Th-He age corresponding
#' to the composition given by \code{uvw}
#'
#' \code{s[t]}: the estimated standard deviation of \code{t}
#'
#' \code{ci[t]}: the 95\% confidence interval of \code{t} for the
#' appropriate degrees of freedom
#'
#' \code{disp[t]}: the 95\% confidence interval for \code{t} augmented
#' by \eqn{\sqrt{MSWD}} to account for overdispersed datasets.  }
#'
#'
#' \item{df}{the degrees of freedom for the data fit.}
#'
#' \item{tfact}{ the t-value for \code{df} degrees of freedom }
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
#' \item{age}{ a three-element vector with the central age, its
#' standard error and a 95\% confidence interval for the appropriate
#' degrees of freedom. }
#'
#' \item{disp}{ the (over)dispersion of the single grain ages beyond
#' the formal analytical uncertainties. }
#'
#' \item{df}{ degrees of freedom for the Chi-square test}
#'
#' }
#'
#' }
#' 
#' @examples
#' data(examples)
#' print(age(examples$UPb))
#' print(age(examples$UPb,type=1))
#' print(age(examples$UPb,type=2))
#' @rdname age
#' @export
age.UPb <- function(x,type=1,wetherill=TRUE,exterr=TRUE,i=NA,
                    sigdig=NA,common.Pb=0,...){
    if (common.Pb %in% c(1,2,3))
        X <- common.Pb.correction(x,option=common.Pb)
    else
        X <- x
    if (type==1)
        out <- UPb.age(X,exterr=exterr,i=i,sigdig=sigdig,...)
    else if (type==2)
        out <- concordia.age(X,wetherill=wetherill,exterr=exterr,...)
    else if (type==3)
        out <- concordia.intersection.ludwig(x,wetherill=wetherill,exterr=exterr)
    out
}
#' @rdname age
#' @export
age.PbPb <- function(x,isochron=TRUE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE,exterr=exterr,sigdig=sigdig,...)
    else out <- PbPb.age(x,exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
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
#'     \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#'     \eqn{^{143}}Nd/\eqn{^{144}}Nd, \eqn{^{187}}Os/\eqn{^{188}}Os or
#'     \eqn{^{176}}Hf/\eqn{^{177}}Hf ratio from an isochron
#'     fit. Setting \code{i2i} to \code{FALSE} uses the default values
#'     stored in \code{settings('iratio',...)}  or zero (for the Pb-Pb
#'     method). When applied to data of class \code{ThU}, setting
#'     \code{i2i} to \code{TRUE} applies a detrital Th-correction.
#' @rdname age
#' @export
age.ArAr <- function(x,isochron=FALSE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE,exterr=exterr,sigdig=sigdig,...)
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
age.ThU <- function(x,isochron=FALSE,i2i=TRUE,exterr=TRUE,i=NA,sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE,exterr=exterr,sigdig=sigdig,...)
    else out <- ThU.age(x,exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
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
    out <- c(0,0)
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
