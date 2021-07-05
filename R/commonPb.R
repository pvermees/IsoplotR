#' @title Common Pb correction
#'
#' @description
#' Applies a common-Pb correction to a U-Pb dataset using either the
#' Stacey-Kramers mantle evolution model, isochron regression, or any
#' nominal inital Pb isotope composition.
#'
#' @details
#'
#' \code{IsoplotR} implements nine different methods to correct for
#' the presence of non-radiogenic (`common') lead. This includes three
#' strategies tailored to datasets that include \eqn{^{204}}Pb
#' measurements, three strategies tailored to datasets that include
#' \eqn{^{208}}Pb measurements, and a further three strategies for
#' datasets that only include \eqn{^{206}}Pb and
#' \eqn{^{207}}Pb.
#'
#' \eqn{^{204}}Pb is the only one of lead's four stable isotopes that
#' does not have a naturally occurring radioactive parent. This makes
#' it very useful for common-Pb correction:
#' 
#' \eqn{\left[\frac{{}^{206|7}Pb}{{}^{204}Pb}\right]_r =
#' \left[\frac{{}^{206|7}Pb}{{}^{204}Pb}\right]_m -
#' \left[\frac{{}^{206|7}Pb}{{}^{204}Pb}\right]_\circ}
#'
#' where \eqn{[{}^{206|7}Pb/^{204}Pb]_r} marks the radiogenic
#' \eqn{{}^{206}}Pb or \eqn{{}^{207}}Pb component;
#' \eqn{[{}^{206|7}Pb/^{204}Pb]_m} is the measured ratio; and
#' \eqn{[{}^{206|7}Pb/^{204}Pb]_\circ} is the non-radiogenic component.
#' 
#' \code{IsoplotR} offers three different ways to determine
#' \eqn{[{}^{206|7}Pb/^{204}Pb]_\circ}. The first and easiest option
#' is to simply use a nominal value such as the
#' \eqn{{}^{206|7}}Pb/\eqn{^{204}}Pb-ratio of a cogenetic feldspar,
#' assuming that this is representative for the common-Pb composition
#' of the entire sample. A second method is to determine the
#' non-radiogenic isotope composition by fitting an isochron line
#' through multiple aliquots of the same sample, using the
#' 3-dimensional regression algorithm of Ludwig (1998).
#'
#' Unfortunately, neither of these two methods is applicable to
#' detrital samples, which generally lack identifiable cogenetic
#' minerals and aliquots. For such samples, \code{IsoplotR} infers the
#' common-Pb composition from the two-stage crustal evolution model of
#' Stacey and Kramers (1975). The second stage of this model is
#' described by:
#'
#' \eqn{\left[\frac{{}^{206}Pb}{{}^{204}Pb}\right]_\circ =
#' \left[\frac{{}^{206}Pb}{{}^{204}Pb}\right]_{3.7Ga} +
#' \left[\frac{{}^{238}U}{{}^{204}Pb}\right]_{sk}
#' \left(e^{\lambda_{238}3.7Ga}-e^{\lambda_{238}t}\right)}
#'
#' where \eqn{\left[{}^{206}Pb/{}^{204}Pb\right]_{3.7Ga} = 11.152} and
#' \eqn{\left[{}^{238}U/{}^{204}Pb\right]_{sk} = 9.74}. These
#' Equations can be solved for \eqn{t} and
#' \eqn{\left[{}^{206}Pb/{}^{204}Pb\right]_\circ} using the method of
#' maximum likelihood. The \eqn{{}^{207}}Pb/\eqn{{}^{204}}Pb-ratio is
#' corrected in exactly the same way, using
#' \eqn{\left[{}^{207}Pb/{}^{204}Pb\right]_{3.7Ga} = 12.998}.
#'
#' In the absence of \eqn{^{204}}Pb measurements, a \eqn{^{208}}Pb-based
#' common lead correction can be used:
#'
#' \eqn{\frac{{}^{206|7}Pb_r}{{}^{208}Pb_\circ} =
#' \frac{{}^{206|7}Pb_m}{{}^{208}Pb_\circ} -
#' \left[\frac{{}^{206|7}Pb}{{}^{208}Pb}\right]_\circ}
#'
#' where \eqn{{}^{208}Pb_\circ} marks the non-radiogenic
#' \eqn{{}^{208}Pb}-component, which is obtained by removing the
#' radiogenic component for any given age.
#'
#' If neither \eqn{{}^{204}}Pb nor \eqn{{}^{208}}Pb were measured,
#' then a \eqn{^{207}} Pb-based common lead correction can be used:
#'
#' \eqn{ \left[\frac{{}^{207}Pb}{{}^{206}Pb}\right]_m = f
#' \left[\frac{{}^{207}Pb}{{}^{206}Pb}\right]_\circ + (1-f)
#' \left[\frac{{}^{207}Pb}{{}^{204}Pb}\right]_r}
#'
#' where \eqn{f} is the fraction of common lead, and
#' \eqn{[{}^{207}Pb/{}^{206}Pb]_r} is obtained by projecting the U-Pb
#' measurements on the concordia line in Tera-Wasserburg space.  Like
#' before, the initial lead composition
#' \eqn{[{}^{207}Pb/{}^{206}Pb]_\circ} can be obtained in three
#' possible ways: by analysing a cogenetic mineral, by isochron
#' regression through multiple aliquots, or from the Stacey and
#' Kramers (1975) model.
#'
#' Besides the common-Pb problem, a second reason for U-Pb discordance
#' is radiogenic Pb-loss during igneous and metamorphic activity.
#' This moves the data away from the concordia line along a linear
#' array, forming an isochron or `discordia' line.  \code{IsoplotR}
#' fits this line using the Ludwig (1998) algorithm. If the data are
#' plotted on a Wetherill concordia diagram, the program will not only
#' report the usual lower intercept with the concordia line, but the
#' upper intercept as well. Both values are geologically meaningful as
#' they constrain both the initial igneous age as well as the timing
#' of the partial resetting event.
#'
#' 
#' @param x an object of class \code{UPb}
#'
#' @param option one of either
#' 
#' \code{1}: nominal common Pb isotope composition
#'
#' \code{2}: isochron regression
#'
#' \code{3}: Stacey-Kramers correction
#' 
#' @param omit4c vector with indices of aliquots that should be
#'     omitted from the isochron regression (only used if
#'     \code{option=2})
#'
#' @return
#' Returns a list in which \code{x.raw} contains the original data and
#' \code{x} the common Pb-corrected compositions. All other items in
#' the list are inherited from the input data.
#'
#' @references
#' Ludwig, K.R., 1998. On the treatment of concordant uranium-lead
#' ages. Geochimica et Cosmochimica Acta, 62(4), pp.665-676.
#' 
#' Stacey, J.T. and Kramers, 1., 1975. Approximation of terrestrial
#' lead isotope evolution by a two-stage model. Earth and Planetary
#' Science Letters, 26(2), pp.207-221.
#' 
#' @examples
#' attach(examples)
#' corrected <- Pb0corr(UPb,option=2)
#' concordia(corrected)
#' # produces identical results as:
#' dev.new()
#' concordia(UPb,common.Pb=2)
#' @export
Pb0corr <- function(x,option=3,omit4c=NULL){
    ns <- length(x)
    out <- x
    out$x.raw <- x$x
    if (option == 1){
        x.corr <- common.Pb.nominal(x)
    } else if (option == 2){
        x.corr <- common.Pb.isochron(x,omit=omit4c)
    } else if (option == 3){
        x.corr <- common.Pb.stacey.kramers(x)
    } else {
        return(out)
    }
    if (x$format==1){
        out$x <- tw2w(x.corr,format=2)
    } else if (x$format==2){
        out$x <- x.corr
    } else if (x$format==3){
        out$x[,1:4] <- tw2w(x.corr,format=2)[,1:4,drop=FALSE]
        out$x[,c('Pb207Pb206','errPb207Pb206')] <- x.corr[,3:4,drop=FALSE]
    } else if (x$format==4){
        out$x[,c('Pb207U235','errPb207U235',
                 'Pb206U238','errPb206U238')] <- x.corr[,1:4,drop=FALSE]
        out$x[,'rhoXY'] <- x.corr[,5]
        out$x[,c('Pb204U238','errPb204U238','rhoXZ','rhoYZ')] <- 0
    } else if (x$format==5){
        tw <- w2tw(x.corr,format=1)
        out$x[,c('U238Pb206','errU238Pb206',
                 'Pb207Pb206','errPb207Pb206')] <- tw[,1:4,drop=FALSE]
        out$x[,'rhoXY'] <- tw[,5]
        out$x[,c('Pb204Pb206','errPb204Pb206','rhoXZ','rhoYZ')] <- 0
    } else if (x$format==6){
        tw <- w2tw(x.corr,format=1)
        out$x[,c('Pb207U235','errPb207U235',
                 'Pb206U238','errPb206U238')] <- x.corr[,1:4,drop=FALSE]
        out$x[,c('Pb204U238','errPb204U238')] <- 0
        out$x[,c('Pb207Pb206','errPb207Pb206')] <- tw[,3:4,drop=FALSE]
        out$x[,c('Pb204Pb207','errPb204Pb207',
                 'Pb204Pb206','errPb204Pb206')] <- 0
    } else if (x$format==7){
        out$x <- x.corr
    }  else if (x$format==8){
        out$x <- w2tw(x.corr,format=7)
    } else {
        stop('Incorrect input format.')
    }
    out
}

correct.common.Pb.without.204 <- function(x,i,c76,tt=NULL){
    tw <- tera.wasserburg(x,i)
    m86 <- tw$x['U238Pb206']
    m76 <- tw$x['Pb207Pb206']
    if (is.null(tt)){
        tt <- project.concordia(m86,m76,c76,d=x$d[i])
        cctw <- age_to_terawasserburg_ratios(tt=tt,st=0,d=x$d[i])
        r86 <- cctw$x['U238Pb206']
        r76 <- cctw$x['Pb207Pb206']
        cnames <- c('U238Pb206','Pb207Pb206')
        E <- tw$cov[cnames,cnames]
        sr86 <- sqrt(E[1,1])
        sr76 <- sqrt(E[2,2])
        rho <- stats::cov2cor(E)[1,2]
        out <- c(r86,sr86,r76,sr76,rho)
        names(out) <- c('U238Pb206','errU238Pb206',
                        'Pb207Pb206','errPb207Pb206','rhoXY')
    } else {
        cctw <- age_to_terawasserburg_ratios(tt=tt,st=0,d=x$d[i])
        r86 <- cctw$x['U238Pb206']
        r76 <- cctw$x['Pb207Pb206']
        slope <- (c76-r76)/r86
        p76 <- m76 + slope*m86
        out <- correct.common.Pb.without.204(x=x,i=i,c76=p76)
    }
    out
}
correct.common.Pb.with.204 <- function(x,i,c46,c47,tt=NULL,cc=FALSE){
    ir <- get.UPb.isochron.ratios.204(x,i=i) # 3806, 0406, 3507, 0407
    Jp <- matrix(0,2,4)
    if (is.null(tt)){ # line through measurement
        p3507 <- ir$x['U235Pb207']*c47/(c47-ir$x['Pb204Pb207'])
        p3806 <- ir$x['U238Pb206']*c46/(c46-ir$x['Pb204Pb206'])
        Jp[1,3] <- c47/(c47-ir$x['Pb204Pb207'])
        Jp[1,4] <- ir$x['U235Pb207']*c47/(c47-ir$x['Pb204Pb207'])^2
        Jp[2,1] <- c46/(c46-ir$x['Pb204Pb206'])
        Jp[2,2] <- ir$x['U238Pb206']*c46/(c46-ir$x['Pb204Pb206'])^2
    } else { # line parallel to isochron
        r3507 <- age_to_U235Pb207_ratio(tt,d=x$d[i])[1]
        p3507 <- ir$x['U235Pb207'] + ir$x['Pb204Pb207']*r3507/c47
        r3806 <- age_to_U238Pb206_ratio(tt,d=x$d[i])[1]
        p3806 <- ir$x['U238Pb206'] + ir$x['Pb204Pb206']*r3806/c46
        Jp[1,3] <- 1
        Jp[1,4] <- r3507/c47
        Jp[2,1] <- 1
        Jp[2,2] <- r3806/c46
    }
    Ep <- Jp %*% ir$cov %*% t(Jp)
    J <- matrix(0,2,2)
    J[1,1] <- -1/p3507^2
    J[2,2] <- -1/p3806^2
    E <- J %*% Ep %*% t(J)
    if (cc){
        out <- list()
        cnames <- c('Pb207U235','Pb206U238')
        out$x <- 1/c(p3507,p3806)
        out$cov <- E
        names(out$x) <- cnames
        rownames(out$cov) <- cnames
        colnames(out$cov) <- cnames
    } else {
        out <- rep(NA,5)
        names(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238','rhoXY')
        out[1] <- 1/p3507
        out[3] <- 1/p3806
        out[c(2,4)] <- sqrt(diag(E))
        out[5] <- stats::cov2cor(E)[1,2]
    }
    out
}
correct.common.Pb.with.208 <- function(x,i,tt,c0806,c0807,cc=FALSE){
    ir <- get.UPb.isochron.ratios.208(x,i,tt=tt) # 3806 08c06 3507 08c07 3238 3208 06c08 07c08
    r3507 <- age_to_U235Pb207_ratio(tt,d=x$d[i])[1]
    p3507 <- ir$x['U235Pb207'] + ir$x['Pb208cPb207']*r3507/c0807
    r3806 <- age_to_U238Pb206_ratio(tt,d=x$d[i])[1]
    p3806 <- ir$x['U238Pb206'] + ir$x['Pb208cPb206']*r3806/c0806
    r3208 <- 1/age_to_Pb208Th232_ratio(tt)[1]
    p3208 <- ir$x['Th232Pb208'] + ir$x['Pb207cPb208']*r3208*c0807
    p <- c(p3507,p3806,p3208,ir$x['Th232U238']) # projected compositions
    Jp <- matrix(0,4,8)
    Jp[1,3] <- 1
    Jp[1,4] <- r3507/c0807
    Jp[2,1] <- 1
    Jp[2,2] <- r3806/c0806
    Jp[3,6] <- 1
    Jp[3,8] <- r3208*c0807
    Jp[4,5] <- 1
    Ep <- Jp %*% ir$cov %*% t(Jp)
    J <- matrix(0,4,4)
    J[1,1] <- -1/p[1]^2
    J[2,2] <- -1/p[2]^2
    J[3,3] <- -1/p[3]^2
    J[4,4] <- 1
    E <- J %*% Ep %*% t(J)
    if (cc){
        out <- list()
        cnames <- c('Pb207U235','Pb206U238')
        out$x <- 1/p[1:2]
        out$cov <- E[1:2,1:2]
        names(out$x) <- cnames
        rownames(out$cov) <- cnames
        colnames(out$cov) <- cnames
    } else {
        out <- rep(NA,14)
        names(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238',
                        'Pb208Th232','errPb208Th232','Th232U238','errTh232U238',
                        'rhoXY','rhoXZ','rhoXW','rhoYZ','rhoYW','rhoZW')
        out[c(1,3,5,7)] <- c(1/p[1:3],p[4])
        cormat <- matrix(0,4,4)
        pos <- which(diag(E)>0)
        cormat[pos,pos] <- stats::cov2cor(E[pos,pos])
        out[c(2,4,6,8)] <- sqrt(diag(E))
        out[9:11] <- cormat[1,2:4]
        out[12:13] <- cormat[2,3:4]
        out[14] <- cormat[3,4]
    }
    out
}

common.Pb.stacey.kramers <- function(x){
    ns <- length(x)
    if (x$format %in% c(1,2,3)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('U238Pb206','errU238Pb206',
                           'Pb207Pb206','errPb207Pb206','rhoXY')
        for (i in 1:ns){
            maxt <- get.Pb207Pb206.age(x,i=i)[1]
            tint <- stats::optimise(SS.SK.without.204,
                                    interval=c(0,maxt),x=x,i=i)$minimum
            i6474 <- stacey.kramers(tint)
            c76 <- i6474[,'i74']/i6474[,'i64']
            out[i,] <- correct.common.Pb.without.204(x=x,i=i,c76=c76,tt=tint)
        }
    } else if (x$format %in% c(4,5,6)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb206U238','errPb206U238','rhoXY')
        for (i in 1:ns){
            maxt <- get.Pb207Pb206.age(x,i=i)[1]
            tint <- stats::optimise(SS.SK.with.204,interval=c(0,maxt),x=x,i=i)$minimum
            c6474 <- stacey.kramers(tint)
            c46 <- 1/c6474[,'i64']
            c47 <- 1/c6474[,'i74']
            out[i,] <- correct.common.Pb.with.204(x,i=i,tt=tint,c46=c46,c47=c47)
        }
    } else if (x$format%in%c(7,8)){
        out <- matrix(0,ns,14)
        colnames(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','Th232U238','errTh232U238',
                           'rhoXY','rhoXZ','rhoXW','rhoYZ','rhoYW','rhoZW')
        for (i in 1:ns){
            maxt <- get.Pb208Th232.age(x,i=i)[1]
            tint <- stats::optimise(SS.SK.with.208,interval=c(0,maxt),x=x,i=i)$minimum
            c678 <- stacey.kramers(tint)
            c86 <- c678[,'i84']/c678[,'i64']
            c87 <- c678[,'i84']/c678[,'i74']
            out[i,] <- correct.common.Pb.with.208(x,i=i,tt=tint,c0806=c86,c0807=c87)
        }
    } else {
        stop('Invalid input format.')
    }
    out
}

common.Pb.isochron <- function(x,omit=NULL){
    ns <- length(x)
    calcit <- (1:ns)%ni%omit
    fit <- ludwig(subset(x,subset=calcit))
    tt <- fit$par[1]
    if (x$format %in% c(1,2,3)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('U238Pb206','errU238Pb206','Pb207Pb206',
                           'errPb207Pb206','rhoXY')
        c76 <- fit$par['76i']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.without.204(x,i=i,c76=c76,tt=tt)
        }
    } else if (x$format %in% c(4,5,6)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb206U238','errPb206U238','rhoXY')
        c46 <- 1/fit$par['64i']
        c47 <- 1/fit$par['74i']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.204(x,i=i,c46=c46,c47=c47,tt=tt)
        }
    } else if (x$format%in%c(7,8)){
        out <- matrix(0,ns,14)
        colnames(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','Th232U238','errTh232U238',
                           'rhoXY','rhoXZ','rhoXW','rhoYZ','rhoYW','rhoZW')
        c0806 <- 1/fit$par['68i']
        c0807 <- 1/fit$par['78i']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.208(x,i,tt=tt,c0806=c0806,c0807=c0807)
        }
    }
    out
}

common.Pb.nominal <- function(x){
    ns <- length(x)
    if (x$format %in% c(1,2,3)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('U238Pb206','errU238Pb206',
                           'Pb207Pb206','errPb207Pb206','rhoXY')
        c76 <- settings('iratio','Pb207Pb206')[1]
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.without.204(x,i=i,c76=c76)
        }
    } else if (x$format %in% c(4,5,6)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb206U238','errPb206U238','rhoXY')
        c46 <- 1/settings('iratio','Pb206Pb204')[1]
        c47 <- 1/settings('iratio','Pb207Pb204')[1]
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.204(x,i=i,c46=c46,c47=c47)
        }
    } else if (x$format%in%c(7,8)){
        out <- matrix(0,ns,14)
        colnames(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','Th232U238','errTh232U238',
                           'rhoXY','rhoXZ','rhoXW','rhoYZ','rhoYW','rhoZW')
        c0806 <- settings('iratio','Pb208Pb206')[1]
        c0807 <- settings('iratio','Pb208Pb207')[1]
        for (i in 1:ns){
            tint <- stats::optimise(SS.with.208,interval=c(0,5000),
                                    c0806=c0806,c0807=c0807,x=x,i=i)$minimum
            out[i,] <- correct.common.Pb.with.208(x,i,tt=tint,c0806=c0806,c0807=c0807)
        }
    }
    out
}

SS.SK.without.204 <- function(tt,x,i){
    tw <- tera.wasserburg(x,i)
    X <- tw$x['U238Pb206']
    Y <- tw$x['Pb207Pb206']
    i6474 <- stacey.kramers(tt)
    cct <- age_to_terawasserburg_ratios(tt,st=0,d=x$d[i])
    a <- i6474[2]/i6474[1] # intercept
    b <- (cct$x['Pb207Pb206']-a)/cct$x['U238Pb206'] # slope
    cnames <- c('U238Pb206','Pb207Pb206')
    covmat <- tw$cov[cnames,cnames]
    omega <- solve(covmat)
    x.fitted <- (X*omega[1,1]+Y*omega[1,2]-a*omega[1,2])/
                (omega[1,1]+b*omega[1,2])
    y.fitted <- a + b*x.fitted
    d <- cbind(X-x.fitted,Y-y.fitted)
    as.numeric(d %*% omega %*% t(d))
}
SS.SK.with.204 <- function(tt,x,i){
    c6474 <- stacey.kramers(tt)
    c46 <- 1/c6474[1,'i64']
    c47 <- 1/c6474[1,'i74']    
    ccw <- correct.common.Pb.with.204(x=x,i=i,c46=c46,c47=c47,tt=tt,cc=TRUE)
    LL.concordia.age(tt,ccw,mswd=TRUE,exterr=FALSE,d=x$d[i])
}
SS.SK.with.208 <- function(tt,x,i){
    i678 <- stacey.kramers(tt)
    c0806 <- i678[,'i84']/i678[,'i64']
    c0807 <- i678[,'i84']/i678[,'i74']
    SS.with.208(tt,x,i,c0806,c0807)
}
SS.with.208 <- function(tt,x,i,c0806,c0807){
    ccw <- correct.common.Pb.with.208(x,i=i,tt=tt,c0806=c0806,c0807=c0807,cc=TRUE)
    LL.concordia.age(tt,ccw,mswd=TRUE,exterr=FALSE,d=x$d[i])
}

stacey.kramers <- function(tt,inverse=FALSE){
    nt <- length(tt)
    sk.206.204 <- rep(0,nt)
    sk.207.204 <- rep(0,nt)
    sk.208.204 <- rep(0,nt)
    sk.238.204 <- rep(0,nt)
    sk.232.204 <- rep(0,nt)
    ti <- rep(0,nt)
    young <- which(tt < 3700)
    old <- which(tt >= 3700)
    sk.206.204[young] <- 11.152
    sk.207.204[young] <- 12.998
    sk.208.204[young] <- 31.23
    sk.238.204[young] <- 9.74
    sk.232.204[young] <- 36.84
    ti[young] <- 3700
    sk.206.204[old] <- 9.307
    sk.207.204[old] <- 10.294
    sk.208.204[old] <- 29.487
    sk.238.204[old] <- 7.19
    sk.232.204[old] <- 33.21
    ti[old] <- 4570
    U238U235 <- settings('iratio','U238U235')[1]
    l2 <- lambda('Th232')[1]
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    i64 <- sk.206.204 + sk.238.204*(exp(l8*ti)-exp(l8*tt))
    i74 <- sk.207.204 + sk.238.204*(exp(l5*ti)-exp(l5*tt))/U238U235
    i84 <- sk.208.204 + sk.232.204*(exp(l2*ti)-exp(l2*tt))
    if (inverse){ # for Pb-Pb data
        out <- cbind(1/i64,i74/i64,i84/64)
        colnames(out) <- c('i46','i76','i86')
    } else {
        out <- cbind(i64,i74,i84)
        colnames(out) <- c('i64','i74','i84')
    }
    out
}
# TODO: add option for Pb207Pb206 ratios to be used with inverse isochrons
sk2t <- function(Pb206Pb204=rep(NA,2),Pb207Pb204=rep(NA,2)){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    ti.young <- 3700
    sk.206.204.young <- 11.152
    sk.207.204.young <- 12.998
    sk.208.204.young <- 31.23
    sk.238.204.young <- 9.74
    sk.232.204.young <- 36.84
    ti.old <- 4570
    sk.206.204.old <- 9.307
    sk.207.204.old <- 10.294
    sk.208.204.old <- 29.487
    sk.238.204.old <- 7.19
    sk.232.204.old <- 33.21
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    out <- c(0,ti.old)
    # 1. 206/204
    min64 <- sk.206.204.old
    max64 <- sk.206.204.young+sk.238.204.young*(exp(l8*ti.young)-1)
    good64 <- !is.na(Pb206Pb204)
    big64 <- good64 & (Pb206Pb204>max64)
    small64 <- good64 & (Pb206Pb204<min64)
    mid64 <- good64 & !big64 & !small64
    out[big64] <- 0
    out[small64] <- ti.old
    out[mid64] <- log( exp(l8*ti.young) +
                       (sk.206.204.young-Pb206Pb204[mid64])/sk.238.204.young )/l8
    if (any(mid64) && out[mid64]>ti.young)
        out[mid64] <- log( exp(l8*ti.old) +
                           (sk.206.204.old-Pb206Pb204[mid64])/sk.238.204.old )/l8
    # 2. 207/204
    min74 <- sk.207.204.old
    max74 <- sk.207.204.young + sk.238.204.young*(exp(l5*ti.young)-1)/U
    good74 <- !is.na(Pb207Pb204)
    big74 <- good74 & (Pb207Pb204>max74)
    small74 <- good74 & (Pb207Pb204<min74)
    mid74 <- good74 & !big74 & !small74
    out[big74] <- 0
    out[small74] <- ti.old
    out[mid74] <- log( exp(l5*ti.young) +
                       U*(sk.207.204.young-Pb207Pb204[mid74])/sk.238.204.young )/l5
    if (any(mid74) && out[mid74]>ti.young)
        out[mid74] <- log( exp(l5*ti.old) +
                           (sk.207.204.old-Pb207Pb204[mid74])/sk.238.204.old )/l5
    out
}

project.concordia <- function(m86,m76,c76,d=diseq()){
    t68 <- get.Pb206U238.age(1/m86,d=d)[1]
    t76 <- get.Pb207Pb206.age(m76,d=d)[1]
    a <- c76
    b <- (m76-c76)/m86
    neg <- (c76>m76) # negative slope?
    pos <- !neg
    above <- (t76>t68) # above concordia?
    below <- !above
    search.range <- c(t68,t76)
    go.ahead <- FALSE
    if (pos & above){
        go.ahead <- TRUE
    } else if (pos & below){
        go.ahead <- TRUE
    } else if (neg & above){
        search.range <- c(0,t68)
        go.ahead <- TRUE
    } else if (neg & below){
        if (neg){
            search.range[1] <- 0
            search.range[2] <- get.search.limit(a=a,b=b,d=d,m=0,M=5000)
            if (!is.na(search.range[2])) go.ahead <- TRUE
        } else {
            tm <- get.search.limit(a=a,b=b,d=d,m=0,M=5000)
            tM <- get.search.limit(a=a,b=b,d=d,m=5000,M=0)
            if (!is.na(tm) | !is.na(tM))
                go.ahead <- TRUE
            if (is.na(tm) & !is.na(tM))
                search.range <- c(tM,5000)
            else if (is.na(tM) & !is.na(tm))
                search.range <- c(0,tm)
            else if (!is.na(tm) & !is.na(tM) & (t68<tm))
                search.range <- c(0,tm)
            else if (!is.na(tm) & !is.na(tM) & (t68>tM))
                search.range <- c(tM,5000)
        }                
    }
    if (go.ahead){
        out <- stats::uniroot(intersection.misfit.york,
                              search.range,a=a,b=b,d=d)$root
    } else {
        out <- t68
    }
    out
}
get.search.limit <- function(a,b,d,m,M){
    ttt <- seq(from=m,to=M,length.out=100)
    for (tt in ttt){
        misfit <- intersection.misfit.york(tt,a=a,b=b,d=d)
        if (misfit<0) return(tt)
    }
    return(NA)
}
