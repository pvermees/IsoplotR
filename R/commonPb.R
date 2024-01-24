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
        out$x[,'rXY'] <- x.corr[,5]
        out$x[,c('Pb204U238','errPb204U238','rXZ','rYZ')] <- 0
    } else if (x$format==5){
        tw <- w2tw(x.corr,format=1)
        out$x[,c('U238Pb206','errU238Pb206',
                 'Pb207Pb206','errPb207Pb206')] <- tw[,1:4,drop=FALSE]
        out$x[,'rXY'] <- tw[,5]
        out$x[,c('Pb204Pb206','errPb204Pb206','rXZ','rYZ')] <- 0
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
    } else if (x$format==9){
        out$x[,'U238Pb206'] <- 1/x.corr[,'Pb206U238']
        out$x[,'errU238Pb206'] <- x.corr[,'errPb206U238']*out$x[,'U238Pb206']^2
        out$x[,c('Pb204Pb206','errPb204Pb206','rXY')] <- 0
    } else if (x$format==10){
        out$x[,'U235Pb207'] <- 1/x.corr[,'Pb207U235']
        out$x[,'errU235Pb207'] <- x.corr[,'errPb207U235']*out$x[,'U235Pb207']^2
        out$x[,c('Pb204Pb207','errPb204Pb207','rXY')] <- 0
    } else if (x$format==11){
        out$x[,'U238Pb206'] <- 1/x.corr[,'Pb206U238']
        out$x[,'errU238Pb206'] <- x.corr[,'errPb206U238']*out$x[,'U238Pb206']^2
        out$x[,c('Pb208Pb206','errPb208Pb206','rXY')] <- 0
    } else if (x$format==12){
        out$x[,'U235Pb207'] <- 1/x.corr[,'Pb207U235']
        out$x[,'errU235Pb207'] <- x.corr[,'errPb207U235']*out$x[,'U235Pb207']^2
        out$x[,c('Pb208Pb207','errPb208Pb207','rXY')] <- 0
    } else if (x$format==85){
        tw <- w2tw(x.corr,format=1)
        out$x[,c('U238Pb206','errU238Pb206',
                 'Pb207Pb206','errPb207Pb206')] <- tw[,1:4,drop=FALSE]
        out$x[,'rXY'] <- tw[,5]
        out$x[,c('Pb208Pb206','errPb208Pb206','rXZ','rYZ')] <- 0
    } else if (x$format==119){
        out$x[,'U238Pb206'] <- 1/x.corr[,'Pb206U238']
        out$x[,'errU238Pb206'] <- x.corr[,'errPb206U238']*out$x[,'U238Pb206']^2
        out$x[,c('Pb208Pb206','errPb208Pb206','rXY')] <- 0
    } else if (x$format==1210){
        out$x[,'U235Pb207'] <- 1/x.corr[,'Pb207U235']
        out$x[,'errU235Pb207'] <- x.corr[,'errPb207U235']*out$x[,'U235Pb207']^2
        out$x[,c('Pb208Pb207','errPb208Pb207','rXY')] <- 0
    } else {
        stop('Incorrect input format.')
    }
    out
}

correct.common.Pb.without.20x <- function(x,i,c76,tt=NULL){
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
        rXY <- stats::cov2cor(E)[1,2]
        out <- c(r86,sr86,r76,sr76,rXY)
        names(out) <- c('U238Pb206','errU238Pb206',
                        'Pb207Pb206','errPb207Pb206','rXY')
    } else {
        cctw <- age_to_terawasserburg_ratios(tt=tt,st=0,d=x$d[i])
        r86 <- cctw$x['U238Pb206']
        r76 <- cctw$x['Pb207Pb206']
        slope <- (c76-r76)/r86
        p76 <- m76 + slope*m86
        out <- correct.common.Pb.without.20x(x=x,i=i,c76=p76)
    }
    out
}
correct.common.Pb.with.20x <- function(x,i,cx6=NULL,cx7=NULL,tt=NULL,cc=FALSE){
    ir <- get.UPb.isochron.ratios.20x(x,i=i) # (3806, 0406), (3507, 0407)
    ni <- ifelse(x$format%in%c(4,5,6),2,1)
    Jp <- matrix(0,ni,2*ni)
    Pbx6label <- ifelse(x$format%in%c(85,119),'Pb208Pb206','Pb204Pb206')
    Pbx7label <- ifelse(x$format%in%c(85,1210),'Pb208Pb207','Pb204Pb207')
    if (is.null(tt)){ # line through measurement
        if (x$format%in%c(4,5,6,10,85,1210)){
            p3507 <- ir$x['U235Pb207']*cx7/(cx7-ir$x[Pbx7label])
            Jp[1,2*ni-1] <- cx7/(cx7-ir$x[Pbx7label])
            Jp[1,2*ni] <- ir$x['U235Pb207']*cx7/(cx7-ir$x[Pbx7label])^2
        }
        if (x$format%in%c(4,5,6,9,85,119)){
            p3806 <- ir$x['U238Pb206']*cx6/(cx6-ir$x[Pbx6label])
            Jp[ni,1] <- cx6/(cx6-ir$x[Pbx6label])
            Jp[ni,2] <- ir$x['U238Pb206']*cx6/(cx6-ir$x[Pbx6label])^2
        }
    } else { # line parallel to isochron
        if (x$format%in%c(4,5,6,10,85,1210)){
            r3507 <- age_to_U235Pb207_ratio(tt,d=x$d[i])[1]
            p3507 <- ir$x['U235Pb207'] + ir$x[Pbx7label]*r3507/cx7
            Jp[1,2*ni-1] <- 1
            Jp[1,2*ni] <- r3507/cx7
        }
        if (x$format%in%c(4,5,6,9,85,119)){
            r3806 <- age_to_U238Pb206_ratio(tt,d=x$d[i])[1]
            p3806 <- ir$x['U238Pb206'] + ir$x[Pbx6label]*r3806/cx6
            Jp[ni,1] <- 1
            Jp[ni,2] <- r3806/cx6
        }
    }
    J <- matrix(0,ni,ni)
    Ep <- Jp %*% ir$cov %*% t(Jp)
    if (x$format%in%c(4,5,6,10,85,1210)) J[1,1] <- -1/p3507^2
    if (x$format%in%c(4,5,6,9,85,119)) J[ni,ni] <- -1/p3806^2
    E <- J %*% Ep %*% t(J)
    if (cc){
        out <- list()
        cnames <- c('Pb207U235','Pb206U238')
        out$x <- 1/c(p3507,p3806)
        out$cov <- E
        names(out$x) <- cnames
        rownames(out$cov) <- cnames
        colnames(out$cov) <- cnames
    } else if (x$format%in%c(9,119)){
        out <- c('Pb206U238'=unname(1/p3806),'errPb206U238'=unname(sqrt(E)))
    } else if (x$format%in%c(10,1210)){
        out <- c('Pb207U235'=unname(1/p3507),'errPb207U235'=unname(sqrt(E)))
    } else {
        out <- rep(NA,5)
        names(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238','rXY')
        out[1] <- 1/p3507
        out[3] <- 1/p3806
        out[c(2,4)] <- sqrt(diag(E))
        out[5] <- stats::cov2cor(E)[1,2]
    }
    out
}
correct.common.Pb.with.208 <- function(x,i,tt,c0608=NULL,c0708=NULL,cc=FALSE){
    # (3806, 08c06), (3507, 08c07), (3238, 3208), (06c08), (07c08):
    ir <- get.UPb.isochron.ratios.208(x,i,tt=tt)
    if (x$format%in%c(7,8,12)){
        r3507 <- age_to_U235Pb207_ratio(tt,d=x$d[i])[1]
        p3507 <- ir$x['U235Pb207'] + ir$x['Pb208cPb207']*r3507*c0708
    }
    if (x$format%in%c(7,8,11)){
        r3806 <- age_to_U238Pb206_ratio(tt,d=x$d[i])[1]
        p3806 <- ir$x['U238Pb206'] + ir$x['Pb208cPb206']*r3806*c0608
    }
    r3208 <- 1/age_to_Pb208Th232_ratio(tt)[1]
    if (x$format==11){
        p3208 <- ir$x['Th232Pb208'] + ir$x['Pb206cPb208']*r3208/c0608
    } else {
        p3208 <- ir$x['Th232Pb208'] + ir$x['Pb207cPb208']*r3208/c0708
    }
    # projected compositions:
    if (x$format%in%c(7,8)){
        ni <- 4
        p <- c(p3507,p3806,p3208,ir$x['Th232U238'])
    } else if (x$format==11){
        ni <- 1
        p <- c(p3806,p3208)
    } else if (x$format==12){
        ni <- 1
        p <- c(p3507,p3208)
    } else {
        stop('Invalid U-Pb format.')
    }
    if (x$format%in%c(7,8)){
        Jp <- matrix(0,4,8)
        Jp[1,3] <- 1
        Jp[1,4] <- r3507*c0708
        Jp[2,1] <- 1
        Jp[2,2] <- r3806*c0608
        Jp[3,6] <- 1
        Jp[3,8] <- r3208/c0708
        Jp[4,5] <- 1
    } else if (x$format==11){
        Jp <- matrix(0,2,4)
        Jp[1,1] <- 1
        Jp[1,2] <- r3806*c0608
        Jp[2,3] <- 1
        Jp[2,4] <- r3208/c0608
    } else if (x$format==12){
        Jp <- matrix(0,2,4)
        Jp[1,1] <- 1
        Jp[1,2] <- r3507*c0708
        Jp[2,3] <- 1
        Jp[2,4] <- r3208/c0708
    }
    Ep <- Jp %*% ir$cov %*% t(Jp)
    if (x$format%in%c(7,8)){
        J <- diag(4)
        diag(J)[1:3] <- -1/p[1:3]^2
    } else {
        J <- diag(-1/p^2)
    }
    E <- J %*% Ep %*% t(J)
    if (cc){
        out <- list()
        if (x$format%in%c(7,8)){
            cnames <- c('Pb207U235','Pb206U238')
            out$x <- 1/p[1:2]
            out$cov <- E[1:2,1:2]
        } else if (x$format==11){
            cnames <- c('Pb206U238','Pb208Th232')
            out$x <- 1/p
            out$cov <- E
        } else {
            cnames <- c('Pb207U235','Pb208Th232')
            out$x <- 1/p
            out$cov <- E
        }
        names(out$x) <- cnames
        rownames(out$cov) <- cnames
        colnames(out$cov) <- cnames
    } else {
        if (x$format%in%c(7,8)){
            out <- rep(NA,14)
            names(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238',
                            'Pb208Th232','errPb208Th232','Th232U238','errTh232U238',
                            'rXY','rXZ','rXW','rYZ','rYW','rZW')
            out[c(1,3,5,7)] <- c(1/p[1:3],p[4])
            cormat <- matrix(0,4,4)
            pos <- which(diag(E)>0)
            cormat[pos,pos] <- stats::cov2cor(E[pos,pos])
            out[c(2,4,6,8)] <- sqrt(diag(E))
            out[9:11] <- cormat[1,2:4]
            out[12:13] <- cormat[2,3:4]
            out[14] <- cormat[3,4]
        } else { # formats 11 and 12
            out <- rep(NA,5)
            if (x$format==11){
                names(out) <- c('Pb206U238','errPb206U238',
                                'Pb208Th232','errPb208Th232','rXY')
            } else {
                names(out) <- c('Pb207U235','errPb207U235',
                                'Pb208Th232','errPb208Th232','rXY')
            }
            out[c(1,3)] <- 1/p
            cormat <- stats::cov2cor(E)
            out[c(2,4)] <- sqrt(diag(E))
            out[5] <- cormat[1,2]
        }
    } 
    out
}

common.Pb.stacey.kramers <- function(x){
    ns <- length(x)
    if (x$format %in% c(1,2,3)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('U238Pb206','errU238Pb206',
                           'Pb207Pb206','errPb207Pb206','rXY')
        for (i in 1:ns){
            maxt <- get.Pb207Pb206.age(x,i=i)[1]
            tint <- stats::optimise(SKmisfit,interval=c(0,maxt),x=x,i=i)$minimum
            i6474 <- stacey.kramers(tint)
            c76 <- i6474[,'i74']/i6474[,'i64']
            out[i,] <- correct.common.Pb.without.20x(x=x,i=i,c76=c76,tt=tint)
        }
    } else if (x$format %in% c(4,5,6,85)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb206U238','errPb206U238','rXY')
        for (i in 1:ns){
            maxt <- get.Pb207Pb206.age(x,i=i)[1]
            tint <- stats::optimise(SKmisfit,interval=c(0,maxt),x=x,i=i)$minimum
            c6784 <- stacey.kramers(tint)
            if (x$format==85){
                cx6 <- c6784[,'i84']/c6784[,'i64']
                cx7 <- c6784[,'i84']/c6784[,'i74']
            } else {
                cx6 <- 1/c6784[,'i64']
                cx7 <- 1/c6784[,'i74']
            }
            out[i,] <- correct.common.Pb.with.20x(x,i=i,tt=tint,cx6=cx6,cx7=cx7)
        }
    } else if (x$format%in%c(7,8)){
        out <- matrix(0,ns,14)
        colnames(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','Th232U238','errTh232U238',
                           'rXY','rXZ','rXW','rYZ','rYW','rZW')
        for (i in 1:ns){
            maxt <- get.Pb208Th232.age(x,i=i)[1]
            tint <- stats::optimise(SKmisfit,interval=c(0,maxt),x=x,i=i)$minimum
            c678 <- stacey.kramers(tint)
            c68 <- c678[,'i64']/c678[,'i84']
            c78 <- c678[,'i74']/c678[,'i84']
            out[i,] <- correct.common.Pb.with.208(x,i=i,tt=tint,c0608=c68,c0708=c78)
        }
    } else if (x$format%in%c(9,119)){
        out <- matrix(0,ns,2)
        colnames(out) <- c('Pb206U238','errPb206U238')
        tmax <- get.Pb206U238.age(x=x)[,1]
        for (i in 1:ns){
            tint <- stats::uniroot(SKmisfit,interval=c(0,tmax[i]),x=x,i=i)$root
            c6784 <- stacey.kramers(tint)
            cx6 <- ifelse(x$format==119,c6784[,'i84'],1)/c6784[,'i64']
            out[i,] <- correct.common.Pb.with.20x(x,i=i,tt=tint,cx6=cx6)
        }
    } else if (x$format%in%c(10,1210)){
        out <- matrix(0,ns,2)
        colnames(out) <- c('Pb207U235','errPb207U235')
        tmax <- get.Pb207U235.age(x=x)[,1]
        for (i in 1:ns){
            tint <- stats::uniroot(SKmisfit,interval=c(0,tmax[i]),x=x,i=i)$root
            c6784 <- stacey.kramers(tint)
            cx7 <- ifelse(x$format==1210,c6784[,'i84'],1)/c6784[,'i74']
            out[i,] <- correct.common.Pb.with.20x(x,i=i,tt=tint,cx7=cx7)
        }
    } else if (x$format==11){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','rXY')
        tmax <- get.Pb208Th232.age(x=x)[,1]
        for (i in 1:ns){
            tint <- stats::optimise(SKmisfit,interval=c(-1,tmax[i]),x=x,i=i)$minimum
            if (tint>0){
                c678 <- stacey.kramers(tint)
                c0608 <- c678[,'i64']/c678[,'i84']
                out[i,] <- correct.common.Pb.with.208(x,i=i,tt=tint,c0608=c0608)
            } else {
                out[i,] <- NA
            }
        }
    } else if (x$format==12){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb208Th232','errPb208Th232','rXY')
        tmax <- get.Pb208Th232.age(x=x)[,1]
        for (i in 1:ns){
            tint <- stats::optimise(SKmisfit,interval=c(-1,tmax[i]),x=x,i=i)$minimum
            if (tint>0){
                c678 <- stacey.kramers(tint)
                c0708 <- c678[,'i74']/c678[,'i84']
                out[i,] <- correct.common.Pb.with.208(x,i=i,tt=tint,c0708=c0708)
            } else {
                out[i,] <- NA
            }
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
                           'errPb207Pb206','rXY')
        c76 <- fit$par['a0']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.without.20x(x,i=i,c76=c76,tt=tt)
        }
    } else if (x$format %in% c(4,5,6,85)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb206U238','errPb206U238','rXY')
        cx6 <- 1/fit$par['a0']
        cx7 <- 1/fit$par['b0']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.20x(x,i=i,cx6=cx6,cx7=cx7,tt=tt)
        }
    } else if (x$format%in%c(7,8)){
        out <- matrix(0,ns,14)
        colnames(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','Th232U238','errTh232U238',
                           'rXY','rXZ','rXW','rYZ','rYW','rZW')
        c0608 <- fit$par['a0']
        c0708 <- fit$par['b0']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.208(x,i,tt=tt,c0608=c0608,c0708=c0708)
        }
    } else if (x$format%in%c(9,119)){
        out <- matrix(0,ns,2)
        colnames(out) <- c('Pb206U238','errPb206U238')
        cx6 <- 1/fit$par['a0']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.20x(x,i=i,cx6=cx6,tt=tt)
        }
    } else if (x$format%in%c(10,1210)){
        out <- matrix(0,ns,2)
        colnames(out) <- c('Pb207U235','errPb207U235')
        cx7 <- 1/fit$par['b0']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.20x(x,i=i,cx7=cx7,tt=tt)
        }
    } else if (x$format==11){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','rXY')
        c0608 <- fit$par['a0']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.208(x,i,tt=tt,c0608=c0608)
        }
    } else if (x$format==12){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb208Th232','errPb208Th232','rXY')
        c0708 <- fit$par['b0']
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.208(x,i,tt=tt,c0708=c0708)
        }
    } else {
        stop('Invalid U-Pb format.')
    }
    out
}

common.Pb.nominal <- function(x){
    ns <- length(x)
    if (x$format %in% c(1,2,3)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('U238Pb206','errU238Pb206',
                           'Pb207Pb206','errPb207Pb206','rXY')
        c76 <- iratio('Pb207Pb206')[1]
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.without.20x(x=x,i=i,c76=c76)
        }
    } else if (x$format %in% c(4,5,6,85)){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb206U238','errPb206U238','rXY')
        if (x$format==85){
            cx6 <- 1/iratio('Pb206Pb208')[1]
            cx7 <- 1/iratio('Pb207Pb208')[1]
        } else {
            cx6 <- 1/iratio('Pb206Pb204')[1]
            cx7 <- 1/iratio('Pb207Pb204')[1]
        }
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.20x(x=x,i=i,cx6=cx6,cx7=cx7)
        }
    } else if (x$format%in%c(7,8)){
        out <- matrix(0,ns,14)
        colnames(out) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','Th232U238','errTh232U238',
                           'rXY','rXZ','rXW','rYZ','rYW','rZW')
        c0608 <- iratio('Pb206Pb208')[1]
        c0708 <- iratio('Pb207Pb208')[1]
        tmax <- get.Pb208Th232.age(x=x)[,1]
        for (i in 1:ns){
            tint <- stats::optimise(SS.Pb0,interval=c(0,tmax),
                                    c0608=c0608,c0708=c0708,x=x,i=i)$minimum
            out[i,] <- correct.common.Pb.with.208(x,i,tt=tint,c0608=c0608,c0708=c0708)
        }
    } else if (x$format%in%c(9,119)){
        out <- matrix(0,ns,2)
        colnames(out) <- c('Pb206U238','errPb206U238')
        cx6 <- 1/ifelse(x$format==119,
                        iratio('Pb206Pb208')[1],
                        iratio('Pb206Pb204')[1])
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.20x(x=x,i=i,cx6=cx6)
        }
    } else if (x$format%in%c(10,1210)){
        out <- matrix(0,ns,2)
        colnames(out) <- c('Pb207U235','errPb207U235')
        cx7 <- 1/ifelse(x$format==1210,
                        iratio('Pb207Pb208')[1],
                        iratio('Pb207Pb204')[1])
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.20x(x=x,i=i,cx7=cx7)
        }
    } else if (x$format==11){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','rXY')
        c0608 <- iratio('Pb206Pb208')[1]
        tmax <- get.Pb208Th232.age(x=x)[,1]
        for (i in 1:ns){
            tint <- stats::optimise(SS.Pb0,interval=c(-1,tmax[i]),
                                    x=x,i=i,c0608=c0608)$minimum
            if (tint>0){
                out[i,] <- correct.common.Pb.with.208(x=x,i=i,tt=tint,c0608=c0608)
            } else {
                out[i,] <- NA
            }
        }
    } else if (x$format==12){
        out <- matrix(0,ns,5)
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb208Th232','errPb208Th232','rXY')
        c0708 <- iratio('Pb207Pb208')[1]
        tmax <- get.Pb208Th232.age(x=x)[,1]
        for (i in 1:ns){
            tint <- stats::optimise(SS.Pb0,interval=c(-1,tmax[i]),
                                    x=x,i=i,c0708=c0708)$minimum
            if (tint>0){
                out[i,] <- correct.common.Pb.with.208(x=x,i=i,tt=tint,c0708=c0708)
            } else {
                out[i,] <- NA
            }
        }
    } else {
        stop('Invalid U-Pb format')
    }
    out
}

SKmisfit <- function(tt,x,i){
    if (x$format%in%c(1,2,3)){ # sum of squares
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
        xy <- cbind('X'=X,'sX'=sqrt(covmat[1,1]),
                    'Y'=Y,'sY'=sqrt(covmat[2,2]),
                    'rXY'=stats::cov2cor(covmat)[1,2])
        xy.fitted <- get.york.xy(xy,a,b)
        d <- cbind(X-xy.fitted[,"x"],Y-xy.fitted[,"y"])
        out <- as.numeric(d %*% omega %*% t(d))
    } else if (x$format%in%c(4,5,6,85)){ # sum of squares
        c6784 <- stacey.kramers(tt)
        if (x$format==85){
            cx6 <- c6784[1,'i84']/c6784[1,'i64']
            cx7 <- c6784[1,'i84']/c6784[1,'i74']
        } else {
            cx6 <- 1/c6784[1,'i64']
            cx7 <- 1/c6784[1,'i74']
        }
        ccw <- correct.common.Pb.with.20x(x=x,i=i,cx6=cx6,cx7=cx7,tt=tt,cc=TRUE)
        out <- LL.concordia.age(stats::setNames(tt,'t'),ccw,
                                mswd=TRUE,exterr=FALSE,d=x$d[i])
    } else if (x$format%in%c(7,8)){ # sum of squares
        i678 <- stacey.kramers(tt)
        c0608 <- i678[,'i64']/i678[,'i84']
        c0708 <- i678[,'i74']/i678[,'i84']
        out <- SS.Pb0(tt=tt,x=x,i=i,c0608=c0608,c0708=c0708)
    } else if (x$format%in%c(9,119)){ # predicted - observed
        c6784 <- stacey.kramers(tt)
        if (x$format==9){
            a <- cx6 <- 1/c6784[1,'i64']
            Pbx6label <- 'Pb204Pb206'
        } else { # format == 119
            a <- cx6 <- c6784[1,'i84']/c6784[1,'i64']
            Pbx6label <- 'Pb208Pb206'
        }
        p0638 <- correct.common.Pb.with.20x(x=x,i=i,cx6=cx6,tt=tt)[1]
        b <- -p0638*cx6
        out <- a + b*x$x[i,'U238Pb206'] - x$x[i,Pbx6label]
    } else if (x$format%in%c(10,1210)){ # predicted - observed
        c6784 <- stacey.kramers(tt)
        if (x$format==10){
            a <- cx7 <- 1/c6784[1,'i74']
            Pbx7label <- 'Pb204Pb207'
        } else { # format == 1210
            a <- cx7 <- c6784[1,'i84']/c6784[1,'i74']
            Pbx7label <- 'Pb208Pb207'
        }
        p0735 <- correct.common.Pb.with.20x(x=x,i=i,cx7=cx7,tt=tt)[1]
        b <- -p0735*cx7
        out <- a + b*x$x[i,'U235Pb207'] - x$x[i,Pbx7label]
    } else if (x$format==11){ # sum of squares
        i678 <- stacey.kramers(tt)
        c0608 <- i678[,'i64']/i678[,'i84']
        out <- SS.Pb0(tt=tt,x=x,i=i,c0608=c0608)
    } else if (x$format==12){ # sum of squares
        i678 <- stacey.kramers(tt)
        c0708 <- i678[,'i74']/i678[,'i84']
        out <- SS.Pb0(tt=tt,x=x,i=i,c0708=c0708)
    } else {
        stop('Invalid U-Pb format for SSmisfit().')
    }
    out
}

SS.Pb0 <- function(tt,x,i,c0608=NULL,c0708=NULL){
    if (x$format%in%c(7,8)){
        cc <- correct.common.Pb.with.208(x,i=i,tt=tt,c0608=c0608,
                                         c0708=c0708,cc=TRUE)
        out <- LL.concordia.age(stats::setNames(tt,'t'),cc=cc,type=1,
                                mswd=TRUE,exterr=FALSE,d=x$d[i])
    } else if (x$format==11){ 
        McL <- mclean(tt=tt,d=x$d)
        Pb6U8 <- 1/x$x[i,'U238Pb206']
        Pb8Th2 <- x$x[i,'Pb208Pb206']/(x$x[i,'U238Pb206']*x$x[i,'Th232U238'])
        obs <- Pb6U8-c0608*x$x[i,'Th232U238']*(Pb8Th2-McL$Pb208Th232)
        pred <- McL$Pb206U238
        out <- (obs-pred)^2
    } else if (x$format==12){
        McL <- mclean(tt=tt,d=x$d)
        U85 <- iratio('U238U235')[1]
        Pb7U5 <- 1/x$x[i,'U235Pb207']
        Pb8Th2 <- x$x[i,'Pb208Pb207']/(x$x[i,'U235Pb207']*U85*x$x[i,'Th232U238'])
        obs <- Pb7U5-c0708*x$x[i,'Th232U238']*U85*(Pb8Th2-McL$Pb208Th232)
        pred <- McL$Pb207U235
        out <- (obs-pred)^2
    } else {
        stop('Invalid U-Pb format for nominalPb0misfit().')
    }
    out
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
    U238U235 <- iratio('U238U235')[1]
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
    get.search.limit <- function(a,b,d,m,M){
        ttt <- seq(from=m,to=M,length.out=100)
        for (tt in ttt){
            misfit <- intersection.misfit.york(tt,a=a,b=b,d=d)
            if (misfit<0) return(tt)
        }
        return(NA)
    }
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
        from <- 0
        to <- 5000
        if (neg){
            search.range[1] <- from
            search.range[2] <- get.search.limit(a=a,b=b,d=d,m=from,M=to)
            if (!is.na(search.range[2])) go.ahead <- TRUE
        } else {
            tm <- get.search.limit(a=a,b=b,d=d,m=from,M=to)
            tM <- get.search.limit(a=a,b=b,d=d,m=to,M=from)
            if (!is.na(tm) | !is.na(tM))
                go.ahead <- TRUE
            if (is.na(tm) & !is.na(tM))
                search.range <- c(tM,to)
            else if (is.na(tM) & !is.na(tm))
                search.range <- c(from,tm)
            else if (!is.na(tm) & !is.na(tM) & (t68<tm))
                search.range <- c(from,tm)
            else if (!is.na(tm) & !is.na(tM) & (t68>tM))
                search.range <- c(tM,to)
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
