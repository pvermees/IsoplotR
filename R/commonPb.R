#' Common Pb correction
#'
#' Applies a common-Pb correction to a U-Pb dataset using either the
#' Stacey-Kramers mantle evolution model, isochron regression, or any
#' nominal inital Pb isotope composition.
#'
#' @details \code{IsoplotR} implements six different methods to
#'     correct for the presence of non-radiogenic (`common')
#'     lead. This includes three strategies tailored to datasets that
#'     include \eqn{^{204}}Pb measurements and a further three
#'     strategies for datasets that do not. \eqn{^{204}}Pb is the only
#'     one of lead's four stable isotopes that does not have a
#'     naturally occurring radioactive parent. This makes it very
#'     useful for common-Pb correction:
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
#' Equations can be solved iteratively for \eqn{t} and
#' \eqn{\left[{}^{206}Pb/{}^{204}Pb\right]_\circ}. The
#' \eqn{{}^{207}}Pb/\eqn{{}^{204}}Pb-ratio is corrected in exactly the
#' same way, using \eqn{\left[{}^{207}Pb/{}^{204}Pb\right]_{3.7Ga} =
#' 12.998}.
#'
#' In the absence of \eqn{^{204}}Pb measurements, a \eqn{^{207}}
#' Pb-based common lead correction can be used:
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
#' \enumerate{
#'    \item Stacey-Kramers correction
#'    \item isochron regression
#'    \item nominal common Pb isotope composition
#' }
#'
#' @param omit vector with indices of aliquots that should be omitted
#'     from the isochron regression (only used if \code{option==2})
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
#' lead isotope evolution by a two-stage model. Earth and planetary
#' science letters, 26(2), pp.207-221.
#' @examples
#' data(examples)
#' UPb <- Pb0corr(examples$UPb,option=1)
#' concordia(UPb)
#' # produces identical results as:
#' dev.new()
#' concordia(examples$UPb,common.Pb=1)
#' @export
Pb0corr <- function(x,option=1,omit=NULL){
    ns <- length(x)
    out <- x
    out$x.raw <- x$x
    if (option == 1){
        x.corr <- common.Pb.stacey.kramers(x)
    } else if (option == 2){
        x.corr <- common.Pb.isochron(x,omit=omit)
    } else if (option == 3){
        x.corr <- common.Pb.nominal(x)
    } else {
        return
    }
    if (x$format==1){
        out$x[,c('Pb207U235','errPb207U235',
                 'Pb206U238','errPb206U238','rhoXY')] <- tw2w(x.corr)
    } else if (x$format==2){
        out$x[,c('U238Pb206','errU238Pb206',
                 'Pb207Pb206','errPb207Pb206','rhoXY')] <- x.corr
    } else if (x$format==3){
        out$x[,c('Pb207U235','errPb207U235',
                 'Pb206U238','errPb206U238')] <- tw2w(x.corr)[,1:4]
        out$x[,c('Pb207Pb206','errPb207Pb206')] <- x.corr[,3:4]
    } else if (x$format==4){
        out$x[,c('Pb207U235','errPb207U235',
                 'Pb206U238','errPb206U238')] <- x.corr[,1:4]
        out$x[,'rhoXY'] <- x.corr[,5]
        out$x[,c('Pb204U238','errPb204U238','rhoXZ','rhoYZ')] <- 0
    } else if (x$format==5){
        tw <- w2tw(x.corr)
        out$x[,c('U238Pb206','errU238Pb206',
                 'Pb207Pb206','errPb207Pb206')] <- tw[,1:4]
        out$x[,'rhoXY'] <- tw[,5]
        out$x[,c('Pb204Pb206','errPb204Pb206','rhoXZ','rhoYZ')] <- 0
    } else if (x$format==6){
        tw <- w2tw(x.corr)
        out$x[,c('Pb207U235','errPb207U235',
                 'Pb206U238','errPb206U238')] <- x.corr[,1:4]
        out$x[,c('Pb204U238','errPb204U238')] <- 0
        out$x[,c('Pb207Pb206','errPb207Pb206')] <- tw[,3:4]
        out$x[,c('Pb204Pb207','errPb204Pb207',
                 'Pb204Pb206','errPb204Pb206')] <- 0
    }
    out
}

correct.common.Pb.without.204 <- function(x,i,c76,lower=TRUE,project.err=TRUE){
    tw <- tera.wasserburg(x,i)
    m86 <- tw$x['U238Pb206']
    m76 <- tw$x['Pb207Pb206']
    tint <- project.concordia(m86,m76,c76,d=x$d,lower=lower)
    cctw <- age_to_terawasserburg_ratios(tt=tint,st=0,d=x$d)
    r86 <- cctw$x['U238Pb206']
    r76 <- cctw$x['Pb207Pb206']
    if (project.err){
        f <- (m76-r76)/(c76-r76)
        E <- tw$cov/((1-f)^2)
    } else {
        E <- tw$cov
    }
    sr86 <- sqrt(E[1,1])
    sr76 <- sqrt(E[2,2])
    rho <- stats::cov2cor(E)[1,2]
    out <- c(r86,sr86,r76,sr76,rho)
    out
}
correct.common.Pb.with.204 <- function(x,i,c46,c47){
    ir <- get.UPb.isochron.ratios(x,i) # 68, 46, 75, 47
    m68 <- ir$x['Pb206U238']
    m46 <- ir$x['Pb204Pb206']
    m75 <- ir$x['Pb207U235']
    m47 <- ir$x['Pb204Pb207']    
    r75 <- m75*(1-m47/c47)
    r68 <- m68*(1-m46/c46)
    J <- matrix(0,2,4)
    J[1,3] <- 1-m47/c47
    J[1,4] <- -m75/c47
    J[2,1] <- 1-m46/c46
    J[2,2] <- -m68/c46
    E <- J %*% ir$cov %*% t(J)
    sr75 <- sqrt(E[1,1])
    sr68 <- sqrt(E[2,2])
    rho <- stats::cov2cor(E)[1,2]
    out <- c(r75,sr75,r68,sr68,rho)
    out
}

common.Pb.stacey.kramers <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format < 4){
        for (i in 1:ns){
            tint <- stats::optimise(SS.SK.without.204,
                                    interval=c(0,5000),x=x,i=i)$minimum
            i6474 <- stacey.kramers(tint)
            c76 <- i6474[2]/i6474[1]
            out[i,] <- correct.common.Pb.without.204(x,i,c76,lower=FALSE)
        }
    } else {
        for (i in 1:ns){
            tint <- stats::optimise(SS.SK.with.204,
                                    interval=c(0,5000),x=x,i=i)$minimum
            c6474 <- stacey.kramers(tint)
            c46 <- 1/c6474[1]
            c47 <- 1/c6474[2]
            out[i,] <- correct.common.Pb.with.204(x,i,c46,c47)
        }
    }
    out
}

common.Pb.isochron <- function(x,omit=NULL){
    ns <- length(x)
    calcit <- (1:ns)%ni%omit
    fit <- ludwig(subset(x,subset=calcit))
    out <- matrix(0,ns,5)
    tt <- fit$par[1]
    if (x$format<4){
        rr <- age_to_terawasserburg_ratios(tt,d=x$d)$x
        slope <- (rr['Pb207Pb206']-fit$par['76i'])/rr['U238Pb206']
        m76 <- get.Pb207Pb206.ratios(x)[,1]
        m86 <- get.U238Pb206.ratios(x)[,1]
        c76 <- m76 - slope*m86
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.without.204(x,i,c76[i],lower=TRUE,
                                                     project.err=FALSE)
        }
    } else {
        r68 <- age_to_Pb206U238_ratio(tt=tt,st=0,d=x$d)
        slope.68 <- -r68[1]/fit$par['64i']
        r75 <- age_to_Pb207U235_ratio(tt=tt,st=0,d=x$d)
        slope.75 <- -r75[1]/fit$par['74i']
        for (i in 1:ns){
            rr <- get.UPb.isochron.ratios(x,i)
            m46 <- rr$x['Pb204Pb206']
            m68 <- rr$x['Pb206U238']
            c46 <- m46 - slope.68/m68
            m47 <- rr$x['Pb204Pb207']
            m75 <- rr$x['Pb207U235']
            c47 <- m47 - slope.75/m75
            out[i,] <- correct.common.Pb.with.204(x,i,c46,c47)
        }
    }
    out
}

common.Pb.nominal <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format < 4){
        c76 <- settings('iratio','Pb207Pb206')[1]
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.without.204(x,i,c76,lower=TRUE)
        }
    } else {
        c46 <- 1/settings('iratio','Pb206Pb204')[1]
        c47 <- 1/settings('iratio','Pb207Pb204')[1]
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.204(x,i,c46,c47)
        }
    }
    out
}

SS.SK.without.204 <- function(tt,x,i){
    tw <- tera.wasserburg(x,i)
    X <- tw$x['U238Pb206']
    Y <- tw$x['Pb207Pb206']
    i6474 <- stacey.kramers(tt)
    cct <- age_to_terawasserburg_ratios(tt,st=0,d=x$d)
    a <- i6474[2]/i6474[1] # intercept
    b <- (cct$x['Pb207Pb206']-a)/cct$x['U238Pb206'] # slope
    omega <- solve(tw$cov)
    x.fitted <- (X*omega[1,1]+Y*omega[1,2]-a*omega[1,2])/
                (omega[1,1]+b*omega[1,2])
    y.fitted <- a + b*x.fitted
    d <- cbind(X-x.fitted,Y-y.fitted)
    as.numeric(d %*% omega %*% t(d))
}
SS.SK.with.204 <- function(tt,x,i){
    wi <- wetherill(x,i=i)
    i6474 <- stacey.kramers(tt)
    i64 <- i6474[1]
    i74 <- i6474[2]
    U <- iratio('U238U235')[1]
    ccw <- list(x=rep(0,2),cov=matrix(0,2,2))
    ccw$x[1] <- wi$x['Pb207U235'] - i74*wi$x['Pb204U238']*U
    ccw$x[2] <- wi$x['Pb206U238'] - i64*wi$x['Pb204U238']
    J <- matrix(0,2,3)
    J[1,1] <- 1
    J[1,3] <- -i74*U
    J[2,2] <- 1
    J[2,3] <- -i64
    ccw$cov <- J %*% wi$cov %*% t(J)
    LL.concordia.age(tt,ccw,mswd=FALSE,exterr=FALSE,d=x$d)
}

stacey.kramers <- function(tt,inverse=FALSE){
    nt <- length(tt)
    sk.206.204 <- rep(0,nt)
    sk.207.204 <- rep(0,nt)
    sk.238.204 <- rep(0,nt)
    ti <- rep(0,nt)
    young <- which(tt < 3700)
    old <- which(tt >= 3700)
    sk.206.204[young] <- 11.152
    sk.207.204[young] <- 12.998
    sk.238.204[young] <- 9.74
    ti[young] <- 3700
    sk.206.204[old] <- 9.307
    sk.207.204[old] <- 10.294
    sk.238.204[old] <- 7.19
    ti[old] <- 4570
    U238U235 <- settings('iratio','U238U235')[1]
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    i64 <- sk.206.204 + sk.238.204*(exp(l8*ti)-exp(l8*tt))
    i74 <- sk.207.204 + sk.238.204*(exp(l5*ti)-exp(l5*tt))/U238U235
    if (inverse) out <- cbind(1/i64,i74/i64)
    else out <- cbind(i64,i74)
    out
}

project.concordia <- function(m86,m76,c76,d=diseq(),lower=TRUE){
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
        if (lower){
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
project.concordia.new <- function(m86,m76,c76,d=diseq(),lower=TRUE){
    a <- 1/m86 - m76/(m86*c76)
    b <- 1/(c76*iratio('U238U235')[1])
    fit <- concordia.intersection.ab(a,b,wetherill=TRUE,d=d)
    fit['t[u]']
}
