#' @title
#' Linear regression of U-Pb data with correlated errors, taking
#' into account decay constant uncertainties.
#'
#' @description
#' Implements the maximum likelihood algorithm for Total-Pb/U isochron
#' regression of Ludwig (1998) and extends the underlying methodology
#' to accommodate U-Th-Pb data and initial U-series disequilibrium.
#'
#' @details
#' The 3-dimensional regression algorithm of Ludwig and Titterington
#' (1994) was modified by Ludwig (1998) to fit so-called `Total Pb-U
#' isochrons'. These are constrained to a radiogenic endmember
#' composition that falls on the \code{\link{concordia}} line. In its
#' most sophisticated form, this algorithm does not only allow for
#' correlated errors between variables, but also between
#' aliquots. \code{IsoplotR} currently uses this algorithm to
#' propagate decay constant uncertainties in the total Pb-U isochron
#' ages.
#'
#' @param x an object of class \code{UPb}
#' 
#' @param exterr propagate external sources of
#' uncertainty (i.e. decay constants)?
#' 
#' @param model one of three regression models:
#'
#' \code{1}: fit a discordia_line through the data using the maximum
#' likelihood algorithm of Ludwig (1998), which assumes that the
#' scatter of the data is solely due to the analytical
#' uncertainties. In this case, \code{IsoplotR} will either calculate
#' an upper and lower intercept age (for Wetherill concordia), or a
#' lower intercept age and common
#' (\eqn{^{207}}Pb/\eqn{^{206}}Pb)\eqn{_\circ}-ratio intercept (for
#' Tera-Wasserburg). If the p-value for the chi-square test is less
#' than \code{alpha()}, then the analytical uncertainties are augmented
#' by a factor \eqn{\sqrt{MSWD}}.
#'
#' \code{2}: fit a discordia_line ignoring the analytical uncertainties
#'
#' \code{3}: fit a discordia_line using a modified maximum likelihood
#' algorithm that includes accounts for any overdispersion by adding a
#' geological (co)variance term.
#' 
#' @param anchor
#' control parameters to fix the intercept age or common Pb
#' composition of the isochron fit. This can be a scalar or a vector.
#'
#' If \code{anchor[1]=0}: do not anchor the isochron.
#'
#' If \code{anchor[1]=1}: fix the common Pb composition at the values
#' stored in \code{settings('iratio',...)}.
#'
#' If \code{anchor[1]=2}: force the isochron line to intersect the
#' concordia line at an age equal to \code{anchor[2]}.
#'
#' If \code{anchor[1]=3}: anchor the isochron line to the
#' Stacey-Kramers mantle evolution model.
#'
#' @param type only relevant if \code{x$format>3}. Can take on the following
#' values:
#'
#' \code{'joint'} or \code{0}: 3-dimensional isochron regression.
#'
#' \code{1}: 2-dimensional regression of
#' \eqn{{}^{204}}Pb/\eqn{{}^{206}}Pb vs.
#' \eqn{{}^{238}}U/\eqn{{}^{206}}Pb (for U-Pb formats 4, 5 and 6), or
#' of \eqn{{}^{208}}Pb/\eqn{{}^{206}}Pb vs.
#' \eqn{{}^{238}}U/\eqn{{}^{206}}Pb (for U-Pb formats 7 and 8).
#'
#' \code{2}: 2-dimensional regression of
#' \eqn{{}^{204}}Pb/\eqn{{}^{207}}Pb vs.
#' \eqn{{}^{235}}U/\eqn{{}^{207}}Pb (for U-Pb formats 4, 5 and 6), or
#' of \eqn{{}^{208}}Pb/\eqn{{}^{207}}Pb vs.
#' \eqn{{}^{235}}U/\eqn{{}^{207}}Pb (for U-Pb formats 7 and 8).
#'
#' \code{3}: 2-dimensional regression of
#' \eqn{{}^{206}}Pb/\eqn{{}^{208}}Pb vs.
#' \eqn{{}^{232}}Th/\eqn{{}^{208}}Pb (only for U-Pb formats 7 and 8).
#'
#' \code{4}: 2-dimensional regression of
#' \eqn{{}^{207}}Pb/\eqn{{}^{208}}Pb vs.
#' \eqn{{}^{232}}Th/\eqn{{}^{208}}Pb (only for U-Pb formats 7 and 8).
#'
#' @param plot logical. Only relevant for datasets with measured
#'     disequilibrium. If \code{TRUE}, plots the posterior
#'     distribution of the age and initial activity ratios.
#'
#' @param ... optional arguments
#' 
#' @return
#' \describe{
#'
#' \item{par}{a vector with the lower concordia intercept, the common
#' Pb ratios, the dispersion parameter (if \code{model=3}), and the
#' initial \eqn{{}^{234}}U/\eqn{{}^{238}}U and
#' \eqn{{}^{230}}Th/\eqn{{}^{238}}U activity ratio (in the presence of
#' initial disequilibrium).}
#'
#' \item{cov}{the covariance matrix of \code{par}}
#'
#' \item{df}{the degrees of freedom of the model fit (\eqn{n-2} if
#' \code{x$format<4} or \eqn{2n-3} if \code{x$format>3}, where \eqn{n}
#' is the number of aliquots).}
#'
#' \item{mswd}{the mean square of weighted deviates (a.k.a. reduced
#' Chi-square statistic) for the fit.}
#'
#' \item{p.value}{p-value of a Chi-square test for the linear fit}
#'
#' }
#'
#' @examples
#' f <- system.file("UPb4.csv",package="IsoplotR")
#' d <- read.data(f,method="U-Pb",format=4)
#' fit <- ludwig(d)
#' 
#' @references
#' Ludwig, K.R., 1998. On the treatment of concordant uranium-lead
#' ages. Geochimica et Cosmochimica Acta, 62(4), pp.665-676.
#'
#' Ludwig, K.R. and Titterington, D.M., 1994. Calculation of
#' \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#' Cosmochimica Acta, 58(22), pp.5031-5042.
#'
#' @seealso \code{\link{concordia}}, \code{\link{titterington}},
#'     \code{\link{isochron}}
#' @rdname ludwig
#' @export
ludwig <- function(x,...){ UseMethod("ludwig",x) }
#' @rdname ludwig
#' @export
ludwig <- function(x,model=1,anchor=0,exterr=FALSE,
                   type='joint',plot=FALSE,...){
    type <- checkIsochronType(x,type)
    X <- x
    X$d <- mediand(x$d)
    init <- init_ludwig(X,model=model,anchor=anchor,type=type,buffer=2)
    fit <- contingencyfit(par=init$par,fn=LL_ludwig,lower=init$lower,
                          upper=init$upper,x=X,anchor=anchor,
                          type=type,model=model,exterr=exterr)
    fit$cov <- inverthess(fit$hessian)
    if (measured_disequilibrium(X$d) & type%in%c('joint',0,1,3)){
        fit$posterior <- bayeslud(fit,x=X,anchor=anchor,type=type,
                                  model=model,plot=plot,...)
    }
    efit <- exponentiate(fit)
    afit <- anchormerge(efit,X,anchor=anchor,type=type)
    out <- mswd_lud(afit,x=X,exterr=exterr,type=type)
    out$model <- model
    if (model==3){
        out$wtype <- anchor[1]
        if ('w'%in%names(fit$par)){
            w <- unname(afit$par['w'])
            sw <- sqrt(afit$cov['w','w'])
        } else {
            w <- fixDispersion(model=model,format=x$format,
                               anchor=anchor,type=type)
            sw <- 0
        }
        out$disp <- c('w'=w,'s[w]'=sw)
    }
    out
}

anchormerge <- function(fit,x,anchor=0,type='joint'){
    out <- fit
    inames <- names(fit$par)
    onames <- c('t','a0','b0','w')
    if (!x$d$equilibrium){
        dinames <- c('U48i','ThUi','RaUi','PaUi')
        dnames <- c('U48','ThU','RaU','PaU')
        onames <- c(onames,dinames)
    }
    non <- length(onames)
    out$par <- rep(0,non)
    out$cov <- matrix(0,non,non)
    names(out$par) <- rownames(out$cov) <- colnames(out$cov) <- onames
    if ('t'%ni%inames & anchor[1]==2){
        out$par['t'] <- anchor[2]
    } else {
        out$par['t'] <- fit$par['t']
        out$cov['t',inames] <- fit$cov['t',]
        out$cov[inames,'t'] <- fit$cov[,'t']
    }
    if (type%in%c('joint',0,1,3)){
        if ('a0'%ni%inames){
            if (anchor[1]==1){
                if (x$format<4){
                    out$par['a0'] <- iratio('Pb207Pb206')[1]
                } else if (x$format%in%c(4,5,6,9)){
                    out$par['a0'] <- iratio('Pb206Pb204')[1]
                } else if (x$format%in%c(7,8,11,85,119)){
                    out$par['a0'] <- iratio('Pb206Pb208')[1]
                } else {
                    stop('invalid format')
                }
            } else if (anchor[1]==3){
                sk <- stacey.kramers(fit$par['t'])
                if (x$format<4){
                    out$par['a0'] <- sk[1,'i74']/sk[1,'i64']
                } else if (x$format%in%c(4,5,6,9)){
                    out$par['a0'] <- sk[1,'i64']
                } else if (x$format%in%c(7,8,11,85,119)){
                    out$par['a0'] <- sk[1,'i64']/sk[1,'i84']
                }
            } else {
                stop("Can't determine a0 for this dataset.")
            }
        } else {
            out$par['a0'] <- fit$par['a0']
            out$cov['a0',inames] <- fit$cov['a0',]
            out$cov[inames,'a0'] <- fit$cov[,'a0']        
        }
    } else {
        out$par['a0'] <- NA
    }
    if (x$format>3 & type%in%c('joint',0,2,4)){
        if ('b0'%ni%inames){
            if (anchor[1]==1){
                if (x$format%in%c(4,5,6,10)){
                    out$par['b0'] <- iratio('Pb207Pb204')[1]
                } else if (x$format%in%c(7,8,12,85,1210)){
                    out$par['b0'] <- iratio('Pb207Pb208')[1]
                } else {
                    stop('invalid format')
                }
            } else if (anchor[1]==3){
                sk <- stacey.kramers(fit$par['t'])
                if (x$format%in%c(4,5,6,10)){
                    out$par['b0'] <- sk[1,'i74']
                } else if (x$format%in%c(7,8,12,85,1210)){
                    out$par['b0'] <- sk[1,'i74']/sk[1,'i84']
                } else {
                    stop('invalid format')
                }
            } else {
                stop("Can't determine b0 for this dataset.")
            }
        } else {
            out$par['b0'] <- fit$par['b0']
            out$cov['b0',inames] <- fit$cov['b0',]
            out$cov[inames,'b0'] <- fit$cov[,'b0']
        }
    } else {
        out$par['b0'] <- NA
    }
    if ('w'%in%inames){
        out$par['w'] <- fit$par['w']
        out$cov['w',inames] <- fit$cov['w',]
        out$cov[inames,'w'] <- fit$cov[,'w']
    } else {
        out$par['w'] <- NA
    }
    if (!x$d$equilibrium){
        for (i in seq_along(dnames)){
            diname <- dinames[i]
            if (diname%in%inames){
                out$par[diname] <- fit$par[diname]
                out$cov[diname,inames] <- fit$cov[diname,]
                out$cov[inames,diname] <- fit$cov[,diname]
            } else {
                dname <- dnames[i]
                if (x$d[[dname]]$option==2){
                    McL <- mclean(tt=fit$par['t'],d=x$d)
                    out$par[diname] <- McL[[diname]]
                } else {
                    out$par[diname] <- x$d[[dname]]$x
                }
            }
        }
    }
    out
}

init_ludwig <- function(x,model=1,anchor=0,type='joint',buffer=1){
    init <- york2ludwig(x,anchor=anchor,buffer=buffer,type=type,model=model)
    if (x$d$U48$option==2 | x$d$ThU$option==2){
        if ('t'%in%names(init$par)) tt <- exp(init$par['t'])
        else if (anchor[1]==2) tt <- anchor[2]
        else tt <- 0
        McL <- mclean(tt=tt,d=x$d)
    }
    if (type%in%c('joint',0,1,3)){
        if (x$d$U48$option==1 & x$d$U48$sx>0){
            init$par['U48i'] <- x$d$U48$x
        } else if (x$d$U48$option==2 & x$d$U48$sx>0){
            init$par['U48i'] <- max(x$d$buffer,McL$U48i)
        }
        if ('U48i'%in%names(init$par)){
            init$lower['U48i'] <- x$d$U48$m + x$d$buffer
            init$upper['U48i'] <- x$d$U48$M - x$d$buffer
        }
        if (x$d$ThU$option==1 & x$d$ThU$sx>0){
            init$par['ThUi'] <- x$d$ThU$x
        } else if (x$d$ThU$option==2 & x$d$ThU$sx>0){
            init$par['ThUi'] <- max(x$d$buffer,McL$ThUi)
        }
        if ('ThUi'%in%names(init$par)){
            init$lower['ThUi'] <- x$d$ThU$m + x$d$buffer
            init$upper['ThUi'] <- x$d$ThU$M - x$d$buffer
        }
        if (x$d$RaU$option==1 & x$d$RaU$sx>0){
            init$par['RaUi'] <- x$d$RaU$x
            init$lower['RaUi'] <- x$d$RaU$m + x$d$buffer
            init$upper['RaUi'] <- x$d$RaU$M - x$d$buffer
        }
    }
    if (type%in%c('joint',0,2,4)){
        if (x$d$PaU$option==1 & x$d$PaU$sx>0){
            init$par['PaUi'] <- x$d$PaU$x
            init$lower['PaUi'] <- x$d$PaU$m + x$d$buffer
            init$upper['PaUi'] <- x$d$PaU$M - x$d$buffer
        }
    }
    init
}

fixDispersion <- function(model,format,anchor,type){
    out <- 0
    if (model==3 & anchor[1]==1){
        if (format<4){
            out <- iratio('Pb207Pb206')[2]
        } else if (format%in%c(4,5,6,9,10)){
            if (type==1){
                out <- iratio('Pb206Pb204')[2]
            } else if (type==2){
                out <- iratio('Pb207Pb204')[2]
            }
        } else if (format%in%c(7,8,11,12,85,119,1210)){
            if (type==1){
                out <- iratio('Pb206Pb208')[2]
            } else if (type==2){
                out <- iratio('Pb207Pb208')[2]
            } else if (type==3){
                out <- iratio('Pb206Pb208')[2]
            } else if (type==4){
                out <- iratio('Pb207Pb208')[2]
            }
        }
    } else if (model==3 & anchor[1]==2 & length(anchor)>2){
        out <- anchor[3]
    }
    out
}

LL_ludwig <- function(par,x,X=x,model=1,exterr=FALSE,anchor=0,type='joint'){
    pnames <- names(par)
    if ('t' %in% pnames){
        tt <- exp(par['t'])
    } else if (anchor[1]==2){
        tt <- anchor[2]
    } else {
        stop('missing t')
    }
    if (anchor[1]==3){
        sk <- stacey.kramers(tt)
    }
    if ('a0' %in% pnames){
        a0 <- exp(par['a0'])
    } else if (x$format<4){
        a0 <- ifelse(anchor[1]==3,
                     sk[1,'i74']/sk[1,'i64'],
                     iratio('Pb207Pb206')[1])
    } else if (x$format%in%c(4,5,6,9) & type%in%c('joint',0,1)){
        a0 <- ifelse(anchor[1]==3,
                     sk[1,'i64'],
                     iratio('Pb206Pb204')[1])
    } else if (x$format%in%c(7,8,11,85,119) & type%in%c('joint',0,1,3)){
        a0 <- ifelse(anchor[1]==3,
                     sk[1,'i64']/sk[1,'i84'],
                     iratio('Pb206Pb208')[1])
    }
    if ('b0' %in% pnames){
        b0 <- exp(par['b0'])
    } else if (x$format%in%c(4,5,6,10) & type%in%c('joint',0,2)){
        b0 <- ifelse(anchor[1]==3,
                     sk[1,'i74'],
                     iratio('Pb207Pb204')[1])
    } else if (x$format%in%c(7,8,12,85,1210) & type%in%c('joint',0,2,4)){
        b0 <- ifelse(anchor[1]==3,
                     sk[1,'i74']/sk[1,'i84'],
                     iratio('Pb207Pb208')[1])
    }
    for (aname in c('U48','ThU','RaU','PaU')){
        pname <- paste0(aname,'i')
        if (pname%in%pnames){
            X$d[[aname]]$x <- par[pname]
            X$d[[aname]]$sx <- 0
            X$d[[aname]]$option <- 1
        }
    }
    LL <- 0
    if ('w'%in%pnames){
        w <- exp(unname(par['w']))
    } else {
        w <- fixDispersion(model=model,format=x$format,anchor=anchor,type=type)
    }
    if (x$format<4){
        if ('a0'%in%pnames & anchor[1]==1 & iratio('Pb207Pb206')[2]>0){
            LL <- LL - stats::dnorm(x=a0,
                                    mean=iratio('Pb207Pb206')[1],
                                    sd=iratio('Pb207Pb206')[2],
                                    log=TRUE)
        } else if ('t'%in%pnames & anchor[1]==2 & (length(anchor)>2 & anchor[3]>0)){
            LL <- LL - stats::dnorm(x=tt,mean=anchor[2],sd=anchor[3],log=TRUE)
        }
        ta0b0w <- c('t'=unname(tt),'a0'=unname(a0),'w'=w)
    } else if (x$format%in%c(4,5,6,9,10)){
        if (anchor[1]==1){
            if ('a0'%in%pnames & type%in%c('joint',0,1) & iratio('Pb206Pb204')[2]>0){
                LL <- LL - stats::dnorm(x=a0,
                                        mean=iratio('Pb206Pb204')[1],
                                        sd=iratio('Pb206Pb204')[2],
                                        log=TRUE)
            }
            if ('b0'%in%pnames & type%in%c('joint',0,2) & iratio('Pb207Pb204')[2]>0){
                LL <- LL - stats::dnorm(x=b0,
                                        mean=iratio('Pb207Pb204')[1],
                                        sd=iratio('Pb207Pb204')[2],
                                        log=TRUE)
            }
        } else if ('t'%in%pnames & anchor[1]==2 & (length(anchor)>2 & anchor[3]>0)){
            LL <- LL - stats::dnorm(x=tt,mean=anchor[2],sd=anchor[3],log=TRUE)
        }
        if (type%in%c('joint',0)){
            ta0b0w <- c('t'=unname(tt),'a0'=unname(a0),'b0'=unname(b0),'w'=w)
        } else if (type==1){
            ta0b0w <- c('t'=unname(tt),'a0'=unname(a0),'w'=w)
        } else if (type==2){
            ta0b0w <- c('t'=unname(tt),'b0'=unname(b0),'w'=w)
        }
    } else if (x$format%in%c(7,8,11,12,85,119,1210)){
        if (anchor[1]==1){
            if ('a0'%in%pnames & type%in%c('joint',0,1,3) & (iratio('Pb206Pb208')[2]>0)){
                LL <- LL - stats::dnorm(x=a0,
                                        mean=iratio('Pb206Pb208')[1],
                                        sd=iratio('Pb206Pb208')[2],
                                        log=TRUE)
            }
            if ('b0'%in%pnames & type%in%c('joint',0,2,4) & iratio('Pb207Pb208')[2]>0){
                LL <- LL - stats::dnorm(x=b0,
                                        mean=iratio('Pb207Pb208')[1],
                                        sd=iratio('Pb207Pb208')[2],
                                        log=TRUE)
            }
        } else if ('t'%in%pnames & anchor[1]==2 & (length(anchor)>2 & anchor[3]>0)){
            LL <- LL - stats::dnorm(x=tt,mean=anchor[2],sd=anchor[3],log=TRUE)
        }
        if (type%in%c('joint',0)){
            ta0b0w <- c('t'=unname(tt),'a0'=unname(a0),'b0'=unname(b0),'w'=w)
        } else if (type%in%c(1,3)){
            ta0b0w <- c('t'=unname(tt),'a0'=unname(a0),'w'=w)
        } else if (type%in%c(2,4)){
            ta0b0w <- c('t'=unname(tt),'b0'=unname(b0),'w'=w)
        }
    } else {
        stop('Invalid U-Pb format.')
    }
    if (model==2){
        if (type%in%c('joint',0) & x$format>3){
            LL <- LL + LL_ludwig_model2(ta0b0w,x=X,exterr=exterr)
        } else {
            LL <- LL + LL_ludwig_model2_2d(ta0b0w,x=X,exterr=exterr,type=type)
        }
    } else {
        if (type%in%c('joint',0) & x$format>3){
            LL <- LL + data2ludwig(X,ta0b0w,exterr=exterr)$LL
        } else {
            LL <- LL + data2ludwig_2d(ta0b0w,x=X,model=model,exterr=exterr,
                                      type=type,anchor=anchor)$LL
        }
    }
    if ('U48i'%in%pnames & x$d$U48$option==1){
        LL <- LL - stats::dnorm(x=par['U48i'],mean=x$d$U48$x,
                                sd=x$d$U48$sx,log=TRUE)
    } else if (x$d$U48$option==2 & x$d$U48$sx>0){
        pred <- mclean(tt=tt,d=X$d)
        LL <- LL - stats::dnorm(x=pred$U48,mean=x$d$U48$x,
                                sd=x$d$U48$sx,log=TRUE) -
            prior(x=X$d$U48$x,a=X$d$U48)
    }
    if ('ThUi'%in%pnames & x$d$ThU$option==1){
        LL <- LL - stats::dnorm(x=par['ThUi'],mean=x$d$ThU$x,
                                sd=x$d$ThU$sx,log=TRUE)
    } else if (x$d$ThU$option==2 & x$d$ThU$sx>0){
        pred <- mclean(tt=tt,d=X$d)
        LL <- LL - stats::dnorm(x=pred$ThU,mean=x$d$ThU$x,
                                sd=x$d$ThU$sx,log=TRUE) -
            prior(x=X$d$ThU$x,a=X$d$ThU)
    }
    if ('RaUi'%in%pnames & x$d$RaU$sx>0){
        LL <- LL - stats::dnorm(x=par['RaUi'],mean=x$d$RaU$x,
                                sd=x$d$RaU$sx,log=TRUE)
    }
    if ('PaUi'%in%pnames & x$d$PaU$sx>0){
        LL <- LL - stats::dnorm(x=par['PaUi'],mean=x$d$PaU$x,
                                sd=x$d$PaU$sx,log=TRUE)
    }
    LL
}

data2ludwig <- function(x,ta0b0w,exterr=FALSE){
    out <- list()
    U <- iratio('U238U235')[1]
    tt <- ta0b0w['t']
    a0 <- ta0b0w['a0']
    b0 <- ta0b0w['b0']
    model3 <- 'w'%in%names(ta0b0w) & !is.na(ta0b0w['w'])
    ns <- length(x)
    zeros <- rep(0,ns)
    X <- zeros
    Y <- zeros
    K0 <- zeros
    McL <- mclean(tt=tt,d=x$d,exterr=exterr)
    if (x$format%in%c(4,5,6,85)){
        Z <- zeros
        L0 <- zeros
        NP <- 3 # tt, a0, b0
        NR <- 3 # X, Y, Z
    } else if (x$format%in%c(7,8)){
        Z <- zeros
        W <- zeros
        L0 <- zeros
        NP <- 3 # tt, a0, b0
        NR <- 4 # X, Y, Z, W
    } else {
        stop('Incorrect input format.')
    }
    E <- matrix(0,NR*ns+7,NR*ns+7)
    J <- matrix(0,NP*ns,NR*ns+7)
    J[1:(NP*ns),1:(NP*ns)] <- diag(NP*ns)
    nc <- length(McL$ThUi) # nc>1 if each aliquot has its own diseq correction
    j <- 1
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        if (x$format%in%c(4,5,6)){
            Z[i] <- wd$x['Pb204U238']
        } else if (x$format==85){
            Z[i] <- wd$x['Pb208U238']
        } else if (x$format%in%c(7,8)){
            Z[i] <- wd$x['Pb208Th232']
            W[i] <- wd$x['Th232U238']
        }
        ii <- (0:(NR-1))*ns+i
        E[ii,ii] <- wd$cov
        if (nc>1) j <- i
        J[i,NR*ns+2] <- -McL$dPb207U235dl35[j]     #dKdl35
        J[i,NR*ns+5] <- -McL$dPb207U235dl31[j]     #dKdl31
        J[ns+i,NR*ns+1] <- -McL$dPb206U238dl38[j]  #dLdl38
        J[ns+i,NR*ns+3] <- -McL$dPb206U238dl34[j]  #dLdl34
        J[ns+i,NR*ns+6] <- -McL$dPb206U238dl30[j]  #dLdl30
        J[ns+i,NR*ns+7] <- -McL$dPb206U238dl26[j]  #dLdl26
        if (x$format>6) J[2*ns+i,NR*ns+4] <- -McL$dPb208Th232dl32[j] # dMdl32
    }
    E[NR*ns+1:7,NR*ns+1:7] <- getEl()
    EE <- J%*%E%*%t(J)
    if (model3){
        Ew <- getEw2(w=ta0b0w['w'],x=x,McL=McL)
        EE <- EE + Ew
    }
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    i3 <- (2*ns+1):(3*ns)
    O <- blockinverse3x3(AA=EE[i1,i1],BB=EE[i1,i2],CC=EE[i1,i3],
                         DD=EE[i2,i1],EE=EE[i2,i2],FF=EE[i2,i3],
                         GG=EE[i3,i1],HH=EE[i3,i2],II=EE[i3,i3])
    if (x$format%in%c(4,5,6,85)){
        C1 <- U*b0
        C2 <- a0
        C3 <- rep(0,ns)
    } else if (x$format%in%c(7,8)){
        C1 <- U*b0*W
        C2 <- a0*W
        C3 <- -rep(McL$Pb208Th232,ns)
    }
    K0 <- X - C1*Z - McL$Pb207U235
    L0 <- Y - C2*Z - McL$Pb206U238
    AA <- C1*O[i1,i1]*C1 + C1*O[i1,i2]*C2 + C1*O[i1,i3] +
        C2*O[i2,i1]*C1 + C2*O[i2,i2]*C2 + C2*O[i2,i3] +
        O[i3,i1]*C1 + O[i3,i2]*C2 + O[i3,i3]
    BB <- C1*O[i1,i1]%*%K0 + C1*O[i1,i2]%*%L0 + C1*O[i1,i3]%*%C3 +
        C2*O[i2,i1]%*%K0 + C2*O[i2,i2]%*%L0 + C2*O[i2,i3]%*%C3 +
        O[i3,i1]%*%K0 + O[i3,i2]%*%L0 + O[i3,i3]%*%C3
    N <- as.vector(solve(-AA,BB))
    K <- as.vector(K0 + C1*N)
    L <- as.vector(L0 + C2*N)
    M <- N + C3
    KLM <- c(K,L,M)
    out$c0 <- as.vector(Z - N)
    out$SS <- KLM%*%O%*%KLM
    detEE <- determinant(EE,logarithm=TRUE)$modulus
    out$LL <- (NP*ns*log(2*pi) + detEE + out$SS)/2
    out
}

# attributes overdispersion to common Pb (only for 2D isochrons)
getEw1 <- function(w=0,x,McL=mclean(),type,a0=NA,b0=NA,c0){
    format <- x$format
    ns <- length(x)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    J <- matrix(0,2*ns,ns)
    U85 <- iratio('U238U235')[1]
    if (format<4){
        diag(J[i2,i1]) <- -c0/(U85*a0^2)   # dLda0
    } else if (format%in%c(4,5,6,9,10,85,119,1210)){
        if (type==1){
            diag(J[i2,i1]) <- -c0          # dLda0
        } else if (type==2){
            diag(J[i2,i1]) <- -c0          # dLdb0
        } else {
            stop('invalid isochron type')
        }
    } else if (format%in%c(7,8,11,12)){
        ThU <- x$x[,'Th232U238']
        if (type==1){
            diag(J[i2,i1]) <- -c0*ThU      # dLda0
        } else if (type==2){
            diag(J[i2,i1]) <- -c0*ThU*U85  # dLdb0
        } else if (type==3){
            diag(J[i2,i1]) <- -c0/a0^2     # dLda0
        } else if (type==4){
            diag(J[i2,i1]) <- -c0/b0^2     # dLdb0
        } else {
            stop('invalid isochron type')
        }
    } else {
        stop('invalid format')
    }
    J%*%diag(w^2,ns,ns)%*%t(J)
}
# attributes overdispersion to t
getEw2 <- function(w=0,x,McL=mclean(),type='joint'){
    format <- x$format
    ns <- length(x)
    joint <- type%in%c('joint',0)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    if (joint & format%in%c(4:8,85)){ # 3D regression
        J <- matrix(0,3*ns,ns)
        i3 <- (2*ns+1):(3*ns)
    } else { # 2D regression
        J <- matrix(0,2*ns,ns)
    }
    if (joint & format%in%c(4:8,85)){
        diag(J[i1,i1]) <- -McL$dPb207U235dt      # dKdt
        diag(J[i2,i1]) <- -McL$dPb206U238dt      # dLdt
        if (format%in%c(7,8)) {
            diag(J[i3,i1]) <- -McL$dPb208Th232dt # dMdt
        }
    } else if (format<4){
        diag(J[i1,i1]) <- -McL$dPb206U238dt      # dKdt
        diag(J[i2,i1]) <- -McL$dPb207U235dt      # dLdt
    } else if (format%in%c(4,5,6,9,10,85,119,1210)){
        if (type==1){
            diag(J[i2,i1]) <- -McL$dPb206U238dt  # dLdt
        } else if (type==2){
            diag(J[i2,i1]) <- -McL$dPb207U235dt  # dLdt
        } else {
            stop('invalid isochron type')
        }
    } else if (format%in%c(7,8,11,12)){
        ThU <- x$x[,'Th232U238']
        U85 <- iratio('U238U235')[1]
        if (type==1){
            diag(J[i1,i1]) <- -McL$dPb208Th232dt*ThU      # dKdt
            diag(J[i2,i1]) <- -McL$dPb206U238dt           # dLdt
        } else if (type==2){
            diag(J[i1,i1]) <- -McL$dPb208Th232dt*ThU*U85  # dKdt
            diag(J[i2,i1]) <- -McL$dPb207U235dt           # dLdt
        } else if (type==3){
            diag(J[i1,i1]) <- -McL$dPb206U238dt/ThU       # dKdt
            diag(J[i2,i1]) <- -McL$dPb208Th232dt          # dLdt
        } else if (type==4){
            diag(J[i1,i1]) <- -McL$dPb206U238dt/(ThU*U85) # dKdt
            diag(J[i2,i1]) <- -McL$dPb208Th232dt          # dLdt
        } else {
            stop('invalid isochron type')
        }
    } else {
        stop('invalid format')
    }
    J%*%diag(w^2,ns,ns)%*%t(J)
}

LL_ludwig_model2 <- function(ta0b0,x,exterr=FALSE){
    tt <- ta0b0['t']
    a0 <- ta0b0['a0']
    b0 <- ta0b0['b0']
    ns <- length(x)
    McL <- mclean(tt=tt,d=x$d,exterr=exterr)
    if (x$format%in%c(7,8,11,12)){
        XY <- get_UPb_isochron_ratios_208(x,tt=tt)[,1:4]
    } else {
        XY <- get_UPb_isochron_ratios_20x(x)
    }
    X6 <- XY[,1,drop=FALSE] # U238Pb206
    Y6 <- XY[,2,drop=FALSE] # Pb204Pb206 or Pb208cPb206
    X7 <- XY[,3,drop=FALSE] # U235Pb207
    Y7 <- XY[,4,drop=FALSE] # Pb204Pb207 or Pb208cPb207
    a6 <- 1/a0
    a7 <- 1/b0
    b6 <- -a6*McL$Pb206U238
    b7 <- -a7*McL$Pb207U235
    N6 <- (Y6-a6-b6*X6)*b6
    D6 <- 1+b6^2
    x6 <- X6 + N6/D6
    N7 <- (Y7-a7-b7*X7)*b7
    D7 <- 1+b7^2
    x7 <- X7 + N7/D7
    D <- cbind(X6-x6,X7-x7)
    E <- stats::cov(D)
    DD <- c(X6-x6,X7-x7)
    EE <- matrix(0,2*ns,2*ns)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    diag(EE[i1,i1]) <- E[1,1]
    diag(EE[i2,i2]) <- E[2,2]
    diag(EE[i1,i2]) <- diag(EE[i2,i1]) <- E[1,2]
    if (exterr){
        db6dl38 <- -a6*McL$dPb206U238dl38
        db6dl34 <- -a6*McL$dPb206U238dl34
        db6dl30 <- -a6*McL$dPb206U238dl30
        db6dl26 <- -a6*McL$dPb206U238dl26
        db7dl35 <- -a7*McL$dPb207U235dl35
        db7dl31 <- -a7*McL$dPb207U235dl31
        dN6dl38 <-  (Y6-a6)*db6dl38 - 2*X6*b6*db6dl38 # N6 = (Y6-a6)*b6 - X6*b6^2
        dN6dl34 <-  (Y6-a6)*db6dl34 - 2*X6*b6*db6dl34
        dN6dl30 <-  (Y6-a6)*db6dl30 - 2*X6*b6*db6dl30
        dN6dl26 <-  (Y6-a6)*db6dl26 - 2*X6*b6*db6dl26
        dD6dl38 <- 2*b6*db6dl38                       # D6 = 1+b6^2
        dD6dl34 <- 2*b6*db6dl34
        dD6dl30 <- 2*b6*db6dl30
        dD6dl26 <- 2*b6*db6dl26 
        dN7dl35 <-  (Y7-a7)*db7dl35 - 2*X7*b7*db7dl35 # N7 = (Y7-a7)*b7 - X7*b7^2
        dN7dl31 <-  (Y7-a7)*db7dl31 - 2*X7*b7*db7dl31
        dD7dl35 <- 2*b7*db7dl35                       # D7 = 1+b7^2
        dD7dl31 <- 2*b7*db7dl31
        dx6dl38 <- (dN6dl38*D6-N6*dD6dl38)/D6^2       # x6 = N6/D6
        dx6dl34 <- (dN6dl34*D6-N6*dD6dl34)/D6^2
        dx6dl30 <- (dN6dl30*D6-N6*dD6dl30)/D6^2
        dx6dl26 <- (dN6dl26*D6-N6*dD6dl26)/D6^2
        dx7dl35 <- (dN7dl35*D7-N7*dD7dl35)/D7^2       # x7 = N7/D7
        dx7dl31 <- (dN7dl31*D7-N7*dD7dl31)/D7^2
        El <- getEl()
        J <- matrix(0,2*ns,7)
        colnames(J) <- colnames(El)
        J[i1,'U238'] <- -dx6dl38
        J[i1,'U234'] <- -dx6dl34
        J[i1,'Th230'] <- -dx6dl30
        J[i1,'Ra226'] <- -dx6dl26
        J[i2,'U235'] <- -dx7dl35
        J[i2,'Pa231'] <- -dx7dl31
        EE <- EE + J%*%El%*%t(J)
    }
    LL_norm(DD,EE)
}

data2ludwig_2d <- function(ta0b0w,x,model=1,exterr=FALSE,type=1,anchor=0){
    out <- list()
    pnames <- names(ta0b0w)
    tt <- ta0b0w['t']
    McL <- mclean(tt=tt,d=x$d,exterr=exterr)
    a0 <- ifelse('a0'%in%pnames,ta0b0w['a0'],NA)
    b0 <- ifelse('b0'%in%pnames,ta0b0w['b0'],NA)
    ns <- length(x)
    l8 <- lambda('U238')[1]
    l5 <- lambda('U235')[1]
    l2 <- lambda('Th232')[1]
    U85 <- iratio('U238U235')[1]
    multiplier <- 1
    if (x$format<4){ # X=07/35, Y=06/38
        yd <- data2york(x,option=1)
        A <- rep(McL$Pb207U235,ns)
        b <- 1/(a0*U85)
        L0 <- yd[,'Y'] - McL$Pb206U238 - b*yd[,'X']
    } else if (x$format%in%c(4,5,6,9,10,85,119,1210)){
        if (type==1){ # X=0x/38, Y=06/38
            yd <- data2york.UPb(x,option=10)
            A <- rep(0,ns)
            b <- a0
            L0 <- yd[,'Y'] - McL$Pb206U238 - b*yd[,'X']
        } else if (type==2){ # X=0x/35, Y=07/35
            yd <- data2york.UPb(x,option=11)
            A <- rep(0,ns)
            b <- b0
            L0 <- yd[,'Y'] - McL$Pb207U235 - b*yd[,'X']
        } else {
            stop('invalid isochron type')
        }
    } else if (x$format%in%c(7,8,11,12)){
        ThU <- x$x[,'Th232U238']
        if (type==1){ # X=08/38, Y=06/38
            multiplier <- ThU
            yd <- data2york.UPb(x,option=12)
            A <- McL$Pb208Th232*multiplier
            b <- a0
            L0 <- yd[,'Y'] - McL$Pb206U238 - b*yd[,'X']
        } else if (type==2){ # X=08/35, Y=07/35
            multiplier <- ThU*U85
            yd <- data2york.UPb(x,option=13)
            A <- McL$Pb208Th232*multiplier
            b <- b0
            L0 <- yd[,'Y'] - McL$Pb207U235 - b*yd[,'X']
        } else if (type==3){ # X=06/32, Y=08/32
            yd <- data2york.UPb(x,option=14)
            A <- McL$Pb206U238/ThU
            b <- 1/a0
            L0 <- yd[,'Y'] - McL$Pb208Th232 - b*yd[,'X']
        } else if (type==4){ # X=07/32, Y=08/32
            yd <- data2york.UPb(x,option=15)
            A <- McL$Pb206U238/(ThU*U85)
            b <- 1/b0
            L0 <- yd[,'Y'] - McL$Pb208Th232 - b*yd[,'X']
        } else {
            stop('invalid isochron type')
        }
    } else {
        stop('Invalid U-Pb format for data2ludwig_2d')
    }
    E <- matrix(0,nrow=2*ns,ncol=2*ns)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    diag(E)[i1] <- yd[,'sX']^2
    diag(E)[i2] <- yd[,'sY']^2
    E[i1,i2] <- E[i2,i1] <- diag(ns)*yd[,'rXY']*yd[,'sX']*yd[,'sY']
    if (model==3){
        if (anchor[1]==1){
            c0 <- getc0SSLL(yd,A=A,b=b,L0=L0,E=E,multiplier=multiplier)$c0
            Ew <- getEw1(w=ta0b0w['w'],x=x,McL=McL,type=type,a0=a0,b0=b0,c0=c0)
        } else {
            Ew <- getEw2(w=ta0b0w['w'],x=x,McL=McL,type=type)
        }
        E <- E + Ew
    }
    if (exterr){
        El <- getEl()
        J <- matrix(0,2*ns,7)
        colnames(J) <- colnames(El)
        if (x$format<4){
            J[i1,'U235'] <- -McL$dPb207U235dl35
            J[i1,'Pa231'] <- -McL$dPb207U235dl31
            J[i2,'U238'] <- -McL$dPb206U238dl38
            J[i2,'U234'] <- -McL$dPb206U238dl34
            J[i2,'Th230'] <- -McL$dPb206U238dl30
            J[i2,'Ra226'] <- -McL$dPb206U238dl26
        } else if (x$format%in%c(4,5,6,9,10,85,119,1210)){
            if (type==1){
                J[i2,'U238'] <- -McL$dPb206U238dl38
                J[i2,'U234'] <- -McL$dPb206U238dl34
                J[i2,'Th230'] <- -McL$dPb206U238dl30
                J[i2,'Ra226'] <- -McL$dPb206U238dl26
            } else if (type==2){
                J[i2,'U235'] <- -McL$dPb207U235dl35
                J[i2,'Pa231'] <- -McL$dPb207U235dl31                
            } else {
                stop('illegal type')
            }
        } else if (x$format%in%c(7,8)){
            if (type==1){
                J[i1,'Th232'] <- -McL$dPb208Th232dl32*ThU
                J[i2,'U238'] <- -McL$dPb206U238dl38
                J[i2,'U234'] <- -McL$dPb206U238dl34
                J[i2,'Th230'] <- -McL$dPb206U238dl30
                J[i2,'Ra226'] <- -McL$dPb206U238dl26
            } else if (type==2){
                J[i1,'Th232'] <- -McL$dPb208Th232dl32*ThU*U85
                J[i2,'U235'] <- -McL$dPb207U235dl35
                J[i2,'Pa231'] <- -McL$dPb207U235dl31                
            } else if (type==3){
                J[i1,'Th232'] <- -McL$dPb208Th232dl32/ThU
                J[i2,'U238'] <- -McL$dPb206U238dl38
                J[i2,'U234'] <- -McL$dPb206U238dl34
                J[i2,'Th230'] <- -McL$dPb206U238dl30
                J[i2,'Ra226'] <- -McL$dPb206U238dl26                
            } else if (type==4){
                J[i1,'Th232'] <- -McL$dPb208Th232dl32/(ThU*U85)
                J[i2,'U235'] <- -McL$dPb207U235dl35
                J[i2,'Pa231'] <- -McL$dPb207U235dl31                
            } else {
                stop('illegal type')
            }
        } else {
            stop('illegal format')
        }
        E <- E + J %*% El %*% t(J)
    }
    out <- append(out,getc0SSLL(yd,A=A,b=b,L0=L0,E=E,multiplier=multiplier))
    out
}
getc0SSLL <- function(yd,A,b,L0,E,multiplier=1){
    ns <- nrow(E)/2
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    O <- blockinverse(AA=E[i1,i1],BB=E[i1,i2],
                      CC=E[i2,i1],DD=E[i2,i2],doall=TRUE)
    BB <- O[i1,i1] + O[i1,i2]*b + b*O[i2,i1] + b*O[i2,i2]*b
    CC <- O[i1,i2]%*%L0 + b*O[i2,i2]%*%L0 - O[i1,i1]%*%A - b*O[i1,i2]%*%A
    DD <- t(CC)
    M <- -solve(BB,CC)
    K <- M-A
    L <- L0+b*M
    KL <- c(K,L)
    out <- list()
    out$c0 <- (yd[,'X']-M)/multiplier
    out$SS <- KL%*%O%*%KL
    out$LL <- LL_norm(KL,E)
    out
}

LL_ludwig_model2_2d <- function(ta0b0,x,exterr=FALSE,type=1){
    pnames <- names(ta0b0)
    tt <- ta0b0['t']
    if ('a0'%in%pnames) a0 <- ta0b0['a0']
    else if ('b0'%in%pnames) b0 <- ta0b0['b0']
    else stop('neither a0 nor b0 were supplied')
    McL <- mclean(tt=tt,d=x$d,exterr=exterr)
    dbdl38 <- dbdl34 <- dbdl30 <- dbdl26 <- dbdl35 <- dbdl31 <- dbdl32 <- 0
    if (x$format<4){ # X=07/06, Y=38/06
        yd <- data2york(x,option=2,tt=tt)
        a <- a0
        b <- -McL$Pb206U238*(a0-McL$Pb207Pb206)
        dbdl38 <- -McL$dPb206U238dl38*(a0-McL$Pb207Pb206) +
            McL$Pb206U238*McL$dPb207Pb206dl38
        dbdl34 <- -McL$dPb206U238dl34*(a0-McL$Pb207Pb206) +
            McL$Pb206U238*McL$dPb207Pb206dl34
        dbdl30 <- -McL$dPb206U238dl30*(a0-McL$Pb207Pb206) +
            McL$Pb206U238*McL$dPb207Pb206dl30
        dbdl26 <- -McL$dPb206U238dl26*(a0-McL$Pb207Pb206) +
            McL$Pb206U238*McL$dPb207Pb206dl26
        dbdl35 <- McL$Pb206U238*McL$dPb207Pb206dl35
        dbdl31 <- McL$Pb206U238*McL$dPb207Pb206dl31
    } else if (x$format%in%c(4,5,6,9,10,85,119,1210)){
        if (type==1){         # X=38/06, Y=0x/06
            yd <- data2york(x,option=3,tt=tt)
            a <- 1/a0
            b <- -McL$Pb206U238/a0
            dbdl38 <- -McL$dPb206U238dl38/a0
            dbdl34 <- -McL$dPb206U238dl34/a0
            dbdl30 <- -McL$dPb206U238dl30/a0
            dbdl26 <- -McL$dPb206U238dl26/a0
        } else if (type==2){  # X=35/07, Y=0x/07
            yd <- data2york(x,option=4,tt=tt)
            a <- 1/b0
            b <- -McL$Pb207U235/b0
            dbdl35 <- -McL$dPb207U235dl35/b0
            dbdl31 <- -McL$dPb207U235dl31/b0
        } else {
            stop('invalid type')
        }
    } else if (x$format%in%c(7,8)){
        if (type==1){         # X=38/06, Y=08c/06
            yd <- data2york(x,option=6,tt=tt)
            a <- 1/a0
            b <- -McL$Pb206U238/a0
            dbdl38 <- -McL$dPb206U238dl38/a0
            dbdl34 <- -McL$dPb206U238dl34/a0
            dbdl30 <- -McL$dPb206U238dl30/a0
            dbdl26 <- -McL$dPb206U238dl26/a0
        } else if (type==2){  # X=35/07, Y=08c/07
            yd <- data2york(x,option=7,tt=tt)
            a <- 1/b0
            b <- -McL$Pb207U235/b0
            dbdl35 <- -McL$dPb207U235dl35/b0
            dbdl31 <- -McL$dPb207U235dl31/b0
        } else if (type==3){  # X=32/08, Y=06c/08
            yd <- data2york(x,option=8,tt=tt)
            a <- a0
            b <- -McL$Pb208Th232*a0
            dbdl32 <- -McL$dPb208Th232dl32*a0
        } else if (type==4){  # X=32/08, Y=07c/08
            yd <- data2york(x,option=9,tt=tt)
            a <- b0
            b <- -McL$Pb208Th232*b0
            dbdl32 <- -McL$dPb208Th232dl32*b0
        } else {
            stop('invalid type')
        }
    } else {
        stop('illegal format')
    }
    X <- yd[,'X',drop=FALSE]
    Y <- yd[,'Y',drop=FALSE]
    # Deming regression:
    num <- Y-a-b*X
    den <- sqrt(1+b^2)
    sigma <- stats::sd(num/den)
    if (exterr){
        ns <- length(x)
        El <- getEl()
        J <- matrix(0,ns,7)
        colnames(J) <- colnames(El)
        dnumdb <- -X
        ddendb <- -b/den
        dxdb <- (dnumdb*den-num*ddendb)/den^2
        J[,'U238'] <- -dxdb*dbdl38
        J[,'U235'] <- -dxdb*dbdl35
        J[,'U234'] <- -dxdb*dbdl34
        J[,'Th232'] <- -dxdb*dbdl32
        J[,'Pa231'] <- -dxdb*dbdl31
        J[,'Th230'] <- -dxdb*dbdl30
        J[,'Ra226'] <- -dxdb*dbdl26
        D <- as.vector(num/den)
        E <- diag(0,ns,ns)
        diag(E) <- stats::var(D)
        E <- E + J%*%El%*%t(J)
        LL <- LL_norm(D,E)
    } else {
        LL <- -sum(stats::dnorm(x=num/den,mean=0,sd=sigma,log=TRUE))
    }
    LL
}

mswd_lud <- function(fit,x,exterr=FALSE,type='joint'){
    out <- fit
    ns <- length(x)
    np <- sum(diag(fit$cov[1:4,1:4])>0)
    if (x$format%in%c(1,2,3,9,10,119,1210) || type%ni%c('joint',0)){
        out$df <- ns-np
    } else {
        out$df <- 2*ns-np
    }
    X <- x
    if (x$d$U48$option==2) X$d$U48 <- list(x=unname(fit$par['U48i']),option=1)
    if (x$d$ThU$option==2) X$d$ThU <- list(x=unname(fit$par['ThUi']),option=1)
    if (type%in%c('joint',0) && x$format%in%c(4,5,6,7,8,11,12,85)){
        SS <- data2ludwig(X,ta0b0w=fit$par)$SS
    } else {
        SS <- data2ludwig_2d(ta0b0w=fit$par,x=X,type=type,exterr=exterr)$SS
    }
    out <- append(out,getMSWD(SS,out$df))
    out$n <- ns
    out
}

exponentiate <- function(fit){
    out <- fit
    pnames <- names(fit$par)
    totransform <- c('t','a0','b0','w')
    lpnames <- totransform[totransform%in%pnames]
    out$par[lpnames] <- exp(fit$par[lpnames])
    np <- length(lpnames)
    J <- diag(out$par[lpnames],np,np)
    out$cov[lpnames,lpnames] <- J %*% fit$cov[lpnames,lpnames] %*% t(J)
    dimnames(out$cov) <- dimnames(fit$cov)
    out
}
