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
#' \code{1}: fit a discordia line through the data using the maximum
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
#' \code{2}: fit a discordia line ignoring the analytical uncertainties
#'
#' \code{3}: fit a discordia line using a modified maximum likelihood
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
ludwig <- function(x,model=1,anchor=0,exterr=FALSE,type='joint',plot=FALSE,...){
    init <- init.ludwig(x,model=model,anchor=anchor,type=type,buffer=2)
    ctrl <- list()
    fit <- stats::optim(init$par,fn=LL.ludwig,method='L-BFGS-B',
                        lower=init$lower,upper=init$upper,
                        hessian=TRUE,x=x,anchor=anchor,type=type,
                        model=model,exterr=exterr)
    if (fit$convergence>0){
        ctrl <- list(fnscale=1e-15,maxit=1000)
        fit <- stats::optim(init$par,fn=LL.ludwig,method='L-BFGS-B',
                            lower=init$lower,upper=init$upper,control=ctrl,
                            hessian=TRUE,x=x,anchor=anchor,type=type,
                            model=model,exterr=exterr)
        NMfit <- stats::optim(init$par,fn=LL.ludwig,hessian=TRUE,x=x,
                              anchor=anchor,type=type,model=model,exterr=exterr)
        if (fit$value>NMfit$value){
            fit <- NMfit
            warning('L-BFGS-B did not converge. Switched to Nelder-Mead.')
        }
        if (fit$convergence>0){
            warning('ludwig() did not converge.')
        }
    }
    if (fit$convergence>0){
    }
    fit$cov <- inverthess(fit$hessian)
    if (measured.disequilibrium(x$d) && type%in%c('joint',0,1,3)){
        fit$posterior <- bayeslud(fit,x=x,anchor=anchor,type=type,
                                  control=ctrl,model=model,plot=plot)
    }
    efit <- exponentiate(fit)
    afit <- anchormerge(efit,x,anchor=anchor,type=type)
    out <- mswd.lud(afit,x=x,exterr=exterr,type=type)
    out$model <- model
    if (model==3){
        disp <- out$par['w']
        sdisp <- sqrt(out$cov['w','w'])
        out$disp <- c('w'=unname(disp),'s[w]'=unname(sdisp))
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
    if ('t'%ni%inames && anchor[1]==2){
        out$par['t'] <- anchor[2]
    } else {
        out$par['t'] <- fit$par['t']
        out$cov['t',inames] <- fit$cov['t',]
        out$cov[inames,'t'] <- fit$cov[,'t']
    }
    if (type%in%c('joint',0,1,3)){
        if ('a0'%ni%inames && anchor[1]==1){
            if (x$format<4) out$par['a0'] <- iratio('Pb207Pb206')[1]
            else if (x$format<7) out$par['a0'] <- iratio('Pb206Pb204')[1]
            else out$par['a0'] <- 1/iratio('Pb208Pb206')[1]
        } else {
            out$par['a0'] <- fit$par['a0']
            out$cov['a0',inames] <- fit$cov['a0',]
            out$cov[inames,'a0'] <- fit$cov[,'a0']        
        }
    } else {
        out$par['a0'] <- NA
    }
    if (x$format>3 && type%in%c('joint',0,2,4)){
        if ('b0'%ni%inames && anchor[1]==1){
            if (x$format<7) out$par['b0'] <- iratio('Pb207Pb204')[1]
            else out$par['b0'] <- 1/iratio('Pb208Pb207')[1]
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

inithelper <- function(yd,x0=NULL,y0=NULL){
    out <- c('a'=NA,'b'=NA,'x0inv'=NA)
    if (is.null(x0) & is.null(y0)){
        fit <- york(yd)
        out['a'] <- abs(fit$a[1])
        out['b'] <- -abs(fit$b[1])
    } else if (is.null(y0)){ # anchor x
        i <- which.min(yd[,'X'])
        if (x0>yd[i,'X']){
            out['b'] <- yd[i,'Y']/(yd[i,'X']-x0)
        } else {
            out['b'] <- yd[i,'Y']/(0-x0)
        }
        out['a'] <- abs(out['b']*x0)
    } else { # anchor y
        i <- which.min(yd[,'Y'])
        if (y0>yd[i,'Y']){
            out['b'] <- (yd[i,'Y']-y0)/yd[i,'X']
        } else {
            out['b'] <- y0/(0-yd[i,'X'])
        }
        out['a'] <- y0
    }
    out['x0inv'] <- -out['b']/out['a']
    out
}

init.ludwig <- function(x,model=1,anchor=0,type='joint',buffer=1,debug=FALSE){
    if (debug) browser()
    par <- vector()
    if (model==3){
        pilot <- ludwig(x=x,model=1,anchor=anchor,type=type)
        w <- ifelse(pilot$cov['t','t']==0,
                    pilot$par['t']/100,
                    sqrt(pilot$cov['t','t']*pilot$mswd))
    }
    if (x$format<4){
        yd <- data2york(x,option=2)
        if (anchor[1]==1){
            Pb76c <- iratio('Pb207Pb206')[1]
            abx <- inithelper(yd=yd,y0=Pb76c)
            tt <- get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1]
            par['t'] <- log(tt)
            if (iratio('Pb207Pb206')[2]>0) par['a0'] <- log(Pb76c)
        } else if (anchor[1]==2 && length(anchor)>1){
            tt <- anchor[2]
            U8Pb6r <- 1/mclean(tt=tt,d=x$d)$Pb206U238
            abx <- inithelper(yd=yd,x0=U8Pb6r)
            if (length(anchor)>2 && anchor[3]>0) par['t'] <- log(tt)
            par['a0'] <- log(abx['a'])
        } else {
            abx <- inithelper(yd)
            tt <- get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1]
            par['t'] <- log(tt)
            par['a0'] <- log(abx['a'])
        }
    } else if (x$format<7){
        yda <- data2york(x,option=3)
        ydb <- data2york(x,option=4)
        if (anchor[1]==1){
            if (type%in%c('joint',0,1)){
                Pb64c <- iratio('Pb206Pb204')[1]
                abxa <- inithelper(yd=yda,y0=1/Pb64c)
                tt <- get.Pb206U238.age(x=abxa['x0inv'],d=x$d)[1]
                par['t'] <- log(tt)
                if (iratio('Pb206Pb204')[2]>0) par['a0'] <- log(Pb64c)
            }
            if (type%in%c('joint',0,2)){
                Pb74c <- iratio('Pb207Pb204')[1]
                abxb <- inithelper(yd=ydb,y0=1/Pb74c)
            }
            if (type==2){
                tt <- get.Pb206U238.age(x=abxb['x0inv'],d=x$d)[1]
                par['t'] <- log(tt)
            }
            if (type%in%c('joint',0,2)){
                if (iratio('Pb207Pb204')[2]>0) par['b0'] <- log(Pb74c)
            }
        } else if (anchor[1]==2 && length(anchor)>1){
            tt <- anchor[2]
            if ((length(anchor)>2 && anchor[3]>0)) par['t'] <- log(tt)
            if (type%in%c('joint',0,1)){
                Pb6U8r <- mclean(tt=tt)$Pb206U238
                abxa <- inithelper(yd=yda,x0=1/Pb6U8r)
                par['a0'] <- log(1/abxa['a'])
            }
            if (type%in%c('joint',0,2)){
                Pb7U5r <- mclean(tt=tt)$Pb207U235
                abxb <- inithelper(yd=ydb,x0=1/Pb7U5r)
                par['b0'] <- log(1/abxb['a'])
            }
        } else {
            if (type%in%c('joint',0,1)){
                abxa <- inithelper(yd=yda)
                tt <- get.Pb206U238.age(x=abxa['x0inv'],d=x$d)[1]
            }
            if (type%in%c('joint',0,2)){
                abxb <- inithelper(yd=ydb)
            }
            if (type==2){
                tt <- get.Pb207U235.age(x=abxb['x0inv'],d=x$d)[1]
            }
            par['t'] <- log(tt)
            if (type%in%c('joint',0,1)){
                par['a0'] <- log(1/abxa['a'])
            }
            if (type%in%c('joint',0,2)){
                par['b0'] <- log(1/abxb['a'])
            }
        }
    } else { # formats 7 and 8
        if (anchor[1]==1){
            yd <- data2york(x,option=2)
            abx <- inithelper(yd=yd)
            pilott <- get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1]
            if (type==1){ # 0806 vs 38/06
                yd <- data2york(x,option=6,tt=pilott)
                y0 <- iratio('Pb208Pb206')[1]
            } else if (type==2){ # 0807 vs 35/07
                yd <- data2york(x,option=7,tt=pilott)
                y0 <- iratio('Pb208Pb207')[1]
            } else if (type==3){ # 0608 vs 32/08
                yd <- data2york(x,option=8,tt=pilott)
                y0 <- 1/iratio('Pb208Pb206')[1]
            } else if (type==4){ # 0708 vs 32/08
                yd <- data2york(x,option=9,tt=pilott)
                y0 <- 1/iratio('Pb208Pb207')[1]
            } else { # joint, 0 or 1
                y0 <- iratio('Pb207Pb206')[1]
            }
            abx <- inithelper(yd=yd,y0=y0)
            if (type==1){
                par['t'] <- log(get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1])
                if (iratio('Pb208Pb206')[2]>0) par['a0'] <- log(y0)
            } else if (type==2){
                par['t'] <- log(get.Pb207U235.age(x=abx['x0inv'],d=x$d)[1])
                if (iratio('Pb208Pb207')[2]>0) par['b0'] <- log(y0)
            } else if (type==3){
                par['t'] <- log(get.Pb208Th232.age(x=abx['x0inv'],d=x$d)[1])
                if (iratio('Pb208Pb206')[2]>0) par['a0'] <- log(1/y0)
            } else if (type==4){
                par['t'] <- log(get.Pb208Th232.age(x=abx['x0inv'],d=x$d)[1])
                if (iratio('Pb208Pb207')[2]>0) par['b0'] <- log(1/y0)
            } else { # joint, 0 or 1
                par['t'] <- log(get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1])
                if (iratio('Pb208Pb206')[2]>0){
                    par['a0'] <- log(iratio('Pb208Pb206')[1])
                }
                if (iratio('Pb208Pb207')[2]>0){
                    par['b0'] <- log(iratio('Pb208Pb207')[1])
                }
            }
        } else if (anchor[1]==2 && length(anchor)>1){
            tt <- anchor[2]
            if (length(anchor)>2 && anchor[3]>0) par['t'] <- log(tt)
            if (type==2){ # 0807 vs 35/07
                yd <- data2york(x,option=7,tt=tt)
                x0 <- age_to_U235Pb207_ratio(tt)[1]
            } else if (type==3){ # 0608 vs 32/08
                yd <- data2york(x,option=8,tt=tt)
                x0 <- 1/age_to_Pb208Th232_ratio(tt)[1]
            } else if (type==4){ # 0708 vs 32/08
                yd <- data2york(x,option=9,tt=tt)
                x0 <- 1/age_to_Pb208Th232_ratio(tt)[1]
            } else { # joint, 0 or 1: 0806 vs 38/06
                yd <- data2york(x,option=6,tt=tt)
                x0 <- age_to_U238Pb206_ratio(tt)[1]
            }
            abx <- inithelper(yd=yd,x0=x0)
            if (type==1){
                par['a0'] <- log(1/abx['a'])
            } else if (type==2){
                par['b0'] <- log(1/abx['a'])
            } else if (type==3){
                par['a0'] <- log(abx['a'])
            } else if (type==4){
                par['b0'] <- log(abx['a'])
            } else {
                par['a0'] <- log(1/abx['a'])
                ydb <- data2york(x,option=7)
                abxb <- inithelper(yd=ydb,x0=x0)
                par['b0'] <- log(1/abxb['a'])
            }            
        } else {
            yd <- data2york(x,option=2)
            abx <- inithelper(yd=yd)
            pilott <- get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1]
            if (type==1){ # 0806 vs 38/06
                yd <- data2york(x,option=6,tt=pilott)
            } else if (type==2){ # 0807 vs 35/07
                yd <- data2york(x,option=7,tt=pilott)
            } else if (type==3){ # 0608 vs 32/08
                yd <- data2york(x,option=8,tt=pilott)
            } else if (type==4){ # 0708 vs 32/08
                yd <- data2york(x,option=9,tt=pilott)
            } else { # joint, 0 or 1
                # keep yd
            }
            abx <- inithelper(yd=yd)
            if (type==1){ # 0806 vs 38/06
                par['t'] <- log(get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1])
                par['a0'] <- log(1/abx['a'])
            } else if (type==2){ # 0807 vs 35/07
                par['t'] <- log(get.Pb207U235.age(x=abx['x0inv'],d=x$d)[1])
                par['b0'] <- log(1/abx['a'])
            } else if (type==3){ # 0608 vs 32/08
                par['t'] <- log(get.Pb208Th232.age(x=abx['x0inv'],d=x$d)[1])
                par['a0'] <- log(abx['a'])
            } else if (type==4){ # 0708 vs 32/08
                par['t'] <- log(get.Pb208Th232.age(x=abx['x0inv'],d=x$d)[1])
                par['b0'] <- log(abx['a'])
            } else { # joint, 0 or 1
                par['t'] <- log(get.Pb206U238.age(x=abx['x0inv'],d=x$d)[1])
                yda <- data2york(x,option=6)
                ydb <- data2york(x,option=7)
                abxa <- inithelper(yd=yda)
                abxb <- inithelper(yd=ydb)
                par['a0'] <- log(1/abxa['a'])
                par['b0'] <- log(1/abxb['a'])
            }            
        }
    }
    if (model==3) par['w'] <- log(w)
    if (x$d$U48$option==2 || x$d$ThU$option==2){
        McL <- mclean(tt=tt,d=x$d)
    }
    lower <- par - buffer
    upper <- par + buffer
    if (type%in%c('joint',0,1,3)){
        if (x$d$U48$option==1 && x$d$U48$sx>0){
            par['U48i'] <- x$d$U48$x
        } else if (x$d$U48$option==2 && x$d$U48$sx>0){
            par['U48i'] <- max(0,McL$U48i)
        }
        if ('U48i'%in%names(par)){
            lower['U48i'] <- x$d$U48$m + x$d$buffer
            upper['U48i'] <- x$d$U48$M - x$d$buffer
        }
        if (x$d$ThU$option==1 && x$d$ThU$sx>0){
            par['ThUi'] <- x$d$ThU$x
        } else if (x$d$ThU$option==2 && x$d$ThU$sx>0){
            par['ThUi'] <- max(0,McL$ThUi)
        }
        if ('ThUi'%in%names(par)){
            lower['ThUi'] <- x$d$ThU$m + x$d$buffer
            upper['ThUi'] <- x$d$ThU$M - x$d$buffer
        }
        if (x$d$RaU$option==1 && x$d$RaU$sx>0){
            par['RaUi'] <- x$d$RaU$x
            lower['RaUi'] <- x$d$RaU$m + x$d$buffer
            upper['RaUi'] <- x$d$RaU$M - x$d$buffer
        }
    }
    if (type%in%c('joint',0,2,4)){
        if (x$d$PaU$option==1 && x$d$PaU$sx>0){
            par['PaUi'] <- x$d$PaU$x
            lower['PaUi'] <- x$d$PaU$m + x$d$buffer
            upper['PaUi'] <- x$d$PaU$M - x$d$buffer
        }
    }
    list(par=par,lower=lower,upper=upper)
}

LL.ludwig <- function(par,x,X=x,model=1,exterr=FALSE,
                      anchor=0,type='joint',debug=FALSE){
    if (debug) browser()
    pnames <- names(par)
    if ('t' %in% pnames){
        tt <- exp(par['t'])
    } else if (anchor[1]==2){
        tt <- anchor[2]
    } else {
        stop('missing t')
    }
    if ('a0' %in% pnames){
        a0 <- exp(par['a0'])
    } else if (x$format<4){
        a0 <- iratio('Pb207Pb206')[1]
    } else if (x$format<7 && type%in%c('joint',0,1)){
        a0 <- iratio('Pb206Pb204')[1]
    } else if (type%in%c('joint',0,1,3)){
        a0 <- 1/iratio('Pb208Pb206')[1]
    }
    if ('b0' %in% pnames){
        b0 <- exp(par['b0'])
    } else if (x$format%in%c(4,5,6) && type%in%c('joint',0,2)){
        b0 <- iratio('Pb207Pb204')[1]
    } else if (x$format%in%c(7,8) && type%in%c('joint',0,2,4)){
        b0 <- 1/iratio('Pb208Pb207')[1]
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
    if (x$format<4){
        if (anchor[1]==1 && iratio('Pb207Pb206')[2]>0){
            LL <- LL - stats::dnorm(x=a0,
                                    mean=iratio('Pb207Pb206')[1],
                                    sd=iratio('Pb207Pb206')[2],
                                    log=TRUE)
        } else if (anchor[1]==2 &&
                   (length(anchor)>2 && anchor[3]>0)){
            LL <- LL - stats::dnorm(x=tt,mean=anchor[2],sd=anchor[3],log=TRUE)
        }
        ta0b0w <- c('t'=unname(tt),'a0'=unname(a0))
    } else if (x$format<7){
        if (anchor[1]==1){
            if (type%in%c('joint',0,1) && iratio('Pb206Pb204')[2]>0){
                LL <- LL - stats::dnorm(x=a0,
                                        mean=iratio('Pb206Pb204')[1],
                                        sd=iratio('Pb206Pb204')[2],
                                        log=TRUE)
            }
            if (type%in%c('joint',0,2) && iratio('Pb207Pb204')[2]>0){
                LL <- LL - stats::dnorm(x=b0,
                                        mean=iratio('Pb207Pb204')[1],
                                        sd=iratio('Pb207Pb204')[2],
                                        log=TRUE)
            }
        } else if (anchor[1]==2 && (length(anchor)>2 && anchor[3]>0)){
            LL <- LL - stats::dnorm(x=tt,mean=anchor[2],sd=anchor[3],log=TRUE)
        }
        if (type%in%c('joint',0)){
            ta0b0w <- c('t'=unname(tt),'a0'=unname(a0),'b0'=unname(b0))
        } else if (type==1){
            ta0b0w <- c('t'=unname(tt),'a0'=unname(a0))
        } else if (type==2){
            ta0b0w <- c('t'=unname(tt),'b0'=unname(b0))
        }
    } else {
        if (anchor[1]==1){
            if (type%in%c('joint',0,1,3) && (iratio('Pb208Pb206')[2]>0)){
                LL <- LL - stats::dnorm(x=1/a0,
                                        mean=iratio('Pb208Pb206')[1],
                                        sd=iratio('Pb208Pb206')[2],
                                        log=TRUE)
            }
            if (type%in%c('joint',0,2,4) && iratio('Pb208Pb207')[2]>0){
                LL <- LL - stats::dnorm(x=1/b0,
                                        mean=iratio('Pb208Pb207')[1],
                                        sd=iratio('Pb208Pb207')[2],
                                        log=TRUE)
            }
        } else if (anchor[1]==2 && (length(anchor)>2 && anchor[3]>0)){
            LL <- LL - stats::dnorm(x=tt,mean=anchor[2],sd=anchor[3],log=TRUE)
        }
        if (type%in%c('joint',0)){
            ta0b0w <- c('t'=unname(tt),'a0'=unname(a0),'b0'=unname(b0))
        } else if (type%in%c(1,3)){
            ta0b0w <- c('t'=unname(tt),'a0'=unname(a0))
        } else if (type%in%c(2,4)){
            ta0b0w <- c('t'=unname(tt),'b0'=unname(b0))
        }
    }
    if (model==3){
        ta0b0w['w'] <- exp(par['w'])
    }
    if (model==2 && (type%in%c('joint',0) || x$format<4)){
        LL <- LL + LL.ludwig.model2(ta0b0w,x=X,exterr=exterr)
    } else if (type%in%c('joint',0) || x$format<4){
        LL <- LL + data2ludwig(X,ta0b0w,exterr=exterr)$LL
    } else {
        LL <- LL + LL.ludwig.2d(ta0b0w,x=X,model=model,exterr=exterr,type=type)
    }
    if (x$d$U48$option==1 && 'U48i'%in%pnames){
        LL <- LL - stats::dnorm(x=par['U48i'],mean=x$d$U48$x,
                                sd=x$d$U48$sx,log=TRUE)
    } else if (x$d$U48$option==2 && x$d$U48$sx>0){
        pred <- mclean(tt=tt,d=X$d)
        LL <- LL - stats::dnorm(x=pred$U48,mean=x$d$U48$x,
                                sd=x$d$U48$sx,log=TRUE) -
            prior(x=X$d$U48$x,a=X$d$U48)
    }
    if (x$d$ThU$option==1 && 'ThUi'%in%pnames){
        LL <- LL - stats::dnorm(x=par['ThUi'],mean=x$d$ThU$x,
                                sd=x$d$ThU$sx,log=TRUE)
    } else if (x$d$ThU$option==2 && x$d$ThU$sx>0){
        pred <- mclean(tt=tt,d=X$d)
        LL <- LL - stats::dnorm(x=pred$ThU,mean=x$d$ThU$x,
                                sd=x$d$ThU$sx,log=TRUE) -
            prior(x=X$d$ThU$x,a=X$d$ThU)
    }
    if ('RaUi'%in%pnames && x$d$RaU$sx>0){
        LL <- LL - stats::dnorm(x=par['RaUi'],mean=x$d$RaU$x,
                                sd=x$d$RaU$sx,log=TRUE)
    }
    if ('PaUi'%in%pnames && x$d$PaU$sx>0){
        LL <- LL - stats::dnorm(x=par['PaUi'],mean=x$d$PaU$x,
                                sd=x$d$PaU$sx,log=TRUE)
    }
    LL
}

data2ludwig <- function(x,ta0b0w,exterr=FALSE,debug=FALSE){
    if (debug) browser()
    out <- list()
    U <- iratio('U238U235')[1]
    tt <- ta0b0w['t']
    a0 <- ta0b0w['a0']
    if (x$format>3) b0 <- ta0b0w['b0']
    ns <- length(x)
    zeros <- rep(0,ns)
    X <- zeros
    Y <- zeros
    K0 <- zeros
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    if (x$format%in%c(1,2,3)){
        NP <- 2 # number of fit parameters (tt, a0)
        NR <- 2 # number of isotopic ratios (X, Y)
    } else if (x$format%in%c(4,5,6)){
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
    nc <- length(D$ThUi) # nc>1 if each aliquot has its own diseq correction
    j <- 1
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        if (x$format%in%c(4,5,6)){
            Z[i] <- wd$x['Pb204U238']
        } else if (x$format>6){
            Z[i] <- wd$x['Pb208Th232']
            W[i] <- wd$x['Th232U238']
        }
        E[(0:(NR-1))*ns+i,(0:(NR-1))*ns+i] <- wd$cov
        if (nc>1) j <- i
        J[i,NR*ns+2] <- -D$dPb207U235dl35[j]     #dKdl35
        J[i,NR*ns+5] <- -D$dPb207U235dl31[j]     #dKdl31
        J[ns+i,NR*ns+1] <- -D$dPb206U238dl38[j]  #dLdl38
        J[ns+i,NR*ns+3] <- -D$dPb206U238dl34[j]  #dLdl34
        J[ns+i,NR*ns+6] <- -D$dPb206U238dl30[j]  #dLdl30
        J[ns+i,NR*ns+7] <- -D$dPb206U238dl26[j]  #dLdl26
        if (x$format>6) J[2*ns+i,NR*ns+4] <- -D$dPb208Th232dl32[j] # dMdl32
    }
    E[NR*ns+1:7,NR*ns+1:7] <- getEl()
    ED <- J%*%E%*%t(J)
    if ('w'%in%names(ta0b0w) && !is.na(ta0b0w['w'])){
        Ew <- get.Ewd(w=ta0b0w['w'],format=x$format,ns=ns,D=D)
        ED <- ED + Ew
    }
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    if (x$format<4){
        O <- blockinverse(AA=ED[i1,i1],BB=ED[i1,i2],
                          CC=ED[i2,i1],DD=ED[i2,i2],doall=TRUE)
    } else {
        i3 <- (2*ns+1):(3*ns)
        O <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                             DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                             GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    }
    if (x$format%in%c(1,2,3)){
        K0 <- X - D$Pb207U235 + a0*U*(D$Pb206U238 - Y)
        A <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*a0*U +
               K0%*%(O[i1,i2]+t(O[i2,i1])))
        B <- -(a0*U*(O[i1,i1]+t(O[i1,i1]))*a0*U +
               (O[i2,i2]+t(O[i2,i2])) +
               a0*U*(O[i1,i2]+t(O[i1,i2])) +
               (O[i2,i1]+t(O[i2,i1]))*a0*U)
        L <- as.vector(solve(B,A))
        c0 <- Y - D$Pb206U238 - L
        K <- X - D$Pb207U235 - a0*U*c0
        KLM <- c(K,L)
    } else if (x$format%in%c(4,5,6)){
        K0 <- X - D$Pb207U235 - U*b0*Z
        L0 <- Y - D$Pb206U238 - a0*Z
        V <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*U*b0 +
               L0%*%(O[i1,i2]+t(O[i2,i1]))*U*b0 +
               K0%*%(O[i1,i2]+t(O[i2,i1]))*a0 +
               L0%*%(O[i2,i2]+t(O[i2,i2]))*a0 +
               K0%*%(O[i1,i3]+t(O[i3,i1])) +
               L0%*%(O[i2,i3]+t(O[i3,i2])))
        W <- -(U*b0*(O[i1,i1]+t(O[i1,i1]))*U*b0 +
               U*b0*(O[i1,i2]+t(O[i1,i2]))*a0 +
               U*b0*(O[i1,i3]+t(O[i1,i3])) +
               a0*(O[i2,i1]+t(O[i2,i1]))*U*b0 +
               a0*(O[i2,i2]+t(O[i2,i2]))*a0 +
               a0*(O[i2,i3]+t(O[i2,i3])) +
               (O[i3,i1]+t(O[i3,i1]))*U*b0 +
               (O[i3,i2]+t(O[i3,i2]))*a0 +
               (O[i3,i3]+t(O[i3,i3])))
        M <- as.vector(solve(W,V))
        c0 <- as.vector(Z - M)
        K <- as.vector(X - D$Pb207U235 - U*b0*c0)
        L <- as.vector(Y - D$Pb206U238 - a0*c0)
        KLM <- c(K,L,M)
    } else if (x$format%in%c(7,8)){
        Wd <- diag(W)
        K0 <- X - D$Pb207U235 - (Z-D$Pb208Th232)*U*W*b0
        L0 <- Y - D$Pb206U238 - (Z-D$Pb208Th232)*W*a0
        AA <- (Wd%*%O[i1,i1]%*%Wd)*(U*b0)^2 +
            (Wd%*%O[i2,i2]%*%Wd)*a0^2 + O[i3,i3] +
            U*a0*b0*Wd%*%(O[i1,i2]+O[i2,i1])%*%Wd +
            U*b0*(Wd%*%O[i1,i3]+O[i3,i1]%*%Wd) +
            a0*(Wd%*%O[i2,i3]+O[i3,i2]%*%Wd)
        BT <- t(U*b0*K0%*%O[i1,i1]%*%Wd +
                a0*L0%*%O[i2,i2]%*%Wd +
                a0*K0%*%O[i1,i2]%*%Wd +
                U*b0*L0%*%O[i2,i1]%*%Wd +
                K0%*%O[i1,i3] + L0%*%O[i2,i3])
        CC <- U*b0*Wd%*%O[i1,i1]%*%K0 +
            a0*Wd%*%O[i2,i2]%*%L0 +
            a0*Wd%*%O[i2,i1]%*%K0 +
            U*b0*Wd%*%O[i1,i2]%*%L0 +
            O[i3,i1]%*%K0 + O[i3,i2]%*%L0
        M <- as.vector(solve(-(AA+t(AA)),(BT+CC)))
        c0 <- as.vector(Z - D$Pb208Th232 - M)
        K <- as.vector(X - D$Pb207U235 - c0*b0*U*W)
        L <- as.vector(Y - D$Pb206U238 - c0*a0*W)
        KLM <- c(K,L,M)
    }
    out$c0 <- c0
    out$SS <- KLM%*%O%*%KLM
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- (NP*ns*log(2*pi) + detED + out$SS)/2
    out
}

get.Ewd <- function(w=0,format=1,ns=1,D=mclean()){
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    if (format<4) {
        J <- matrix(0,2*ns,ns)
    } else {
        J <- matrix(0,3*ns,ns)
        i3 <- (2*ns+1):(3*ns)
    }
    diag(J[i1,i1]) <- -D$dPb207U235dt      # dKdt
    diag(J[i2,i1]) <- -D$dPb206U238dt      # dLdt
    if (format>6) {
        diag(J[i3,i1]) <- -D$dPb208Th232dt # dMdt
    }
    dEdx <- w^2
    dEdx*J%*%t(J)
}

LL.ludwig.model2 <- function(ta0b0,x,exterr=FALSE){
    pnames <- names(ta0b0)
    tt <- ta0b0['t']
    a0 <- ta0b0['a0']
    if ('b0'%in%pnames) b0 <- ta0b0['b0']
    nn <- length(x)
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    if (x$format<4){
        xy <- data2york(x,option=2)[,c('X','Y'),drop=FALSE]
        xr <- 1/D$Pb206U238
        yr <- D$Pb207Pb206
        xx <- xy[,'X',drop=FALSE]
        yy <- xy[,'Y',drop=FALSE]
        dem <- deming(a=a0,b=(yr-a0)/xr,x=xx,y=yy)
        SS <- sum(dem$d^2)
        if (exterr){
            dbdxr <- (a0-yr)/xr^2
            dbdyr <- 1/xr
            dxrdPbU <- -xr^2
            dyrdPbPb <- 1
            dddPbU <- dem$dddb*dbdxr*dxrdPbU
            dddPbPb <- dem$dddb*dbdyr*dyrdPbPb
            dPbUdl <- dPbPbdl <- rep(0,7)
            dPbUdl[1] <- D$dPb206U238dl38
            dPbUdl[3] <- D$dPb206U238dl34
            dPbUdl[6] <- D$dPb206U238dl30
            dPbUdl[7] <- D$dPb206U238dl26
            dPbPbdl[1] <- D$dPb207Pb206dl38
            dPbPbdl[2] <- D$dPb207Pb206dl35
            dPbPbdl[3] <- D$dPb207Pb206dl34
            dPbPbdl[5] <- D$dPb207Pb206dl31
            dPbPbdl[6] <- D$dPb207Pb206dl30
            dPbPbdl[7] <- D$dPb207Pb206dl26
            J <- dddPbU%*%dPbUdl + dddPbPb%*%dPbPbdl
            covmat <- diag(SS/(nn-2),nn,nn) + J%*%getEl()%*%t(J)
            LL <- LL.norm(as.vector(dem$d),covmat)
        } else {
            LL <- SS2LL(SS=SS,nn=nn)
        }
    } else {
        if (x$format<7){
            xy <- get.UPb.isochron.ratios.204(x)
        } else {
            xy <- get.UPb.isochron.ratios.208(x,tt=tt)[,1:4]
        }
        x6 <- xy[,1,drop=FALSE] # U238Pb206
        y6 <- xy[,2,drop=FALSE] # Pb204Pb206 or Pb208cPb206
        x7 <- xy[,3,drop=FALSE] # U235Pb207
        y7 <- xy[,4,drop=FALSE] # Pb204Pb207 or Pb208cPb207
        r86 <- 1/D$Pb206U238
        r57 <- 1/D$Pb207U235
        dem6 <- deming(a=1/a0,b=-1/(a0*r86),x=x6,y=y6)
        dem7 <- deming(a=1/b0,b=-1/(b0*r57),x=x7,y=y7)
        SS6 <- sum(dem6$d^2)
        SS7 <- sum(dem7$d^2)
        if (exterr){
            dbd86 <- 1/(a0*r86^2)
            dbd57 <- 1/(b0*r57^2)
            dr68d86 <- -r86^2
            dr57d75 <- -r57^2
            dd6d68 <- dem6$dddb*dbd86*dr68d86
            dd7d75 <- dem7$dddb*dbd57*dr57d75
            d68dl <- matrix(0,1,7)
            d68dl[1] <- D$dPb206U238dl38
            d68dl[3] <- D$dPb206U238dl34
            d68dl[6] <- D$dPb206U238dl30
            d68dl[7] <- D$dPb206U238dl26
            d75dl <- matrix(0,1,7)
            d75dl[2] <- D$dPb207U235dl35
            d75dl[5] <- D$dPb207U235dl31
            J6 <- dd6d68%*%d68dl
            J7 <- dd7d75%*%d75dl
            E6 <- diag(SS6/(nn-2),nn,nn) + J6%*%getEl()%*%t(J6)
            E7 <- diag(SS7/(nn-2),nn,nn) + J7%*%getEl()%*%t(J7)
            LL <- LL.norm(as.vector(dem6$d),E6) + LL.norm(as.vector(dem7$d),E7)
        } else {
            LL <- SS2LL(SS=SS6+SS7,nn=2*nn)
        }
    }
    LL
}

LL.ludwig.2d <- function(ta0b0w,x,model=1,exterr=FALSE,type=1,LL=TRUE){
    pnames <- names(ta0b0w)
    tt <- ta0b0w['t']
    McL <- mclean(tt=tt,d=x$d,exterr=exterr)
    if (x$format %in% (4:6)){
        if (type==1){
            O <- data2york(x,option=3)
            a <- 1/ta0b0w['a0']
            r68 <- McL$Pb206U238
            b <- -a*r68
            dbdt <- -a*McL$dPb206U238dt
        } else {
            O <- data2york(x,option=4)
            a <- 1/ta0b0w['b0']
            r75 <- McL$Pb207U235
            b <- -a*r75
            dbdt <- -a*McL$dPb207U235dt
        }
    } else if (x$format %in% (7:8)){
        if (type==1){
            O <- data2york(x,option=6,tt=tt)
            a <- 1/ta0b0w['a0']
            r68 <- McL$Pb206U238
            b <- -a*r68
            dbdt <- -a*McL$dPb206U238dt
        } else if (type==2){
            O <- data2york(x,option=7,tt=tt)
            a <- 1/ta0b0w['b0']
            r75 <- McL$Pb207U235
            b <- -a*r75
            dbdt <- -a*McL$dPb207U235dt
        } else if (type==3){
            O <- data2york(x,option=8,tt=tt)
            a <- ta0b0w['a0']
            r82 <- age_to_Pb208Th232_ratio(tt)[1]
            b <- -a*r82
            dbdt <- -a*McL$dPb208Th232dt
        } else if (type==4){
            O <- data2york(x,option=9,tt=tt)
            a <- ta0b0w['b0']
            r82 <- age_to_Pb208Th232_ratio(tt)[1]
            b <- -a*r82
            dbdt <- -a*McL$dPb208Th232dt
        }
    } else {
        stop('LL.ludwig.2d is only relevant to U-Pb formats 4-8')
    }
    ns <- length(x)
    if (model==2){
        dem <- deming(a=a,b=b,x=O[,'X'],y=O[,'Y'])
        D <- as.vector(dem$d)
        SS <- sum(D^2)
        E <- diag(SS,ns,ns)/(ns-2)
        if (exterr){
            J <- matrix(0,ns,7)
            El <- getEl()
            colnames(J) <- colnames(El)
            McL <- mclean(tt=tt,d=x$d,exterr=exterr)
            if (type==1){
                J[,'U238'] <- -dem$dddb*a*McL$dPb206U238dl38
                J[,'U234'] <- -dem$dddb*a*McL$dPb206U238dl34
                J[,'Th230'] <- -dem$dddb*a*McL$dPb206U238dl30
                J[,'Ra226'] <- -dem$dddb*a*McL$dPb206U238dl26
            } else if (type==2){
                J[,'U235'] <- -dem$dddb*a*McL$dPb207U235dl35
                J[,'Pa231'] <- -dem$dddb*a*McL$dPb207U235dl31
            } else {
                J[,'Th232'] <- -dem$dddb*a*McL$dPb208Th232dl32
            }
            E <- E + J %*% El %*% t(J)
        }
    } else {
        if ('w'%in%names(ta0b0w) && !is.na(ta0b0w['w'])){
            # convert w from age to slope:
            ww <- ta0b0w['w']*abs(dbdt)
            DE <- york2DE(XY=O,a=a,b=b,w=ww)
        } else {
            DE <- york2DE(XY=O,a=a,b=b)
        }
        D <- DE$D
        if (exterr){
            ix <- 1:ns
            iy <- (ns+1):(2*ns)
            El <- getEl()
            J <- matrix(0,2*ns,7)
            colnames(J) <- colnames(El)
            P <- get.york.xy(XY=O,a=a,b=b)
            if (type==1){
                J[iy,'U238'] <- a*McL$dPb206U238dl38*P[,'x']
                J[iy,'U234'] <- a*McL$dPb206U238dl34*P[,'x']
                J[iy,'Th230'] <- a*McL$dPb206U238dl30*P[,'x']
                J[iy,'Ra226'] <- a*McL$dPb206U238dl26*P[,'x']
            } else if (type==2){
                J[iy,'U235'] <- a*McL$dPb207U235dl35*P[,'x']
                J[iy,'Pa231'] <- a*McL$dPb207U235dl31*P[,'x']
            } else {
                J[iy,'Th232'] <- a*McL$dPb208Th232dl32*P[,'x']
            }
            E <- DE$E + J %*% El %*% t(J)
        } else {
            E <- DE$E
        }
        if (!LL) SS <- stats::mahalanobis(D,center=FALSE,cov=E)
    }
    if (LL) out <- LL.norm(D,E)
    else out <- SS
    out
}

mswd.lud <- function(fit,x,exterr=FALSE,type='joint'){
    out <- fit
    ns <- length(x)
    np <- sum(diag(fit$cov)>0)
    if (x$format<4 || type%ni%c('joint',0)){
        out$df <- ns-np
    } else {
        out$df <- 2*ns-np
    }
    X <- x
    if (x$d$U48$option==2) X$d$U48 <- list(x=unname(fit$par['U48i']),option=1)
    if (x$d$ThU$option==2) X$d$ThU <- list(x=unname(fit$par['ThUi']),option=1)
    if (type%in%c('joint',0)){
        SS <- data2ludwig(X,ta0b0w=fit$par)$SS
    } else {
        SS <- LL.ludwig.2d(ta0b0w=fit$par,x=X,type=type,exterr=exterr,LL=FALSE)
    }
    if (out$df>0){
        out$mswd <- as.vector(SS/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    } else {
        out$mswd <- 1
        out$p.value <- 1
    }
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
