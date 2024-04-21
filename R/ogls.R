#' @title
#' Omnivariant Generalised Least-Squares Regression
#'
#' @description
#' Linear regression with inter-sample error correlations.
#'
#' @param x either a \code{n x (n+1)} matrix obtained by prepending a
#'     vector of alternating \code{X,Y}-values to its covariance
#'     matrix OR an object of class \code{other} with
#'     \code{x$format=6}.
#' @param random.effects logical. If \code{TRUE}, quantifies the
#'     overdispersion associated with the y-intercept of the data.
#' @param ... optional arguments
#' @return a list of the slope and intercept of the best fit line as
#'     well as their standard errors and covariance.
#' @examples
#' fn <- system.file('UW137.csv',package='IsoplotR')
#' UW137 <- read.data(fn,method='other',format=6)
#' fit <- ogls(UW137)
#' @author Pieter Vermeesch and Mathieu Daëron
#' 
#' @references Daëron, M., 2023. Making the Case for Reconciled
#'     \eqn{\Delta}47 Calibrations Using Omnivariant Generalized
#'     Least-Squares Regression (No. EGU23-10066). Copernicus
#'     Meetings.
#'
#' Daëron & Vermeesch, in prep. Omnivariant Generalized Least Squares
#' Regression: Theory, Geochronological Applications, and Making the
#' Case for Reconciled \eqn{\Delta}47 calibrations, Chemical Geology.
#' 
#' @rdname ogls
#' @export
ogls <- function(x,...){ UseMethod("ogls",x) }
#' @rdname ogls
#' @export
ogls.default <- function(x,random.effects=FALSE,...){
    out <- list()
    ydat <- data2york(x=x,format=6)
    yfit <- york(ydat)
    ns <- nrow(x)/2
    if (random.effects){
        init <- c(yfit$a[1],yfit$b[1],lw=unname(log(yfit$a[2])))
        out <- stats::optim(init,LLw_ogls,dat=x,hessian=TRUE)
        SS <- LLw_ogls(out$par,dat=x,LL=FALSE)
        out$df <- ns-3
    } else {
        init <- c(yfit$a[1],yfit$b[1])
        omega <- solve(x[,-1])
        out <- stats::optim(init,LL_ogls,dat=x,omega=omega,hessian=TRUE)
        SS <- LL_ogls(out$par,dat=x,omega=omega,LL=FALSE)
        out$df <- ns-2
    }
    covmat <- solve(out$hessian)
    out$a <- c('a'=unname(out$par[1]),'s[a]'=unname(sqrt(covmat[1,1])))
    out$b <- c('b'=unname(out$par[2]),'s[b]'=unname(sqrt(covmat[2,2])))
    if (random.effects){
        w <- exp(out$par[3])
        sw <- w*sqrt(covmat[3,3])
        out$w <- c('w'=unname(w),'s[w]'=unname(sw))
    }
    out$cov.ab <- covmat[1,2]
    out$type <- "ogls"
    out <- append(out,getMSWD(X2=SS,df=out$df))
    out
}
#' @rdname ogls
#' @export
ogls.other <- function(x,random.effects=FALSE,...){
    ogls(x=x$x,random.effects=random.effects)
}

LL_ogls <- function(ab,dat,omega,LL=TRUE){
    a <- ab[1]
    b <- ab[2]
    ns <- nrow(dat)/2
    iX <- 1:ns
    iY <- (ns+1):(2*ns)
    X <- dat[iX,1,drop=FALSE]
    Y <- dat[iY,1,drop=FALSE]
    rY <- Y - a - b * X
    AA <- omega[iX,iX] + b*omega[iX,iY] + b*omega[iY,iX] + omega[iY,iY]*b^2
    BB <- (omega[iX,iY]+b*omega[iY,iY])%*%rY
    rx <- - solve(AA,BB)
    x <- X - rx
    ry <- dat[iY,1] - a - b*x
    v <- matrix(c(rx,ry),nrow=1,ncol=2*ns)
    SS <- (v %*% omega %*% t(v))
    if (LL){
        out <- (log(2*pi) - determinant(omega,logarithmic=TRUE)$modulus + SS)/2
    } else {
        out <- SS
    }
    out
}
LLw_ogls <- function(ablw,dat,LL=TRUE){
    w <- exp(ablw[3])
    ns <- nrow(dat)/2
    iX <- 1:ns
    iY <- (ns+1):(2*ns)
    E <- dat[,-1]
    diag(E[iY,iY]) <- diag(E[iY,iY]) + w^2
    LL_ogls(ab=ablw[-3],dat=dat,omega=solve(E),LL=LL)
}
