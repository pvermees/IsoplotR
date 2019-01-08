#' Linear regression of X,Y,Z-variables with correlated errors
#'
#' Implements the maximum likelihood algorithm of Ludwig and
#' Titterington (1994) for linear regression of three dimensional data
#' with correlated uncertainties.
#'
#' @details
#' Ludwig and Titterington (1994)'s 3-dimensional linear regression
#' algorithm for data with correlated uncertainties is an extension of
#' the 2-dimensional algorithm by Titterington and Halliday (1979),
#' which itself is equivalent to the algorithm of York et al. (2004).
#' Given \eqn{n} triplets of (approximately) collinear measurements
#' \eqn{X_i}, \eqn{Y_i} and \eqn{Z_i} (for \eqn{1 \leq i \leq n}),
#' their uncertainties \eqn{s[X_i]}, \eqn{s[Y_i]} and \eqn{s[Z_i]},
#' and their covariances cov[\eqn{X_i,Y_i}], cov[\eqn{X_i,Z_i}] and
#' cov[\eqn{Y_i,Z_i}], the \code{titterington} function fits two
#' slopes and intercepts with their uncertainties. It computes the
#' MSWD as a measure of under/overdispersion.  Overdispersed datasets
#' (MSWD>1) can be dealt with in the same three ways that are
#' described in the documentation of the \code{\link{isochron}}
#' function.
#'
#' @param x an \code{[n x 9]} matrix with the following columns:
#'     \code{X, sX, Y, sY, Z, sZ}, \code{rhoXY, rhoXZ, rhoYZ}.
#' @param alpha cutoff value for confidence intervals
#' @return A four-element list of vectors containing:
#'
#'     \describe{
#'     \item{par}{4-element vector \code{c(a,b,A,B)} where \code{a} is
#'               the intercept of the \code{X-Y} regression, \code{b}
#'               is the slope of the \code{X-Y} regression, \code{A}
#'               is the intercept of the \code{X-Z} regression, and
#'               \code{B} is the slope of the \code{X-Z} regression.}
#'
#'     \item{cov}{\code{[4 x 4]}-element covariance matrix of \code{par}}
#'
#'     \item{mswd}{the mean square of the residuals (a.k.a `reduced
#'                 Chi-square') statistic}
#'
#'     \item{p.value}{p-value of a Chi-square test for linearity}
#'
#'     \item{df}{the number of degrees of freedom for the
#'               Chi-square test (3\eqn{n}-3)}
#'
#'     \item{tfact}{the \eqn{100(1-\alpha/2)\%} percentile of the
#'                  t-distribution with \eqn{(n-2k+1)} degrees of freedom}
#' }
#' @examples
#' d <- matrix(c(0.1677,0.0047,1.105,0.014,0.782,0.015,0.24,0.51,0.33,
#'               0.2820,0.0064,1.081,0.013,0.798,0.015,0.26,0.63,0.32,
#'               0.3699,0.0076,1.038,0.011,0.819,0.015,0.27,0.69,0.30,
#'               0.4473,0.0087,1.051,0.011,0.812,0.015,0.27,0.73,0.30,
#'               0.5065,0.0095,1.049,0.010,0.842,0.015,0.27,0.76,0.29,
#'               0.5520,0.0100,1.039,0.010,0.862,0.015,0.27,0.78,0.28),
#'             nrow=6,ncol=9)
#' colnames(d) <- c('X','sX','Y','sY','Z','sZ','rXY','rXZ','rYZ')
#' titterington(d)
#' @seealso \code{\link{york}}, \code{\link{isochron}}, \code{\link{ludwig}}
#' @references
#' Ludwig, K.R. and Titterington, D.M., 1994. Calculation
#' of \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#' Cosmochimica Acta, 58(22), pp.5031-5042.
#'
#' Titterington, D.M. and Halliday, A.N., 1979. On the fitting of
#' parallel isochrons and the method of maximum likelihood. Chemical
#' Geology, 26(3), pp.183-195.
#'
#' York, D., Evensen, N.M., Martinez, M.L. and De Basebe Delgado, J., 2004.
#' Unified equations for the slope, intercept, and standard
#' errors of the best straight line. American Journal of Physics,
#' 72(3), pp.367-375.
#' @export
titterington <- function(x,alpha=0.05){
    ns <- nrow(x)
    fitXY <- york(subset(x,select=c(1,2,3,4,7)))
    a <- fitXY$a[1]
    b <- fitXY$b[1]
    fitXZ <- york(subset(x,select=c(1,2,5,6,8)))
    A <- fitXZ$a[1]
    B <- fitXZ$b[1]
    init <- c(a,b,A,B)
    dat <- matrix2covlist(x)
    fit <- stats::optim(init,fn=S.tit,gr=gr.tit,method="BFGS",dat)
    fish <- fisher.tit(fit$par,dat)
    covmat <- solve(fish)
    out <- list()
    out$par <- fit$par
    out$cov <- covmat[(ns+1):(ns+4),(ns+1):(ns+4)]
    mswd <- mswd.tit(fit$par,dat)
    out <- c(out,mswd)
    parnames <- c('a','b','A','B')
    names(out$par) <- parnames
    rownames(out$cov) <- parnames
    colnames(out$cov) <- parnames
    out$type <- 'titterington'
    out
}

mswd.tit <- function(abAB,dat){
    a <- abAB[1]
    b <- abAB[2]
    A <- abAB[3]
    B <- abAB[4]
    XYZ <- dat$XYZ
    omega <- dat$omega
    ns <- length(XYZ)
    S <- 0
    for (i in 1:ns){
        X <- XYZ[[i]][1]
        Y <- XYZ[[i]][2]
        Z <- XYZ[[i]][3]
        abg <- alpha.beta.gamma(a,b,A,B,XYZ[[i]],omega[[i]])
        alpha <- abg[1]
        beta <- abg[2]
        gamma <- abg[3]
        x <- X + beta/alpha
        S <- S + alpha*(X-x)^2 + 2*beta*(X-x) + gamma
    }
    out <- list()
    out$df <- 2*ns-4
    out$mswd <- S/out$df
    out$p.value <- as.numeric(1-stats::pchisq(S,out$df))
    out
}

matrix2covlist <- function(x){
    ns <- nrow(x)
    out <- list()
    out$XYZ <- list()
    out$omega <- list()
    covmat <- matrix(0,3,3)
    for (i in 1:ns){
        covmat <- cor2cov3(x[i,2],x[i,4],x[i,6],x[i,7],x[i,8],x[i,9])
        out$XYZ[[i]] <- x[i,c(1,3,5)]
        out$omega[[i]] <- solve(covmat)
    }
    out
}

alpha.beta.gamma <- function(a,b,A,B,XYZ,omega){
    r <- XYZ[2]-a-b*XYZ[1]
    R <- XYZ[3]-A-B*XYZ[1]
    alpha <- omega[1,1] + omega[2,2]*b^2 + omega[3,3]*B^2 +
        2*omega[1,2]*b + 2*omega[1,3]*B + 2*omega[2,3]*b*B
    beta <- r*(omega[1,2] + b*omega[2,2] + B*omega[2,3]) +
        R*(omega[1,3] + b*omega[2,3] + B*omega[3,3])
    gamma <- omega[2,2]*r^2 + omega[3,3]*R^2 + 2*omega[2,3]*r*R
    c(alpha,beta,gamma)
}

fisher.tit <- function(abAB,dat){
    a <- abAB[1] # y intercept
    b <- abAB[2] # y slope
    A <- abAB[3] # z intercept
    B <- abAB[4] # z slope
    ns <- length(dat$XYZ)
    out <- matrix(0,ns+4,ns+4)
    d2L.da2 <- 0
    d2L.dadb <- 0
    d2L.db2 <- 0
    d2L.dA2 <- 0
    d2L.dAdB <- 0
    d2L.dB2 <- 0
    d2L.dadA <- 0
    d2L.dadB <- 0
    d2L.dbdA <- 0
    d2L.dbdB <- 0
    for (i in 1:ns){
        XYZ <- dat$XYZ[[i]]
        X <- XYZ[1]
        Y <- XYZ[2]
        Z <- XYZ[3]
        omega <- dat$omega[[i]]
        abg <- alpha.beta.gamma(a,b,A,B,XYZ,omega)
        alpha <- abg[1]
        beta <- abg[2]
        gamma <- abg[3]
        x <- X + beta/alpha
        d2L.da2 <- d2L.da2 - omega[2,2]
        d2L.dadb <- d2L.dadb - x*omega[2,2]
        d2L.db2 <- d2L.db2 - omega[2,2]*x^2
        d2L.dA2 <- d2L.dA2 - omega[3,3]
        d2L.dAdB <- d2L.dAdB - x*omega[3,3]
        d2L.dB2 <- d2L.dB2 - omega[3,3]*x^2
        d2L.dadA <- d2L.dadA - omega[2,3]
        d2L.dadB <- d2L.dadB - x*omega[2,3]
        d2L.dbdA <- d2L.dbdA - x*omega[2,3]
        d2L.dbdB <- d2L.dbdB - omega[2,3]*x^2
        
        d2L.dadx <- -(omega[1,2] + b*omega[2,2] + B*omega[2,3])
        d2L.dbdx <- -x*(omega[1,2] + b*omega[2,2] + B*omega[2,3]) +
            omega[2,2]*(Y-a-b*x) + omega[1,2]*(X-x) + omega[2,3]*(Z-A-B*x)
        d2L.dAdx <- -(omega[1,3] + b*omega[2,3] + B*omega[3,3])
        d2L.dBdx <- -x*(omega[1,3] + b*omega[2,3] + B*omega[3,3]) +
            omega[2,3]*(Y-a-b*x) + omega[1,3]*(X-x) + omega[3,3]*(Z-A-B*x)
        d2L.dx2 <- -(omega[1,1] + 2*b*omega[1,2] + 2*B*omega[1,3] +
                     omega[2,2]*b^2 + omega[2,3]*B^2 + 2*omega[3,3]*b*B)
        out[i,i] <- -d2L.dx2
        out[i,ns+1] <- -d2L.dadx
        out[i,ns+2] <- -d2L.dbdx
        out[i,ns+3] <- -d2L.dAdx
        out[i,ns+4] <- -d2L.dBdx
        out[ns+1,i] <- -d2L.dadx
        out[ns+2,i] <- -d2L.dbdx
        out[ns+3,i] <- -d2L.dAdx
        out[ns+4,i] <- -d2L.dBdx
    }
    out[ns+1,ns+1] <- -d2L.da2
    out[ns+2,ns+2] <- -d2L.db2
    out[ns+3,ns+3] <- -d2L.dA2
    out[ns+4,ns+4] <- -d2L.dB2
    out[ns+1,ns+2] <- -d2L.dadb
    out[ns+1,ns+3] <- -d2L.dadA
    out[ns+1,ns+4] <- -d2L.dadB
    out[ns+2,ns+1] <- -d2L.dadb
    out[ns+2,ns+3] <- -d2L.dbdA
    out[ns+2,ns+4] <- -d2L.dbdB
    out[ns+3,ns+1] <- -d2L.dadA
    out[ns+3,ns+2] <- -d2L.dbdA
    out[ns+3,ns+4] <- -d2L.dAdB
    out[ns+4,ns+1] <- -d2L.dadB
    out[ns+4,ns+2] <- -d2L.dbdB
    out[ns+4,ns+3] <- -d2L.dAdB
    out
}

gr.tit <- function(abAB,dat){
    ns <- length(dat)
    a <- abAB[1] # y intercept
    b <- abAB[2] # y slope
    A <- abAB[3] # z intercept
    B <- abAB[4] # z slope
    dS.da <- 0
    dS.db <- 0
    dS.dA <- 0
    dS.dB <- 0
    dS.dS <- -1/2
    dr.da <- -1
    dR.dA <- -1
    dR.da <- 0
    dr.dA <- 0
    dR.db <- 0
    dr.dB <- 0
    dalpha.da <- 0
    dalpha.dA <- 0
    for (i in 1:ns){
        XYZ <- dat$XYZ[[i]]
        X <- XYZ[1]
        Y <- XYZ[2]
        Z <- XYZ[3]
        omega <- dat$omega[[i]]
        r <- Y-a-b*X
        R <- Z-A-B*X
        abg <- alpha.beta.gamma(a,b,A,B,XYZ,omega)
        alpha <- abg[1]
        beta <- abg[2]
        gamma <- abg[3]
        dr.db <- -X
        dR.dB <- -X
        dbeta.da <- dr.da*(omega[1,2] + b*omega[2,2] + B*omega[2,3]) +
            dR.da*(omega[1,3] + b*omega[2,3] + B*omega[3,3])
        dgamma.da <- 2*(omega[2,2]*r*dr.da + omega[3,3]*R*dR.da +
                        omega[2,3]*R*dr.da + omega[2,3]*r*dR.da)
        dalpha.db <- 2*(omega[2,2]*b + omega[1,2] + omega[2,3]*B)
        dbeta.db <- dr.db*(omega[1,2] + b*omega[2,2] + B*omega[2,3]) +
            dR.db*(omega[1,3] + b*omega[2,3] + B*omega[3,3]) +
            r*omega[2,2] + R*omega[2,2]
        dgamma.db <- 2*(omega[2,2]*r*dr.db + omega[3,3]*R*dR.db +
                        omega[2,3]*R*dr.db + omega[2,3]*r*dR.db)
        dbeta.dA <- dr.dA*(omega[1,2] + b*omega[2,2] + B*omega[2,3]) +
            dR.dA*(omega[1,3] + b*omega[2,3] + B*omega[3,3])
        dgamma.dA <- 2*(omega[2,2]*r*dr.dA + omega[3,3]*R*dR.dA +
                        omega[2,3]*R*dr.dA + omega[2,3]*r*dR.dA)
        dalpha.dB <- 2*(omega[3,3]*B + omega[1,3] + omega[2,3]*b)
        dbeta.dB <- dr.dB*(omega[1,2] + b*omega[2,2] + B*omega[2,3]) +
            dR.dB*(omega[1,3] + b*omega[2,3] + B*omega[3,3]) +
            r*omega[2,3] + R*omega[3,3]
        dgamma.dB <- 2*(omega[2,2]*r*dr.dB + omega[3,3]*R*dR.dB +
                        omega[2,3]*R*dr.dB + omega[2,3]*r*dR.dB)
        dS.da <- dgamma.da - 2*(beta/alpha)*dbeta.da + dalpha.da*(beta/alpha)^2
        dS.db <- dgamma.db - 2*(beta/alpha)*dbeta.db + dalpha.db*(beta/alpha)^2
        dS.dA <- dgamma.dA - 2*(beta/alpha)*dbeta.dA + dalpha.dA*(beta/alpha)^2
        dS.dB <- dgamma.dB - 2*(beta/alpha)*dbeta.dB + dalpha.dB*(beta/alpha)^2
    }
    c(dS.da,dS.db,dS.dA,dS.dB)
}

S.tit <- function(abAB,dat){
    ns <- length(dat)
    a <- abAB[1] # y intercept
    b <- abAB[2] # y slope
    A <- abAB[3] # z intercept
    B <- abAB[4] # z slope
    S <- 0       # initialise sum of squares
    for (i in 1:ns){
        abg <- alpha.beta.gamma(a,b,A,B,dat$XYZ[[i]],dat$omega[[i]])
        alpha <- abg[1]
        beta <- abg[2]
        gamma <- abg[3]
        S <- S + gamma - (beta^2)/alpha
    }
    S
}

# dat is the output of matrix2covlist
get.titterington.xy <- function(dat,abAB){
    ns <- length(dat$omega)
    out <- matrix(NA,ns,3)
    colnames(out) <- c('X','Y','Z')
    for (i in 1:ns){
        abg <- alpha.beta.gamma(abAB[1],abAB[2],abAB[3],abAB[4],
                                dat$XYZ[[i]],dat$omega[[i]])
        XYZ <- dat$XYZ[[i]]
        out[i,1] <- XYZ[1] + abg[2]/abg[1]
        out[i,2] <- XYZ[2] + abAB[1] + abAB[2]*XYZ[1]
        out[i,2] <- XYZ[3] + abAB[3] + abAB[4]*XYZ[1]
    }
    out
}

data2tit <- function(x,...){ UseMethod("data2tit",x) }
data2tit.default <- function(x,...){ stop('default function undefined') }
# osmond = FALSE: 8/2 - 4/2 - 0/2
# osmond = TRUE : 2/8 - 4/8 - 0/8
data2tit.ThU <- function(x,osmond=TRUE,generic=TRUE,...){
    ns <- length(x)
    out <- matrix(0,ns,9)
    if ((x$format == 1 & !osmond) | (x$format == 2 & osmond)){
        out <- x$x
    } else {
        for (i in 1:ns){
            out[i,] <- ThU.convert(x$x[i,])
        }
    }
    if (generic){
        colnames(out) <- c('X','sX','Y','sY','Z','sZ','rXY','rXZ','rYZ')
    } else if (osmond){
        colnames(out) <- c('Th232U238','sTh232U238','U234U238','sU234U238',
                           'Th230U238','sTh230U238','rXY','rXZ','rYZ')
    } else {
        colnames(out) <- c('U238Th232','sU238Th232','U234Th232','su234Th232',
                           'Th230Th232','sTh230Th232','rXY','rXZ','rYZ')
    }
    out
}

# ia = index for the chosen intercept
# ib = index for the chosen slope
tit2york <- function(fit,ia,ib){
    out <- list()
    out$a <- c(fit$par[ia],sqrt(fit$cov[ia,ia]))
    out$b <- c(fit$par[ib],sqrt(fit$cov[ib,ib]))
    out$cov.ab <- fit$cov[ia,ib]
    out$y0 <- c(0,0)
    out$p.value <- 0
    out$mswd <- fit$mswd
    out
}
