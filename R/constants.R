#' Decay constants
#'
#' Returns the decay constants of radioactive istopes
#'
#' @param nuclide the nuclide name
#' @return a two-item list containing:
#'
#' \code{x}: the decay constant  [in Ma-1]
#'
#' \code{e}: the standard error of the decay constant
#' @examples
#' lambda('U238')
#' @export
lambda <- function(nuclide){
    out <- list()
    if (nuclide == 'U238'){
        out$x <- 1.55125e-4
        out$e <- out$x*0.0008
    }
    if (nuclide == 'U235'){
        out$x <- 9.8485e-4
        out$e <- out$x*0.0010
    }
    if (nuclide == 'Th232'){
        out$x <- 4.9475e-5
        out$e <- 0
    }
    out
}

R238235 <- function(){
    out <- list()
    out$x <- 137.818
    out$e <- 0.0225
    out
}

#' Isotope abundance
#'
#' Returns the natural abundance of isotopes
#'
#' @param nuclide one of either \code{'U'}, \code{'U238'},
#'     \code{'U235'}, or \code{'Th232'}
#' @return a two element list containing:
#'
#' \code{x}: a number or a vector of numbers between 0 (absent) and 1 (dominant)
#'
#' and
#'
#' \code{e}: the standard error or covariance matrix of x
#'
#' or, if \code{nuclide} = \code{'U'}:
#'
#' \code{cov}: the covariance matrix of all naturally occurring
#' isotopes
#'
#' @examples
#' I.A('U238')
#' @export
I.A <- function(nuclide){
    out <- list()
    if (nuclide == 'U'){
        R <- R238235()
        A238 <- R$x/(1+R$x)
        A235 <- 1/(1+R$x)
        out$x <- c(A238,A235)
        J <- matrix(c(1/(R$x+1)^2,-1/(R$x+1)^2),nrow=2,ncol=1)
        E <- R$e^2
        out$cov <- J %*% E %*% t(J)
        names(out$x) <- c('U238','U235')
        rownames(out$cov) <- names(out$x)
        colnames(out$cov) <- names(out$x)
    }
    if (nuclide == 'U238'){
        out$x <- I.A('U')$x[1]
        out$e <- sqrt(I.A('U')$cov[1,1])
    }
    if (nuclide == 'U235'){
        out$x <- I.A('U')$x[2]
        out$e <- sqrt(I.A('U')$cov[2,2])
    }
    if (nuclide == 'Th232'){
        out$x <- 1
        out$e <- 0
    }
    out
}
