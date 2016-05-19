.IsoplotR <- new.env(parent = emptyenv())

#' Load settings to and from json
#'
#' Get and set preferred values for decay constants and isotopic
#' abundances from and to a \code{.json} file format
#'
#' @param fname the path of a \code{.json} file
#' @return if fname==NULL, returns a \code{.json} string
#' @examples
#' json <- system.file("defaults.json",package="IsoplotR")
#' settings(json)
#' print(settings())
#' @export
settings <- function(fname=NULL){
    if (is.null(fname)){
        preferences <- as.list(.IsoplotR)
        return(toJSON(preferences))
    } else {
        prefs <- fromJSON(file=fname)
        .IsoplotR$lambda <- prefs$lambda
        .IsoplotR$U238U235 <- prefs$U238U235
        .IsoplotR$I.A <- prefs$I.A
    }
}

#' Decay constants
#'
#' Gets or sets the decay constants of radioactive istopes
#'
#' @param nuclide the nuclide name
#' @param x new value for the decay constant
#' @param e new value for the decay constant uncertainty
#' @return if x == e == NULL, returns a two-item list containing:
#'
#' \code{x}: the decay constant  [in Ma-1]
#'
#' \code{e}: the standard error of the decay constant [in Ma-1]
#'
#' @examples
#' print(lambda('U238')$x)
#' # use the decay constant of Kovarik and Adams (1932)
#' lambda('U238',0.0001537,0.0000068)
#' print(lambda('U238')$x)
#' @export
lambda <- function(nuclide,x=NULL,e=NULL){
    if (is.null(x) & is.null(e)){
        return(.IsoplotR$lambda[[nuclide]])
    }
    if (is.numeric(x)){
        .IsoplotR$lambda[[nuclide]]$x <- x
    }
    if (is.numeric(e)){
        .IsoplotR$lambda[[nuclide]]$e <- e
    }
}

#' 238U/235 ratio
#'
#' Gets or sets the natural 238U/235 ratio. The default value of
#' 137.818 is taken from Hiess et al. (2012)
#'
#' @param x new value for 238U/235U ratio
#' @param e new value for its standard error
#' @return if x == e == NULL, returns a two-item list containing:
#'
#' \code{x}: the 238U/235U ratio
#'
#' \code{e}: the standard error of the 238U/235U ratio
#'
#' @examples
#' print(U238U235()$x)
#' # use the 238U/235U ratio of Steiger and Jaeger (1977)
#' U238U235(138.88,0)
#' print(U238U235()$x)
#' @export
U238U235 <- function(x=NULL,e=NULL){
    if (is.null(x) & is.null(e)){
        return(.IsoplotR$U238U235)
    }
    if (is.numeric(x)){
        .IsoplotR$U238U235$x <- x
    }
    if (is.numeric(e)){
        .IsoplotR$U238U235$e <- e
    }
}

#' Isotope abundance
#'
#' Gets or sets the natural abundance of isotopes
#'
#' @param nuclide one of either \code{'U'}, \code{'U238'},
#'     \code{'U235'}, or \code{'Th232'}
#' @param x new value for the isotope abundance
#' @param e new value for the standard error of the abundance
#' @return if x == e == NULL, returns a two element list containing:
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
#' print(I.A('U238')$x)
#' # use the 238U/235U ratio of Steiger and Jaeger (1977)
#' U238U235(138.88,0)
#' print(I.A('U238')$x)
#' @export
I.A <- function(nuclide,x=NULL,e=NULL){
    if (is.null(x) & is.null(e)){
        return(get.I.A(nuclide))        
    } else {
        set.I.A(nuclide,x,e)
    }
}

# helper function for I.A
get.I.A <- function(nuclide){
    out <- list()
    if (nuclide == 'U'){
        R <- U238U235()
        A238 <- R$x/(1+R$x)
        A235 <- 1/(1+R$x)
        out$x <- c(A238,A235)
        J <- matrix(c(1/(R$x+1)^2,-1/(R$x+1)^2),nrow=2,ncol=1)
        E <- R$e^2
        out$cov <- J %*% E %*% t(J)
        names(out$x) <- c('U238','U235')
        rownames(out$cov) <- names(out$x)
        colnames(out$cov) <- names(out$x)
    } else if (nuclide == 'U238'){
        out$x <- I.A('U')$x[1]
        out$e <- sqrt(I.A('U')$cov[1,1])
    } else if (nuclide == 'U235'){
        out$x <- I.A('U')$x[2]
        out$e <- sqrt(I.A('U')$cov[2,2])
    } else {
        out <- .IsoplotR$I.A[[nuclide]]
    }
    out
}

# helper function for I.A
set.I.A <- function(nuclide,x=NULL,e=NULL){
    if (is.numeric(x)){
        .IsoplotR$I.A[[nuclide]]$x <- x
    }
    if (is.numeric(e)){
        .IsoplotR$I.A[[nuclide]]$e <- e
    }
}
