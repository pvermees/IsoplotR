.IsoplotR <- new.env(parent = emptyenv())

#' Load settings to and from json
#'
#' Get and set preferred values for decay constants and isotopic
#' abundances from and to a \code{.json} file format
#'
#' @param fname the path of a \code{.json} file
#' @return if \code{fname==NULL}, returns a \code{.json} string
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
        .IsoplotR$iratio <- prefs$iratio
    }
}

#' Decay constants
#'
#' Gets or sets the decay constants of radioactive isotopes
#'
#' @param nuclide the nuclide name
#' @param x new value for the decay constant
#' @param e new value for the decay constant uncertainty
#' @return if \code{x==e==NULL}, returns a two-item vector containing
#'     the decay constant [in Ma-1] and its standard error,
#'     respectively.
#' @examples print(lambda('U238'))
#' # use the decay constant of Kovarik and Adams (1932)
#' lambda('U238',0.0001537,0.0000068)
#' print(lambda('U238'))
#' @references
#' U: Jaffey, A. H., et al. "Precision measurement of half-lives and
#' specific activities of U235 and U238." Physical Review C 4.5 (1971): 1889.
#'
#' Th: Le Roux, L. J., and L. E. Glendenin. "Half-life of 232Th.
#' "Proceedings of the National Meeting on Nuclear Energy, Pretoria, South Africa. 1963.
#' 
#' Ar: Renne, Paul R., et al. "Response to the comment by WH Schwarz et al. on
#' "Joint determination of 40K decay constants and 40Ar*/40K for the Fish Canyon
#' sanidine standard, and improved accuracy for 40 Ar/39 Ar geochronology"
#' by PR Renne et al.(2010)." Geochimica et Cosmochimica Acta 75.17 (2011): 5097-5100.
#' @export
lambda <- function(nuclide,x=NULL,e=NULL){
    if (is.null(x) & is.null(e)) return(.IsoplotR$lambda[[nuclide]])
    if (is.numeric(x)) .IsoplotR$lambda[[nuclide]][1] <- x
    if (is.numeric(e)) .IsoplotR$lambda[[nuclide]][2] <- e
}

#' Isotopic ratios
#'
#' Gets or sets natural isotopic ratios.
#'
#' @param ratio one of either \code{'U238U235'}, \code{'Ar40Ar36'},
#'     \code{'Ar38Ar36'}, \code{'Rb85Rb87'}, \code{'Sr88Sr86'},
#'     \code{'Sr87Sr86'}, \code{'Sr84Sr86'}, \code{'Re185Re187'},
#'     \code{'Os184Os192'}, \code{'Os186Os192'}, \code{'Os187Os192'},
#'     \code{'Os188Os192'}, \code{'Os189Os192'}
#' @param x new value for ratio
#' @param e new value for its standard error
#' @return if \code{x=e=NULL}, returns a two-item vector containing
#'     the mean value of the requested ratio and its standard error,
#'     respectively.
#' @examples
#' # returns the 238U/235U ratio of Hiess et al. (2012):
#' print(iratio('U238U235'))
#' # use the 238U/235U ratio of Steiger and Jaeger (1977):
#' iratio('U238U235',138.88,0)
#' print(iratio('U238U235'))
#' @references
#' Ar: Lee, Jee-Yon, et al. "A redetermination of the isotopic abundances
#' of atmospheric Ar." Geochimica et Cosmochimica Acta 70.17 (2006): 4507-4512.
#'
#' Rb: Catanzaro, E. J., et al. "Absolute isotopic abundance ratio and
#' atomic weight of terrestrial rubidium." J. Res. Natl. Bur. Stand. A 73
#' (1969): 511-516.
#'
#' Sr: Moore, L. J., et al. "Absolute isotopic abundance ratios and atomic
#' weight of a reference sample of strontium." J. Res. Natl.Bur. Stand.
#' 87.1 (1982): 1-8.
#'
#' Re: Gramlich, John W., et al. "Absolute isotopic abundance ratio and
#' atomic weight of a reference sample of rhenium." J. Res. Natl. Bur.
#' Stand. A 77 (1973): 691-698.
#'
#' Os: Voelkening, Joachim, Thomas Walczyk, and Klaus G. Heumann.
#' "Osmium isotope ratio determinations by negative thermal ionization
#' mass spectrometry." Int. J. Mass Spect. Ion Proc. 105.2 (1991): 147-159.
#' 
#' U: Hiess, Joe, et al. "\eqn{^{238}}U/\eqn{^{235}}U systematics in terrestrial
#' uranium-bearing minerals." Science 335.6076 (2012): 1610-1614.
#' @export
iratio <- function(ratio,x=NULL,e=NULL){
    if (is.null(x) & is.null(e)) return(.IsoplotR$iratio[[ratio]])
    if (is.numeric(x)) .IsoplotR$iratio[[ratio]][1] <- x
    if (is.numeric(e)) .IsoplotR$iratio[[ratio]][2] <- e
}
