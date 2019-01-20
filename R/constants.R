.IsoplotR <- new.env(parent = emptyenv())

lambda <- function(nuclide,x=NULL,e=NULL){
    if (is.null(x) & is.null(e)) return(.IsoplotR$lambda[[nuclide]])
    if (is.numeric(x)) .IsoplotR$lambda[[nuclide]][1] <- x
    if (is.numeric(e)) .IsoplotR$lambda[[nuclide]][2] <- e
}
iratio <- function(ratio,x=NULL,e=NULL){
    if (is.null(x) & is.null(e)) return(.IsoplotR$iratio[[ratio]])
    if (is.numeric(x)) .IsoplotR$iratio[[ratio]][1] <- x
    if (is.numeric(e)) .IsoplotR$iratio[[ratio]][2] <- e
}
imass <- function(nuclide,x=NULL,e=NULL){
    if (is.null(x) & is.null(e)) return(.IsoplotR$imass[[nuclide]])
    if (is.numeric(x)) .IsoplotR$imass[[nuclide]][1] <- x
    if (is.numeric(e)) .IsoplotR$imass[[nuclide]][2] <- e
}
etchfact <- function(mineral,x=NULL){
    if (is.null(x)) return(.IsoplotR$etchfact[[mineral]])
    if (is.numeric(x)) .IsoplotR$etchfact[[mineral]] <- x
}
tracklength <- function(mineral,x=NULL){
    if (is.null(x)) return(.IsoplotR$tracklength[[mineral]])
    if (is.numeric(x)) .IsoplotR$tracklength[[mineral]] <- x
}
mindens <- function(mineral,x=NULL){
    if (is.null(x)) return(.IsoplotR$mindens[[mineral]])
    if (is.numeric(x)) .IsoplotR$mindens[[mineral]] <- x
}

#' Load settings to and from json
#'
#' Get and set preferred values for decay constants, isotopic
#' abundances, molar masses, fission track etch efficiences, and
#' etchable lengths, and mineral densities, either individually or via
#' a \code{.json} file format.
#'
#' @param setting unless \code{fname} is provided, this should be one
#'     of either:
#'
#' \code{'lambda'}: to get and set decay constants
#'
#' \code{'iratio'}: isotopic ratios
#'
#' \code{'imass'}: isotopic molar masses
#'
#' \code{'mindens'}: mineral densities
#'
#' \code{'etchfact'}: fission track etch efficiency factors
#'
#' \code{'tracklength'}: equivalent isotropic fission track length
#' @param ... depends on the value for \code{setting}:
#'
#' \itemize{
#'
#' \item for \code{'lambda'}: the isotope of interest (one of either
#' \code{"fission"}, \code{"U238"}, \code{"U235"}, \code{"U234"},
#' \code{"Th232"}, \code{"Th230"}, \code{"Re187"}, \code{"Sm147"},
#' \code{"Rb87"}, \code{"Lu176"}, or \code{"K40"}) PLUS (optionally)
#' the decay constant value and its analytical error.  Omitting the
#' latter two numbers simply returns the existing values.
#'
#' \item for \code{'iratio'}: the isotopic ratio of interest (one of
#' either \code{"Ar40Ar36"}, \code{"Ar38Ar36"},\code{"Ca40Ca44"},
#' \code{"Rb85Rb87"}, \code{"Sr88Sr86"}, \code{"Sr87Sr86"},
#' \code{"Sr84Sr86"}, \code{"Re185Re187"}, \code{"Os184Os192"}
#' \code{"Os186Os192"}, \code{"Os187Os192"}, \code{"Os188Os192"},
#' \code{"Os189Os192"}, \code{"Os190Os192"}, \code{"U238U235"},
#' \code{"Sm144Sm152"}, \code{"Sm147Sm152"}, \code{"Sm148Sm152"},
#' \code{"Sm149Sm152"}, \code{"Sm150Sm152"}, \code{"Sm154Sm152"},
#' \code{"Nd142Nd144"}, \code{"Nd143Nd144"}, \code{"Nd145Nd144"},
#' \code{"Nd146Nd144"}, \code{"Nd148Nd144"}, \code{"Nd150Nd144"},
#' \code{"Lu176Lu175"}, \code{"Hf174Hf177"}, \code{"Hf176Hf177"},
#' \code{"Hf178Hf177"}, \code{"Hf179Hf177"}, \code{"Hf180Hf177"}) PLUS
#' (optionally) the isotopic ratio and its analytical error.  Omitting
#' the latter two numbers simply returns the existing values.
#'
#' \item for \code{'imass'}: the (isotopic) molar mass of interest
#' (one of either \code{"U"}, \code{"Rb"}, \code{"Rb85"},
#' \code{"Rb87"}, \code{"Sr84"}, \code{"Sr86"}, \code{"Sr87"},
#' \code{"Sr88"}, \code{"Re"}, \code{"Re185"}, \code{"Re187"},
#' \code{"Os"}, \code{"Os184"}, \code{"Os186"}, \code{"Os187"},
#' \code{"Os188"}, \code{"Os189"}, \code{"Os190"}, \code{"Os192"},
#' \code{"Sm"}, \code{"Nd"}, \code{"Lu"}, \code{"Hf"}) PLUS
#' (optionally) the molar mass and its analytical error.  Omitting the
#' latter two numbers simply returns the existing values.
#'
#' \item for \code{'mindens'}: the mineral of interest (one of either
#' \code{"apatite"} or \code{"zircon"}) PLUS the mineral
#' density. Omitting the latter number simply returns the existing
#' value.
#'
#' \item \code{'etchfact'}: the mineral of interest (one of either
#' \code{"apatite"} or \code{"zircon"}) PLUS the etch efficiency
#' factor. Omitting this number simply returns the existing value.
#'
#' \item \code{'tracklength'}: the mineral of interest (one of either
#' \code{"apatite"} or \code{"zircon"}) PLUS the equivalent isotropic
#' fission track length. Omitting this number simply returns the
#' existing value.
#' }
#' @param fname the path of a \code{.json} file
#' @return if \code{setting=NA} and \code{fname=NA}, returns a
#'     \code{.json} string
#'
#' if \code{...} contains only the name of an isotope, isotopic ratio,
#' element, or mineral and no new value, then \code{settings} returns
#' either a scalar with the existing value, or a two-element vector
#' with the value and its uncertainty.
#'
#' @references
#' \enumerate{
#' \item Decay constants:
#' \itemize{
#'
#' \item \eqn{^{238}}U, \eqn{^{235}}U: Jaffey, A. H., et
#' al. "Precision measurement of half-lives and specific
#' activities of U\eqn{^{235}} and U\eqn{^{238}}."
#' Physical Review C 4.5 (1971): 1889.
#'
#' \item \eqn{^{232}}Th: Le Roux, L. J., and
#' L. E. Glendenin. "Half-life of \eqn{^{232}}Th.", Proceedings of the
#' National Meeting on Nuclear Energy, Pretoria, South Africa. 1963.
#'
#' \item \eqn{^{234}}U, \eqn{^{230}}Th: Cheng, H., Edwards, R.L.,
#' Shen, C.C., Polyak, V.J., Asmerom, Y., Woodhead, J., Hellstrom, J.,
#' Wang, Y., Kong, X., Spotl, C. and Wang, X., 2013. Improvements in
#' \eqn{^{230}}Th dating, \eqn{^{230}}Th and \eqn{^{234}}U half-life
#' values, and U-Th isotopic measurements by multi-collector
#' inductively coupled plasma mass spectrometry. Earth and Planetary
#' Science Letters, 371, pp.82-91.
#'
#' \item \eqn{^{231}}Pa: Audi, G., Bersillon, O., Blachot, J. and
#' Wapstra, A.H., 2003. The NUBASE evaluation of nuclear and decay
#' properties. Nuclear Physics A, 729(1), pp.3-128.
#'
#' \item \eqn{^{226}}Ra: Audi, G., Bersillon, O., Blachot, J. and
#' Wapstra, A.H., 2003. The NUBASE evaluation of nuclear and decay
#' properties. Nuclear Physics A, 729(1), pp.3-128.
#'
#' \item Sm: Lugmair, G. W., and
#' K. Marti. "Lunar initial \eqn{^{143}}Nd/\eqn{^{144}}Nd: differential
#' evolution of the lunar crust and mantle."
#' Earth and Planetary Science Letters 39.3 (1978): 349-357.
#'
#' \item Nd: Zhao, Motian, et
#' al. "Absolute measurements of neodymium isotopic abundances and atomic
#' weight by MC-ICPMS." International Journal of Mass Spectrometry 245.1
#' (2005): 36-40.
#'
#' \item Re: Selby, D., Creaser, R.A., Stein, H.J., Markey, R.J. and Hannah,
#' J.L., 2007.  Assessment of the 187Re decay constant by cross
#' calibration of Re-Os molybdenite and U-Pb zircon chronometers in
#' magmatic ore systems. Geochimica et Cosmochimica Acta, 71(8),
#' pp.1999-2013.
#'
#' \item Ar: Renne, Paul R., et
#' al. "Response to the comment by WH Schwarz et al. on "Joint
#' determination of \eqn{^{40}}K decay constants and
#' \eqn{^{40}}Ar*/\eqn{^{40}}K for the Fish Canyon sanidine standard,
#' and improved accuracy for \eqn{^{40}}Ar/\eqn{^{39}}Ar
#' geochronology" by PR Renne et al.(2010)." Geochimica et
#' Cosmochimica Acta 75.17 (2011): 5097-5100.
#'
#' \item Rb: Villa, I.M., De Bievre, P., Holden, N.E. and Renne, P.R., 2015.
#' "IUPAC-IUGS recommendation on the half life of \eqn{^{87}}Rb".
#' Geochimica et Cosmochimica Acta, 164, pp.382-385.
#'
#' \item Lu: Soederlund, Ulf, et al. "The \eqn{^{176}}Lu decay constant
#' determined by Lu-Hf and U-Pb isotope systematics of Precambrian
#' mafic intrusions." Earth and Planetary Science Letters 219.3 (2004): 311-324.
#'
#' }
#' \item Isotopic ratios:
#' \itemize{
#' \item Ar: Lee, Jee-Yon, et
#' al. "A redetermination of the isotopic abundances of atmospheric Ar."
#' Geochimica et Cosmochimica Acta 70.17 (2006): 4507-4512.
#'
#' \item Ca: Moore, L.J. and Machlan, L.A., 1972. High-accuracy
#' determination of calcium in blood serum by isotope dilution mass
#' spectrometry. Analytical chemistry, 44(14), pp.2291-2296.
#'
#' \item Rb: Catanzaro, E. J., et al. "Absolute isotopic abundance ratio and
#' atomic weight of terrestrial rubidium." J. Res. Natl. Bur. Stand. A 73
#' (1969): 511-516.
#'
#' \item Sr: Moore, L. J., et al. "Absolute isotopic abundance ratios and atomic
#' weight of a reference sample of strontium." J. Res. Natl. Bur. Stand.
#' 87.1 (1982): 1-8.
#'
#' and (for \eqn{^{87}}Sr\eqn{^{86}}Sr):
#'
#' Compston, W., Berry, H., Vernon, M.J., Chappell, B.W. and Kaye,
#' M.J., 1971. Rubidium-strontium chronology and chemistry of lunar
#' material from the Ocean of Storms. In Lunar and Planetary Science
#' Conference Proceedings (Vol. 2, p. 1471).
#'
#' \item Sm: Chang, Tsing-Lien, et al. "Absolute isotopic composition and atomic
#' weight of samarium." International Journal of Mass Spectrometry 218.2
#' (2002): 167-172.
#'
#' \item Re: Gramlich, John W., et al. "Absolute isotopic abundance ratio and
#' atomic weight of a reference sample of rhenium." J. Res. Natl. Bur.
#' Stand. A 77 (1973): 691-698.
#'
#' \item Os: Voelkening, Joachim, Thomas Walczyk, and Klaus G. Heumann.
#' "Osmium isotope ratio determinations by negative thermal ionization
#' mass spectrometry." Int. J. Mass Spect. Ion Proc. 105.2 (1991): 147-159.
#'
#' \item Lu: De Laeter, J. R., and N. Bukilic.
#' "Solar abundance of \eqn{^{176}}Lu and s-process nucleosynthesis."
#' Physical Review C 73.4 (2006): 045806.
#'
#' \item Hf: Patchett, P. Jonathan.
#' "Importance of the Lu-Hf isotopic system in studies of planetary
#' chronology and chemical evolution." Geochimica et Cosmochimica
#' Acta 47.1 (1983): 81-91.
#'
#' \item U: Hiess, Joe, et al. "\eqn{^{238}}U/\eqn{^{235}}U systematics in terrestrial
#' uranium-bearing minerals." Science 335.6076 (2012): 1610-1614.
#' }
#' }
#' @seealso \code{\link{read.data}}
#' @examples
#' # load and show the default constants that come with IsoplotR
#' json <- system.file("constants.json",package="IsoplotR")
#' settings(fname=json)
#' print(settings())
#'
#' # use the decay constant of Kovarik and Adams (1932)
#' settings('lambda','U238',0.0001537,0.0000068)
#' print(settings('lambda','U238'))
#'
#' # returns the 238U/235U ratio of Hiess et al. (2012):
#' print(settings('iratio','U238U235'))
#' # use the 238U/235U ratio of Steiger and Jaeger (1977):
#' settings('iratio','U238U235',138.88,0)
#' print(settings('iratio','U238U235'))
#' @export
settings <- function(setting=NA,...,fname=NA){
    args <- list(...)
    nargs <- length(args)
    if (!is.na(fname)){
        prefs <- fromJSON(file=fname)
        .IsoplotR$lambda <- prefs$lambda
        .IsoplotR$iratio <- prefs$iratio
        .IsoplotR$imass <- prefs$imass
        .IsoplotR$etchfact <- prefs$etchfact
        .IsoplotR$tracklength <- prefs$tracklength
        .IsoplotR$mindens <- prefs$mindens
    } else if (!is.na(setting) & (nargs>0)){
        if (nargs==1){
            Rcommand <- paste0(setting,"('",args[[1]],"')")
            return(eval(parse(text=Rcommand)))
        } else if (nargs==2) {
            Rcommand <- paste0(setting,"('",args[[1]],"',",args[[2]],")")
        } else if (nargs==3) {
            Rcommand <- paste0(setting,"('",args[[1]],"',",args[[2]],",",args[[3]],")")
        }
        eval(parse(text=Rcommand))
    } else {
        preferences <- as.list(.IsoplotR)
        return(toJSON(preferences))
    }
}
