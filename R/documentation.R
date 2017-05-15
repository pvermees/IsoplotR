#' Example datasets for testing \code{IsoplotR}
#'
#' U-Pb, Ar-Ar, Re-Os, Sm-Nd, Rb-Sr, U-Th-He, fission track and
#' detrital datasets
#'
#' \code{examples} a 14-item list containing:
#'
#' \code{UPb}: an object of class \code{UPb} containing a high
#' precision U-Pb dataset of Kamo et al. (1996) packaged with Ken
#' Ludwig's \code{Isoplot} program.
#'
#' \code{DZ}: an object of class \code{detrital} containing a detrital
#' zircon U-Pb dataset from Namibia (Vermeesch et al., 2015).
#'
#' \code{ArAr}: an object of class \code{ArAr} containing a
#' \eqn{^{40}}Ar/\eqn{^{39}}Ar spectrum of Skye basalt produced by Sarah
#' Sherlock (Open University).
#'
#' \code{UThHe}: an object of class \code{UThHe} containing a
#' U-Th-Sm-He dataset of Fish Lake apatite produced by Daniel Stockli
#' (UT Austin).
#'
#' \code{FT1}: an object of class \code{fissiontracks} containing a
#' synthetic external detector dataset.
#'
#' \code{FT2}: an object of class \code{fissiontracks} containing a
#' synthetic LA-ICP-MS-based fission track dataset using the zeta
#' calibration method.
#'
#' \code{FT3}: an object of class \code{fissiontracks} containing a
#' synthetic LA-ICP-MS-based fission track dataset using the absolute
#' dating approach.
#'
#' \code{ReOs}: an object of class \code{ReOs} containing a
#' \eqn{^{187}}Os/\eqn{^{187}}Re-dataset from Selby (2007).
#'
#' \code{SmNd}: an object of class \code{SmNd} containing a
#' \eqn{^{143}}Nd/\eqn{^{147}}Sm-dataset from Lugmair et al. (1975).
#'
#' \code{RbSr}: an object of class \code{RbSr} containing an
#' \eqn{^{87}}Rb/\eqn{^{86}}Sr-dataset from Compston et al. (1971).
#'
#' \code{Namib}: an object of class \code{detritals} containing a
#' detrital zircon U-Pb dataset of Vermeesch and Garzanti (2015)
#'
#' \code{average}: an object of class \code{other} containing the
#' \eqn{^{206}}Pb/\eqn{^{238}}U-ages and errors of dataset \code{UPb}.
#'
#' \code{KDE}: an object of class \code{'other'} containing the
#' \eqn{^{206}}Pb/\eqn{^{238}}U-ages (but not the errors) of dataset
#' \code{UPb}.
#'
#' \code{spectrum}: an object of class \code{'other'} containing the
#' \eqn{^{39}}Ar abundances, \eqn{^{40}}Ar/\eqn{^{39}}Ar-ages and errors of
#' dataset \code{ArAr}.
#'
#' \code{MountTom}: an object of class \code{'other'} containing a
#' dataset of dispersed zircon fission track ages from Brandon and
#' Vance (1992).
#' 
#' @name examples
#' @docType data
#' @examples
#' data(examples)
#' concordia(examples$UPb)
#'
#' agespectrum(examples$ArAr)
#'
#' isochron(examples$ReOs)
#'
#' radialplot(examples$FT1)
#'
#' helioplot(examples$UThHe)
#'
#' kde(examples$Namib)
#'
#' radialplot(examples$MountTom)
#'
#' agespectrum(examples$spectrum)
#'
#' weightedmean(examples$average)
#'
#' @references
#' Brandon, M.T. and Vance, J.A., 1992. Tectonic evolution of the
#' Cenozoic Olympic subduction complex, Washington State, as deduced
#' from fission track ages for detrial zircons. American Journal of
#' Science, 292, pp.565-565.
#'
#' Compston, W., Berry, H., Vernon, M.J., Chappell, B.W. and Kaye,
#' M.J., 1971. Rubidium-strontium chronology and chemistry of lunar
#' material from the Ocean of Storms. In Lunar and Planetary Science
#' Conference Proceedings (Vol. 2, p. 1471).
#'
#' Galbraith, R. F. and Green, P. F., 1990: Estimating the component
#' ages in a finite mixture, Nuclear Tracks and Radiation
#' Measurements, 17, 197-206.
#'
#' Kamo, S.L., Czamanske, G.K. and Krogh, T.E., 1996. A minimum U-Pb
#' age for Siberian flood-basalt volcanism. Geochimica et Cosmochimica
#' Acta, 60(18), 3505-3511.
#'
#' Ludwig, K. R., 2003. User's manual for Isoplot 3.00: a
#'     geochronological toolkit for Microsoft Excel. No. 4.
#'
#' Lugmair, G.W., Scheinin, N.B. and Marti, K., 1975. Sm-Nd age and
#' history of Apollo 17 basalt 75075-Evidence for early
#' differentiation of the lunar exterior. In Lunar and Planetary
#' Science Conference Proceedings (Vol. 6, pp. 1419-1429).
#'
#' Selby, D., 2007. Direct Rhenium-Osmium age of the
#' Oxfordian-Kimmeridgian boundary, Staffin bay, Isle of Skye, UK, and
#' the Late Jurassic time scale.  Norsk Geologisk Tidsskrift, 87(3),
#' p.291.
#'
#' Vermeesch, P. and Garzanti, E., 2015. Making geological sense of
#' 'Big Data' in sedimentary provenance analysis. Chemical Geology, 409,
#' pp.20-27.
#'
#' Vermeesch, P., 2008. Three new ways to calculate average (U-Th)/He ages.
#' Chemical Geology, 249(3),pp.339-347.
NULL
