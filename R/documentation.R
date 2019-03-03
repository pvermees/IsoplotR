#' Example datasets for testing \code{IsoplotR}
#'
#' U-Pb, Pb-Pb, Ar-Ar, K-Ca, Re-Os, Sm-Nd, Rb-Sr, Lu-Hf, U-Th-He, Th-U,
#' fission track and detrital datasets
#'
#' \code{examples} an 18-item list containing:
#'
#' \code{UPb}: an object of class \code{UPb} containing a high
#' precision U-Pb dataset of Kamo et al. (1996) packaged with Ken
#' Ludwig (2003)'s \code{Isoplot} program.
#'
#' \code{PbPb}: an object of class \code{PbPb} containing a Pb-Pb
#' dataset from Connelly et al. (2017).
#'
#' \code{DZ}: an object of class \code{detrital} containing a detrital
#' zircon U-Pb dataset from Namibia (Vermeesch et al., 2015).
#'
#' \code{ArAr}: an object of class \code{ArAr} containing a
#' \eqn{^{40}}Ar/\eqn{^{39}}Ar spectrum of Skye basalt produced by Sarah
#' Sherlock (Open University).
#'
#' \code{KCa}: an object of class \code{KCa} containing a
#' \eqn{^{40}}K/\eqn{^{40}}Ca dataset for sample 140025 grain h spot 5
#' of Harrison et al. (2010).
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
#' \code{LuHf}: an object of class \code{LuHf} containing an
#' \eqn{^{176}}Lu/\eqn{^{177}}Hf-dataset from Barfod et al. (2002).
#'
#' \code{ThU}: an object of class \code{ThU} containing a synthetic
#' `Osmond-type' dataset from Titterington and Ludwig (1994).
#'
#' \code{LudwigMean}: an object of class \code{other} containing a
#' collection of \eqn{^{206}}Pb/\eqn{^{238}}U-ages and errors of the
#' example dataset by Ludwig (2003).
#'
#' \code{LudwigKDE}: an object of class \code{'other'} containing the
#' \eqn{^{206}}Pb/\eqn{^{238}}U-ages (but not the errors) of the
#' example dataset by Ludwig (2003).
#'
#' \code{LudwigSpectrum}: an object of class \code{'other'} containing
#' the \eqn{^{39}}Ar abundances, \eqn{^{40}}Ar/\eqn{^{39}}Ar-ages and
#' errors of the example dataset by Ludwig (2003).
#'
#' \code{LudwigMixture}: an object of class \code{'other'} containing
#' a dataset of dispersed zircon fission track ages of the example
#' dataset by Ludwig (2003).
#'
#' @name examples
#' @docType data
#' @examples
#' data(examples)
#'
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
#' evolution(examples$ThU)
#'
#' kde(examples$DZ)
#'
#' radialplot(examples$LudwigMixture)
#'
#' agespectrum(examples$LudwigSpectrum)
#'
#' weightedmean(examples$LudwigMean)
#'
#' @references
#' Barfod, G.H., Albarede, F., Knoll, A.H., Xiao, S., Telouk, P.,
#' Frei, R. and Baker, J., 2002. New Lu-Hf and Pb-Pb age constraints on
#' the earliest animal fossils. Earth and Planetary Science Letters, 201(1),
#' pp.203-212.
#'
#' Compston, W., Berry, H., Vernon, M.J., Chappell, B.W. and Kaye,
#' M.J., 1971. Rubidium-strontium chronology and chemistry of lunar
#' material from the Ocean of Storms. In Lunar and Planetary Science
#' Conference Proceedings (Vol. 2, p. 1471).
#'
#' Connelly, J.N., Bollard, J. and Bizzarro, M., 2017. Pb-Pb
#' chronometry and the early Solar System. Geochimica et Cosmochimica
#' Acta, 201, pp.345-363.
#'
#' Galbraith, R. F. and Green, P. F., 1990: Estimating the component
#' ages in a finite mixture, Nuclear Tracks and Radiation
#' Measurements, 17, 197-206.
#'
#' Harrison, T.M., Heizler, M.T., McKeegan, K.D. and Schmitt, A.K.,
#' 2010. In situ \eqn{^{40}}K-\eqn{^{40}}Ca `double-plus' SIMS dating
#' resolves Klokken feldspar \eqn{^{40}}K-\eqn{^{40}}Ar paradox. Earth
#' and Planetary Science Letters, 299(3-4), pp.426-433.
#'
#' Kamo, S.L., Czamanske, G.K. and Krogh, T.E., 1996. A minimum U-Pb
#' age for Siberian flood-basalt volcanism. Geochimica et Cosmochimica
#' Acta, 60(18), 3505-3511.
#'
#' Ludwig, K. R., and D. M. Titterington., 1994.
#' "Calculation of \eqn{^{230}}Th/U isochrons, ages, and errors."
#' Geochimica et Cosmochimica Acta 58.22, 5031-5042.
#'
#' Ludwig, K. R., 2003. User's manual for Isoplot 3.00: a
#' geochronological toolkit for Microsoft Excel. No. 4.
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
#' `Big Data' in sedimentary provenance analysis. Chemical Geology, 409,
#' pp.20-27.
#'
#' Vermeesch, P., 2008. Three new ways to calculate average (U-Th)/He ages.
#' Chemical Geology, 249(3),pp.339-347.
NULL
