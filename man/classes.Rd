% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{classes}
\alias{classes}
\alias{UPb}
\alias{PbPb}
\alias{ThPb}
\alias{ThU}
\alias{ArAr}
\alias{KCa}
\alias{PD}
\alias{RbSr}
\alias{SmNd}
\alias{LuHf}
\alias{ReOs}
\alias{UThHe}
\alias{fissiontracks}
\alias{detritals}
\alias{is.UPb}
\alias{is.PbPb}
\alias{is.ThPb}
\alias{is.ThU}
\alias{is.ArAr}
\alias{is.KCa}
\alias{is.PD}
\alias{is.RbSr}
\alias{is.SmNd}
\alias{is.LuHf}
\alias{is.ReOs}
\alias{is.UThHe}
\alias{is.fissiontracks}
\alias{is.detritals}
\alias{is.diseq}
\title{Geochronological data classes}
\usage{
is.UPb(x)

is.PbPb(x)

is.ThPb(x)

is.ArAr(x)

is.KCa(x)

is.PD(x)

is.RbSr(x)

is.SmNd(x)

is.LuHf(x)

is.ReOs(x)

is.ThU(x)

is.UThHe(x)

is.fissiontracks(x)

is.detritals(x)

is.diseq(x)
}
\arguments{
\item{x}{a data object returned by \code{\link{read.data}} or
\code{\link{diseq}}.}
}
\value{
logical
}
\description{
S3 classes to store geochronological data generated by
    \code{\link{read.data}} or \code{\link{diseq}}.
}
\details{
\code{IsoplotR} uses the following S3 classes to store
    geochronological data: \code{UPb}, \code{PbPb}, \code{ThPb},
    \code{KCa}, \code{UThHe}, \code{fissiontracks},
    \code{detritals} and \code{PD}, where the latter is the parent
    class of the simple parent-daughter chronometers \code{RbSr},
    \code{SmNd}, \code{LuHf} and \code{ReOs}. All these classes
    have overloaded versions of the generic \code{length()}
    function and \code{`[`} subsetting method.

\itemize{

\item\code{UPb}: a list containing:

\describe{

\item{\code{x}}{a matrix containing the isotopic measurements}

\item{\code{format}}{a number between 1 and 8}

\item{\code{d}}{an object of class \code{diseq}, i.e. the output of
the \code{\link{diseq}} function}

}

\item\code{ArAr}: a list containing

\describe{

\item{\code{x}}{a matrix containing the isotopic measurements}

\item{\code{J}}{a two-element vector with the J-factor and its
uncertainty}

\item{\code{format}}{a number between 1 and 3}

}

\item\code{ThU}: a list containing

\describe{

\item{\code{x}}{a matrix containing the isotopic measurements}

\item{\code{format}}{a number between 1 and 4}

\item{\code{Th02}}{a two element vector with the assumed initial
    \eqn{^{230}}Th/\eqn{^{232}}Th-ratio of Th-bearing detritus (for
    formats 1 and 2) or rock (for formats 3 and 4) and its standard
    error}

\item{\code{Th02U48}}{9-element vector with the measured
    composition of Th-bearing detritus}

}

\item\code{PbPb}, \code{ThPb}, \code{KCa}, \code{PD}, \code{RbSr},
\code{SmNd}, \code{LuHf}, or \code{ReOs}: a list containing

\describe{

\item{\code{x}}{a matrix containing the isotopic measurements}

\item{\code{format}}{a number between 1 and 3}

}

\item\code{UThHe}: a matrix of He, U, Th (and Sm) measurements

\item\code{fissiontracks}: a list containing

\describe{

\item{\code{format}}{a number between 1 and 3}

\item{\code{x}}{a matrix of spontaneous and induced fission track
counts (only included if \code{format=1})}

\item{\code{rhoD}}{the track density of the dosimeter glass,
extracted from the input data (only included if \code{format=1})}

\item{\code{zeta}}{the zeta calibration constant extracted from the
input data (only included if \code{format<3})}

\item{\code{Ns}}{a list containing the spontaneous fission track
counts (only included if \code{format>1})}

\item{\code{U}}{a list of lists containing the U-concentration or
U/Ca-ratio measurements for each of the analysed grains (only
included if \code{format>1})}

\item{\code{sU}}{a list of lists containing the standard errors of
the U-concentration or U/Ca-ratio measurements for each of the
analysed grains (only include if \code{format>1})}

\item{\code{spotSize}}{the laser ablation spot size (only included
if \code{format>1})}

}

\item\code{detritals}: a list of named vectors, one for each
detrital sample.

\item\code{diseq}: is a class that contains the output of the
\code{\link{diseq}} function, which stores initial disequilibrium
data for U--Pb geochronology.

}
}
\examples{
attach(examples)
ns <- length(UPb)
concordia(UPb[-ns,])
if (is.PD(RbSr)) print('RbSr has class PD')

}
\seealso{
read.data diseq
}
