% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ludwig.R
\name{ludwig}
\alias{ludwig}
\title{Linear regression of U-Pb data with correlated errors, taking
into account decay constant uncertainties.}
\usage{
ludwig(
  x,
  model = 1,
  anchor = 0,
  exterr = FALSE,
  type = "joint",
  plot = FALSE,
  nsteps = NULL,
  ...
)

ludwig(
  x,
  model = 1,
  anchor = 0,
  exterr = FALSE,
  type = "joint",
  plot = FALSE,
  nsteps = NULL,
  ...
)
}
\arguments{
\item{x}{an object of class \code{UPb}}

\item{model}{one of three regression models:

\code{1}: fit a discordia_line through the data using the maximum
likelihood algorithm of Ludwig (1998), which assumes that the
scatter of the data is solely due to the analytical
uncertainties. In this case, \code{IsoplotR} will either calculate
an upper and lower intercept age (for Wetherill concordia), or a
lower intercept age and common
(\eqn{^{207}}Pb/\eqn{^{206}}Pb)\eqn{_\circ}-ratio intercept (for
Tera-Wasserburg). If the p-value for the chi-square test is less
than \code{alpha()}, then the analytical uncertainties are augmented
by a factor \eqn{\sqrt{MSWD}}.

\code{2}: fit a discordia_line ignoring the analytical uncertainties

\code{3}: fit a discordia_line using a modified maximum likelihood
algorithm that includes accounts for any overdispersion by adding a
geological (co)variance term.}

\item{anchor}{control parameters to fix the intercept age or common Pb
composition of the isochron fit. This can be a scalar or a vector.

If \code{anchor[1]=0}: do not anchor the isochron.

If \code{anchor[1]=1}: fix the common Pb composition at the values
stored in \code{settings('iratio',...)}.

If \code{anchor[1]=2}: force the isochron line to intersect the
concordia line at an age equal to \code{anchor[2]}.

If \code{anchor[1]=3}: anchor the isochron line to the
Stacey-Kramers mantle evolution model.}

\item{exterr}{propagate external sources of
uncertainty (i.e. decay constants)?}

\item{type}{only relevant if \code{x$format>3}. Can take on the following
values:

\code{'joint'} or \code{0}: 3-dimensional isochron regression.

\code{1}: 2-dimensional regression of
\eqn{{}^{204}}Pb/\eqn{{}^{206}}Pb vs.
\eqn{{}^{238}}U/\eqn{{}^{206}}Pb (for U-Pb formats 4, 5 and 6), or
of \eqn{{}^{208}}Pb/\eqn{{}^{206}}Pb vs.
\eqn{{}^{238}}U/\eqn{{}^{206}}Pb (for U-Pb formats 7 and 8).

\code{2}: 2-dimensional regression of
\eqn{{}^{204}}Pb/\eqn{{}^{207}}Pb vs.
\eqn{{}^{235}}U/\eqn{{}^{207}}Pb (for U-Pb formats 4, 5 and 6), or
of \eqn{{}^{208}}Pb/\eqn{{}^{207}}Pb vs.
\eqn{{}^{235}}U/\eqn{{}^{207}}Pb (for U-Pb formats 7 and 8).

\code{3}: 2-dimensional regression of
\eqn{{}^{206}}Pb/\eqn{{}^{208}}Pb vs.
\eqn{{}^{232}}Th/\eqn{{}^{208}}Pb (only for U-Pb formats 7 and 8).

\code{4}: 2-dimensional regression of
\eqn{{}^{207}}Pb/\eqn{{}^{208}}Pb vs.
\eqn{{}^{232}}Th/\eqn{{}^{208}}Pb (only for U-Pb formats 7 and 8).}

\item{plot}{logical. Only relevant for datasets with measured
disequilibrium. If \code{TRUE}, plots the posterior
distribution of the age and initial activity ratios.}

\item{nsteps}{resolution of the posterior distributions for
Bayesian credible intervals of disequilibrium-corrected U-Pb
isochrons}

\item{...}{optional arguments}
}
\value{
\describe{

\item{par}{a vector with the lower concordia intercept, the common
Pb ratios, the dispersion parameter (if \code{model=3}), and the
initial \eqn{{}^{234}}U/\eqn{{}^{238}}U and
\eqn{{}^{230}}Th/\eqn{{}^{238}}U activity ratio (in the presence of
initial disequilibrium).}

\item{cov}{the covariance matrix of \code{par}}

\item{df}{the degrees of freedom of the model fit (\eqn{n-2} if
\code{x$format<4} or \eqn{2n-3} if \code{x$format>3}, where \eqn{n}
is the number of aliquots).}

\item{mswd}{the mean square of weighted deviates (a.k.a. reduced
Chi-square statistic) for the fit.}

\item{p.value}{p-value of a Chi-square test for the linear fit}

}
}
\description{
Implements the maximum likelihood algorithm for Total-Pb/U isochron
regression of Ludwig (1998) and extends the underlying methodology
to accommodate U-Th-Pb data and initial U-series disequilibrium.
}
\details{
The 3-dimensional regression algorithm of Ludwig and Titterington
(1994) was modified by Ludwig (1998) to fit so-called `Total Pb-U
isochrons'. These are constrained to a radiogenic endmember
composition that falls on the \code{\link{concordia}} line. In its
most sophisticated form, this algorithm does not only allow for
correlated errors between variables, but also between
aliquots. \code{IsoplotR} currently uses this algorithm to
propagate decay constant uncertainties in the total Pb-U isochron
ages.
}
\examples{
f <- system.file("UPb4.csv",package="IsoplotR")
d <- read.data(f,method="U-Pb",format=4)
fit <- ludwig(d)

}
\references{
Ludwig, K.R., 1998. On the treatment of concordant uranium-lead
ages. Geochimica et Cosmochimica Acta, 62(4), pp.665-676.

Ludwig, K.R. and Titterington, D.M., 1994. Calculation of
\eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
Cosmochimica Acta, 58(22), pp.5031-5042.
}
\seealso{
\code{\link{concordia}}, \code{\link{titterington}},
    \code{\link{isochron}}
}
