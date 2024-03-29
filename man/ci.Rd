% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci.R
\name{ci}
\alias{ci}
\title{Confidence intervals}
\usage{
ci(x = 0, sx, oerr = 3, df = NULL, absolute = FALSE)
}
\arguments{
\item{x}{scalar estimate}

\item{sx}{scalar or vector with the standard error of x without and
(optionally) with \eqn{\sqrt{MSWD}} overdispersion multiplier.}

\item{oerr}{indicates if the confidence intervals should be
    reported as:

\code{1}: 1\eqn{\sigma} absolute uncertainties.

\code{2}: 2\eqn{\sigma} absolute uncertainties.

\code{3}: absolute (1-\eqn{\alpha})\% confidence intervals, where
\eqn{\alpha} equales the value that is stored in
\code{settings('alpha')}.

\code{4}: 1\eqn{\sigma} relative uncertainties (\eqn{\%}).

\code{5}: 2\eqn{\sigma} relative uncertainties (\eqn{\%}).

\code{6}: relative (1-\eqn{\alpha})\% confidence intervals, where
\eqn{\alpha} equales the value that is stored in
\code{settings('alpha')}.}

\item{df}{(optional) number of degrees of freedom. Only used if
\code{sx} is a vector.}

\item{absolute}{logical. Returns absolute uncertainties even if
\code{oerr} is greater than 3. Used for some internal
\code{IsoplotR} functions.}
}
\value{
A scalar or vector of the same size as \code{sx}.
}
\description{
Given a parameter estimate and its standard error,
    calculate the corresponding 1-sigma, 2-sigma or
    \eqn{100(1-\alpha)\%} confidence interval, in absolute or
    relative units.
}
\details{
Several of \code{IsoplotR}'s plotting functions (including
    \code{\link{isochron}}, \code{\link{weightedmean}},
    \code{\link{concordia}}, \code{\link{radialplot}} and
    \code{\link{helioplot}}) return lists of parameters and their
    standard errors. For `model-1' fits, if the data pass a
    Chi-square test of homogeneity, then just one estimate for the
    standard error is reported.  This estimate can be converted to
    a confidence interval by multiplication with the appropriate
    quantile of a Normal distribution. Datasets that fail the
    Chi-square test are said to be `overdispersed' with respect to
    the analytical uncertainties. The simplest way (`model-1') to
    deal with overdispersion is to inflate the standard error with
    a \eqn{\sqrt{MSWD}} premultiplier. In this case,
    \code{IsoplotR} returns two estimates of the standard error.
    To convert the second estimate to a confidence interval
    requires multiplication with the desired quantile of a
    t-distribution with the appropriate degrees of freedom.
}
\examples{
attach(examples)
fit <- isochron(PbPb,plot=FALSE,exterr=FALSE)
err <- ci(x=fit$age[1],sx=fit$age[-1],oerr=5,df=fit$df)
message('age=',signif(fit$age[1],4),'Ma, ',
        '2se=',signif(err[1],2),'\%, ',
        '2se(with dispersion)=',signif(err[2],2),'\%')
}
