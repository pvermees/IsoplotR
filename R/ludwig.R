# Linear regression of X,Y,Z-variables with correlated errors, taking
# into account decay constant uncertainties
#
# Implements the maximum likelihood algorithm of Ludwig (1998)
#
# @param x a \eqn{3n}-element vector \eqn{[X Y Z]}, where \eqn{X},
#     \eqn{Y} and \eqn{Z} are three \eqn{n}-element vectors of
#     (isotopic ratio) values.
# @param covmat a \eqn{[3n x 3n]}-element covariance matrix of
#     \code{x}
# @references
# Ludwig, K.R., 1998. On the treatment of concordant uranium-lead
# ages. Geochimica et Cosmochimica Acta, 62(4), pp.665-676.
#
ludwig <- function(x,covmat){
    
}

data2ludwig <- function(x,...){ UseMethod("data2ludwig",x) }
data2ludwig.default <- function(x,...){ stop('default function undefined') }
data2ludwig.UPb <- function(x,...){
}
data2ludwig.ThU <- function(x,...){
    
}
