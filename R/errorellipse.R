#' Get coordinates of error ellipse for plotting
#'
#' Construct an error ellipse age a given confidence level from its
#' centre and covariance matrix
#'
#' @param x x-coordinate (scalar) for the centre of the ellipse
#' @param y y-coordinate (scalar) for the centre of the ellipse
#' @param covmat covariance matrix of the x-y coordinates
#' @param alpha the probability cutoff for the error ellipses
#' @return a [50x2] matrix of plot coordinates
#' @examples
#' x = 99; y = 101;
#' covmat <- matrix(c(1,0.9,0.9,1),nrow=2)
#' ell <- get.ellipse(x,y,covmat)
#' plot(c(90,110),c(90,110),type='l')
#' polygon(ell,col=rgb(0,1,0,0.5))
#' points(x,y,pch=21,bg='black')
#' @export
get.ellipse <- function(x,y,covmat,alpha=0.05){
    nn <- 50
    cutoff <- stats::qchisq(1-alpha,2)
    e <- eigen(covmat)
    a <- sqrt(cutoff*e$values[1]) # major axis
    b <- sqrt(cutoff*e$values[2]) # minor axis
    v <- e$vectors[,1] # largest eigenvector
    alpha <- atan(v[2]/v[1]) # rotation angle of the ellipse
    theta <- seq(0, 2 * pi, length=nn)
    out <- matrix(0,nrow=nn,ncol=2)
    out[,1] <- x + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
    out[,2] <- y + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
    colnames(out) <- c('x','y')
    out
}
