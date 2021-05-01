#' @export
is.UPb <- function(x) inherits(x,"UPb")
#' @export
is.PbPb <- function(x) inherits(x,"PbPb")
#' @export
is.ThPb <- function(x) inherits(x,"ThPb")
#' @export
is.ArAr <- function(x) inherits(x,"ArAr")
#' @export
is.KCa <- function(x) inherits(x,"KCa")
#' @export
is.PD <- function(x) inherits(x,"PD")
#' @export
is.RbSr <- function(x) inherits(x,"RbSr")
#' @export
is.SmNd <- function(x) inherits(x,"SmNd")
#' @export
is.LuHf <- function(x) inherits(x,"LuHf")
#' @export
is.ReOs <- function(x) inherits(x,"ReOs")
#' @export
is.ThU <- function(x) inherits(x,"ThU")
#' @export
is.UThHe <- function(x) inherits(x,"UThHe")
#' @export
is.fissiontracks <- function(x) inherits(x,"fissiontracks")
#' @export
is.detritals <- function(x) inherits(x,"detritals")
#' @export
is.diseq <- function(x) inherits(x,"diseq")

#' @export
length.UPb  <- function(x){ nrow(x$x) }
#' @export
length.PbPb <- function(x){ nrow(x$x) }
#' @export
length.ArAr <- function(x){ nrow(x$x) }
#' @export
length.PD <- function(x){ nrow(x$x) }
#' @export
length.ThPb <- function(x){ nrow(x$x) }
#' @export
length.KCa <- function(x){ nrow(x$x) }
#' @export
length.ThU <- function(x){ nrow(x$x) }
#' @export
length.UThHe <- function(x){ nrow(x) }
#' @export
length.KDE <- function(x){ length(x$ages) }
#' @export
length.KDEs <- function(x){ length(x$kdes) }
#' @export
length.fissiontracks <- function(x){
    if (x$format==1) return(nrow(x$x))
    else return(length(x$Ns))
}

#' @export
`[.UPb` <- function(x,...){
    out <- x
    out$x <- x$x[...]
    if ('x.raw' %in% names(x)){
        out$x.raw <- x$x.raw[...]
    }
    out    
}
#' @export
`[.ArAr` <- function(x,...){ `[helper`(x,...) }
#' @export
`[.PbPb` <- function(x,...){ `[helper`(x,...) }
#' @export
`[.ThPb` <- function(x,...){ `[helper`(x,...) }
#' @export
`[.KCa` <- function(x,...){ `[helper`(x,...) }
#' @export
`[.PD` <- function(x,...){ `[helper`(x,...) }
#' @export
`[.ThU` <- function(x,...){ `[helper`(x,...) }
#' @export
`[.fissiontracks` <- function(x,...){
    if (x$format==1){
        out <- `[helper`(x,...)
    } else {
        out <- x
        out$Ns <- x$Ns[...]
        out$A <- x$A[...]
        out$U <- x$U[...]
        out$sU <- x$sU[...]
    }
    out
}
#' @export
`[.diseq` <- function(x,i){
    out <- x
    for (ratio in c('U48','ThU','RaU','PaU')){
        j <- min(length(out[[ratio]]$x),i)
        out[[ratio]]$x <- out[[ratio]]$x[j]
    }
    out
}
`[helper` <- function(x,...){
    out <- x
    out$x <- x$x[...]
    out
}

#' @export
subset.UPb  <- function(x,...){
    out <- x
    out$x <- subset(x$x,...)
    if ('x.raw' %in% names(x)){
        out$x.raw <- subset(x$x.raw,...)
    }
    out
}
#' @export
subset.PbPb <- function(x,...){ subset_helper(x,...) }
#' @export
subset.ArAr <- function(x,...){ subset_helper(x,...) }
#' @export
subset.ThPb <- function(x,...){ subset_helper(x,...) }
#' @export
subset.KCa <- function(x,...){ subset_helper(x,...) }
#' @export
subset.PD <- function(x,...){ subset_helper(x,...) }
#' @export
subset.ThU <- function(x,...){ subset_helper(x,...) }
#' @export
subset.UThHe <- function(x,...){
    out <- subset.matrix(x,...)
    class(out) <- class(x)
    out
}
#' @export
subset.fissiontracks <- function(x,...){
    if (x$format==1){
        out <- subset_helper(x,...)
    } else {
        out <- x
        out$Ns <- subset(x$Ns,...)
        out$A <- subset(x$A,...)
        out$U <- subset(x$U,...)
        out$sU <- subset(x$sU,...)
    }
    out
}
subset_helper <- function(x,...){
    out <- x
    out$x <- subset.matrix(x$x,...)
    out
}
