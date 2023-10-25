anchor2anchor <- function(x,anchor=0,inverse=TRUE){
    if (anchor[1]==1 && length(anchor)>1){
        out <- common2anchor(x,inverse=inverse)
    } else if (anchor[1]==2 && length(anchor)>1){
        out <- age2anchor(x,anchor=anchor,inverse=inverse)
    } else {
        out <- anchor
    }
    out
}

common2anchor <- function(x,...){ UseMethod("common2anchor",x) }
common2anchor.default <- function(x,...){ stop("Not implemented") }
common2anchor.PbPb <- function(x,inverse,...){
    Pb74 <- settings('iratio','Pb207Pb204')[1]
    if (inverse) out <- c(2,Pb74)
    else out <- c(1,Pb74)
    out
}

age2anchor <- function(x,...){ UseMethod("age2anchor",x) }
age2anchor.default <- function(x,anchor,inverse,...){ anchor }
age2anchor.PbPb <- function(x,anchor,inverse,...){
    Pb76 <- age_to_Pb207Pb206_ratio(anchor[2])
    if (inverse) out <- c(1,Pb76)
    else out <- c(2,Pb76)
    out
}
