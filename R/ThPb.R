get.ThPb.ratio <- function(tt,st,exterr=FALSE){
    get.PD.ratio(tt=tt,st=st,nuclide='Th232',exterr=exterr)
}

get.ThPb.age <- function(Pb208Th232,sPb208Th232,exterr=FALSE){
    get.PD.age(DP=Pb208Th232,sDP=sPb208Th232,nuclide='Th232',exterr=exterr)
}

ThPb.age <- function(x,exterr=FALSE,i=NULL,i2i=TRUE,...){
    PD.age(x,nuclide='Th232',exterr=exterr,i=i,i2i=i2i,...)
}
