get.ThPb.ratio <- function(tt,st,exterr=TRUE){
    get.PD.ratio(tt=tt,st=st,nuclide='Th232',exterr=exterr)
}

get.ThPb.age <- function(Th232Pb208,sTh232Pb208,exterr=TRUE){
    get.PD.age(DP=Th232Pb208,sDP=sTh232Pb208,nuclide='Th232',exterr=exterr)
}

ThPb.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,i2i=TRUE,...){
    PD.age(x,nuclide='Th232',exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
}
