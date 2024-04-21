get_ThPb_ratio <- function(tt,st,exterr=FALSE){
    getDPratio(tt=tt,st=st,nuclide='Th232',exterr=exterr)
}

get_ThPb_age <- function(Pb208Th232,sPb208Th232,exterr=FALSE){
    getPDage(DP=Pb208Th232,sDP=sPb208Th232,nuclide='Th232',exterr=exterr)
}

ThPb_age <- function(x,exterr=FALSE,i=NULL,i2i=TRUE,...){
    PD_age(x,nuclide='Th232',exterr=exterr,i=i,i2i=i2i,...)
}
