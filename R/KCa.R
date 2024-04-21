get_KCa_ratio <- function(tt,st,exterr=FALSE,bratio=0.895){
    getDPratio(tt,st,'K40',exterr,bratio=bratio)
}

get_KCa_age <- function(K40Ca40,sK40Ca40,exterr=FALSE,bratio=0.895){
    getPDage(K40Ca40,sK40Ca40,'K40',exterr=exterr,bratio=bratio)
}

KCa_age <- function(x,exterr=FALSE,i=NULL,i2i=TRUE,
                    omit4c=NULL,projerr=FALSE,bratio=0.895){
    PD_age(x,'K40',exterr=exterr,i=i,i2i=i2i,
           bratio=bratio,omit4c=omit4c,projerr=projerr)
}
