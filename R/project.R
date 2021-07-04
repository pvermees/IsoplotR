project <- function(x,...){ UseMethod("project",x) }
project.default <- function(x,...){stop('No default project method implemented.')}
# type = 1: nominal
# type = 2: isochron
# type = 3: Stacey-Kramers
project.UPb <- function(x,type=1,inverse=TRUE,projerr=FALSE,tt=0,...)
    for (i in 1:length(x)){
        if (x.format%in%c(7,8)){
            X <- get.UPb.isochron.ratios.208(x,i=i,tt=tt)
        }
    }
}
