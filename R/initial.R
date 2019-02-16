get.initial.ratio <- function(x,...){ UseMethod("get.initial.ratio",x) }
get.initial.ratio.default <- function(x,...){
    stop( "No default method available (yet)." )
}
get.initial.ratio.UPb <- function(x){
    lud <- ludwig(x)
    if (x$format<4){
        tl <- lud$par['t[l]']
        y0 <- lud$par['76i']
        rr <- age_to_terawasserburg_ratios(tl,d=x$d)$x
        slope <- (rr['Pb207Pb206']-y0)/rr['U238Pb206']
        yp <- y0 + slope*get.U238Pb206.ratios(x)[,1]
        dy <- get.Pb207Pb206.ratios(x)[,1] - yp
        out <- y0 + dy
    } else {
        out <- lud$par[c('64i','74i')]
    }
    out
}
