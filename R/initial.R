get.initial.ratio <- function(x,...){ UseMethod("get.initial.ratio",x) }
get.initial.ratio.default <- function(x,...){
    stop( "No default method available (yet)." )
}
get.initial.ratio.UPb <- function(x){
    lud <- ludwig(x)
    tt <- lud$par[1]
    if (x$format<4){
        rr <- age_to_terawasserburg_ratios(tt,d=x$d)$x
        tl <- lud$par['t[l]']
        y0 <- lud$par['76i']
        slope <- (rr['Pb207Pb206']-y0)/rr['U238Pb206']
        yp <- y0 + slope*get.U238Pb206.ratios(x)[,1]
        dy <- get.Pb207Pb206.ratios(x)[,1] - yp
        out <- y0 + dy
    } else {
        out <- matrix(0,length(x),2)
        # first determine 6/4-intercepts:
        y0 <- 1/lud$par['64i'] # 4/6-intercept
        x0 <- age_to_U238Pb206_ratio(tt,d=x$d)[1]   # 8/6-intercept
        X <- get.U238Pb206.ratios(x)[,1]            # 8/6 measurement
        yp <- y0 - X*y0/x0     # predicted 6/4-value
        Y <- X*get.Pb204U238.ratios(x)[,1]          # 4/6 measurement
        dy <- Y - yp           # 6/4 misfit
        out[,1] <- 1/(y0 + dy) # initial 6/4-ratio
        # then determine 7/4-intercepts:
        y0 <- 1/lud$par['74i'] # 4/7-intercept
        x0 <- 1/age_to_Pb207U235_ratio(tt,d=x$d)[1] # 5/7-intercept
        X <- 1/get.Pb207U235.ratios(x)[,1]          # 5/7 measurement
        yp <- y0 - X*y0/x0     # predicted 4/7-value
        U <- settings('iratio','U238U235')[1]
        Y <- X*U*get.Pb204U238.ratios(x)[,1]        # 4/7 measurement
        dy <- Y - yp           # 7/4 misfit
        out[,2] <- 1/(y0 + dy) # initial 7/4-ratio
    }
    out
}
