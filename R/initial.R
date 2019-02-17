get.initial.ratio <- function(x,...){ UseMethod("get.initial.ratio",x) }
get.initial.ratio.default <- function(x,...){
    stop( "No default method available (yet)." )
}
get.initial.ratio.UPb <- function(x){
    lud <- ludwig(x)
    tt <- lud$par[1]
    if (x$format<4){
        rr <- age_to_terawasserburg_ratios(tt,d=x$d)$x
        slope <- (rr['Pb207Pb206']-lud$par['76i'])/rr['U238Pb206']
        out <- get.Pb207Pb206.ratios(x)[,1] -
               slope*get.U238Pb206.ratios(x)[,1]
    } else {
        out <- matrix(0,length(x),2)
        # first determine 6/4-intercepts:
        y0 <- 1/lud$par['64i']     # 4/6-intercept
        x0 <- age_to_U238Pb206_ratio(tt,d=x$d)[1]   # 8/6-intercept
        X <- get.U238Pb206.ratios(x)[,1]            # 8/6 measurement
        Y <- X*get.Pb204U238.ratios(x)[,1]          # 4/6 measurement
        out[,1] <- 1/(Y + X*y0/x0) # initial 6/4-ratio
        # then determine 7/4-intercepts:
        y0 <- 1/lud$par['74i'] # 4/7-intercept
        x0 <- 1/age_to_Pb207U235_ratio(tt,d=x$d)[1] # 5/7-intercept
        X <- 1/get.Pb207U235.ratios(x)[,1]          # 5/7 measurement
        U <- settings('iratio','U238U235')[1]
        Y <- X*U*get.Pb204U238.ratios(x)[,1]        # 4/7 measurement
        out[,2] <- 1/(Y + X*y0/x0) # initial 7/4-ratio
    }
    out
}
get.initial.ratio.ArAr <- function(x){
    get.initial.ratio_helper(x)
}
get.initial.ratio.KCa <- function(x){
    get.initial.ratio.PD(x)
}
get.initial.ratio.PD <- function(x){
    out <- get.initial.ratio_helper(x)
    list(y0=out[,1],sy0=out[,2])
}
get.initial.ratio_helper <- function(x){
    y <- data2york(x,inverse=TRUE)
    fit <- regression(y,model=1)
    yi <- 1/(y[,'Y'] - fit$b[1]*y[,'X'])
    syi <- fit$b[2]*y[,'X']*yi^2
    cbind(yi,syi)
}
