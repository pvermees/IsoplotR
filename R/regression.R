regression <- function(d,model=1){
    if (model==1){
        out <- york(d)
    } else if (model==2){
        out <- list()
        fit <- lm(d[,'Y'] ~ d[,'X'])
        E <- vcov(fit)
        out$a <- c(coef(fit)[1],sqrt(E[1,1]))
        out$b <- c(coef(fit)[2],sqrt(E[2,2]))
        out$cov.ab <- E[1,2]
    } else {
        stop('model must be 1 or 2')
    }
    out$model <- model
    out
}
