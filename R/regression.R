regression <- function(d,model=1,type='york'){
    if (model==1 && identical(type,'york')){
        out <- york(d)
    } else if (model==1 && identical(type,'titterington')){
        out <- titterington(d)
    } else if (model==1 && identical(type,'ludwig')){
        out <- ludwig(d)
    } else if (model==2 && identical(type,'york')){
        out <- list()
        fit <- stats::lm(d[,'Y'] ~ d[,'X'])
        E <- stats::vcov(fit)
        out$a <- c(stats::coef(fit)[1],sqrt(E[1,1]))
        out$b <- c(stats::coef(fit)[2],sqrt(E[2,2]))
        out$cov.ab <- E[1,2]        
    } else if (model==2 && identical(type,'titterington')){
        out <- list()
        fit <- stats::lm(d[,c('Y','Z')] ~ d[,'X'])
        out$par <- c(stats::coef(fit))
        out$cov <- stats::vcov(fit)
        parnames <- c('a','b','A','B')
        names(out$par) <- parnames
        colnames(out$cov) <- parnames
        rownames(out$cov) <- parnames
    } else {
        stop('invalid input to regression(d,model,type)')
    }
    out$model <- model
    out
}
