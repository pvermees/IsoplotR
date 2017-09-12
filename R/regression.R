regression <- function(d,model=1,type='york',alpha=0.05){
    if (model==1 && identical(type,'york')){
        out <- york(d,alpha=alpha)
    } else if (model==1 && identical(type,'titterington')){
        out <- titterington(d,alpha=alpha)
    } else if (model==1 && identical(type,'ludwig')){
        out <- ludwig(d)
    } else if (model==2 && identical(type,'york')){
        out <- list()
        fit <- stats::lm(d[,'Y'] ~ d[,'X'])
        E <- stats::vcov(fit)
        tfact <- qt(1-alpha/2,fit$df.residual)
        out$a <- c(stats::coef(fit)[1],sqrt(E[1,1]),tfact*sqrt(E[1,1]))
        out$b <- c(stats::coef(fit)[2],sqrt(E[2,2]),tfact*sqrt(E[1,1]))
        names(out$a) <- c('a','s[a]','ci[a]')
        names(out$b) <- c('b','s[b]','ci[b]')
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
        out$err <- sqrt(diag(out$cov))
        tfact <- qt(1-alpha/2,out$df.residual)
        out$err <- rbind(tfact*out$err)
        colnames(out$err) <- parnames
        rownames(out$err) <- c('s','ci')
    } else {
        stop('invalid input to regression(d,model,type)')
    }
    out$model <- model
    out
}
