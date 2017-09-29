regression <- function(d,model=1,type='york'){
    if (model==1) out <- model1regression(d,type=type)
    else if (model==2) out <- model2regression(d,type=type)
    else if (model==3) out <- model3regression(d,type=type)
    else stop('invalid regression model')
    out$model <- model
    out
}

model1regression <- function(d,type='york'){
    if (identical(type,'york')){
        out <- york(d)
    } else if (identical(type,'titterington')){
        out <- titterington(d)
    } else if (identical(type,'ludwig')){
        out <- ludwig(d)
    } else {
        stop('invalid output type for model 1 regression')
    }
    out
}

model2regression <- function(d,type='york'){
    if (identical(type,'york')){
        out <- list()
        fit <- stats::lm(d[,'Y'] ~ d[,'X'])
        E <- stats::vcov(fit)
        out$df <- fit$df.residual
        out$a <- c(stats::coef(fit)[1],sqrt(E[1,1]))
        out$b <- c(stats::coef(fit)[2],sqrt(E[2,2]))
        names(out$a) <- c('a','s[a]')
        names(out$b) <- c('b','s[b]')
        out$cov.ab <- E[1,2]        
    } else if (identical(type,'titterington')){
        out <- list()
        fit <- stats::lm(d[,c('Y','Z')] ~ d[,'X'])
        out$df <- fit$df.residual
        out$par <- c(stats::coef(fit))
        out$cov <- stats::vcov(fit)
        parnames <- c('a','b','A','B')
        names(out$par) <- parnames
        colnames(out$cov) <- parnames
        rownames(out$cov) <- parnames
    } else {
        stop('invalid output type for model 2 regression')
    }
    out
}

model3regression <- function(d,type='york'){
    if (identical(type,'york')){
        init <- york(d)$mswd
        w <- optimize(york_disp_misfit,interval=c(0,init),d=d)$minimum
        dd <- augment_errors(d,w)
        out <- york(dd)
    } else if (identical(type,'titterington')){
        out <- titterington(d)
    } else if (identical(type,'ludwig')){
        out <- ludwig(d)
    } else {
        stop('invalid output type for model 3 regression')
    }
    out
}

york_disp_misfit <- function(w2,d){
    dd <- augment_errors(d,w2)
    abs(york(dd)$mswd - 1)
}

augment_errors <- function(d,w){
    out <- d
    out[,'sY'] <- sqrt(d[,'sY']^2 + w^2)
    out[,'rXY'] <- d[,'rXY']*d[,'sY']/out[,'sY']
    out
}

foo <- function(d,interval){
    w2 <- seq(from=interval[1],to=interval[2],length.out=100)
    y <- rep(0,100)
    for (i in 1:length(w2)){
        y[i] <- york_disp_misfit(w2[i],d)
    }
    plot(w2,y,type='l')
}
