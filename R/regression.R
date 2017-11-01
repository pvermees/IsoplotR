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
        fit <- stats::lm(subset(d,select=c('Y','Z')) ~ d[,'X'])
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
    out <- list()
    if (identical(type,'york')){
        init <- york(d)
        out$w <- stats::optimize(york_disp_misfit,
                                 interval=c(0,init$mswd),
                                 d=d)$minimum
        dd <- augment_york_errors(d,out$w)
        out <- c(out,york(dd))
    } else if (identical(type,'titterington')){
        init <- titterington(d)
        out$w <- stats::optimize(titterington_disp_misfit,
                                 interval=c(0,init$mswd),
                                 d=d)$minimum
        dd <- augment_titterington_errors(d,out$w)
        out <- c(out,titterington(dd))
    } else {
        stop('invalid output type for model 3 regression')
    }
    out
}

york_disp_misfit <- function(w2,d){
    dd <- augment_york_errors(d,w2)
    abs(york(dd)$mswd - 1)
}
titterington_disp_misfit <- function(w2,d){
    dd <- augment_titterington_errors(d,w2)
    abs(titterington(dd)$mswd - 1)
}

augment_york_errors <- function(d,w){
    out <- d
    out[,'sY'] <- sqrt(d[,'sY']^2 + w^2)
    out[,'rXY'] <- d[,'rXY']*d[,'sY']/out[,'sY']
    out
}
augment_titterington_errors <- function(d,w){
    out <- d
    out[,'sZ'] <- sqrt(d[,'sZ']^2 + w^2)
    out[,'rXZ'] <- d[,'rXZ']*d[,'sZ']/out[,'sZ']
    out[,'rYZ'] <- d[,'rYZ']*d[,'sZ']/out[,'sZ']
    out
}
