regression <- function(xyz,model=1,type='york',omit=NULL,wtype='a'){
    xyz2calc <- clear(xyz,omit)
    if (model==1){
        out <- model1regression(xyz2calc,type=type)
    } else if (model==2){
        out <- model2regression(xyz2calc,type=type)
    } else if (model==3){
        out <- model3regression(xyz2calc,type=type,wtype=wtype)
        out$wtype <- wtype
    } else {
        stop('invalid regression model')
    }
    out$xyz <- xyz
    out$model <- model
    out$n <- nrow(xyz2calc)
    out$omit <- omit
    out
}

model1regression <- function(xyz,type='york'){
    if (identical(type,'york')){
        out <- york(xyz)
    } else if (identical(type,'titterington')){
        out <- titterington(xyz)
    } else {
        stop('invalid output type for model 1 regression')
    }
    out
}

model2regression <- function(xyz,type='york'){
    if (identical(type,'york')){
        out <- tls(xyz[,c('X','Y')])
        out$df <- nrow(xyz)-2
        out$a <- c(out$par['a'],'s[a]'=unname(sqrt(out$cov['a','a'])))
        out$b <- c(out$par['b'],'s[b]'=unname(sqrt(out$cov['b','b'])))
        out$cov.ab <- out$cov['a','b']
    } else if (identical(type,'titterington')){
        out <- tls(xyz[,c('X','Y','Z')])
        out$df <- 2*nrow(xyz)-4
    } else {
        stop('invalid output type for model 2 regression')
    }
    out
}

# fixes signs and uses logs for numerical stability:
model3regression <- function(xyz,type='york',model=3,wtype='a'){
    pilot <- model1regression(xyz,type=type)
    fact <- max(1,sqrt(pilot$mswd))
    if (identical(type,'york')){
        err <- ifelse(wtype%in%c('intercept',0,'a'),
                      pilot$a['s[a]'],pilot$b['s[b]'])
        lw <- log(fact*err)
        init <- c(pilot$a['a'],pilot$b['b'],lw=unname(lw))
        lower <- init - c(20*fact*pilot$a['s[a]'],20*fact*pilot$b['s[b]'],20)
        upper <- init + c(20*fact*pilot$a['s[a]'],20*fact*pilot$b['s[b]'],2)
        ps <- getparscale(fn=LL.york,args=list(ablw=init,XY=xyz,wtype=wtype))
        out <- contingencyfit(init,LL.york,lower=lower,upper=upper,
                              XY=xyz,wtype=wtype,control=list(parscale=ps))
        out$cov <- inverthess(out$hessian)
        out$a <- c('a'=unname(out$par['a']),'s[a]'=unname(sqrt(out$cov['a','a'])))
        out$b <- c('b'=unname(out$par['b']),'s[b]'=unname(sqrt(out$cov['b','b'])))
        out$cov.ab <- out$cov['a','b']
    } else if (identical(type,'titterington')){
        spar <- fact*sqrt(diag(pilot$cov))
        if (wtype%in%c('intercept',0,'a')) lw <- log(spar['a'])
        else if (wtype%in%c(1,'b')) lw <- log(spar['b'])
        else if (wtype%in%c(2,'A')) lw <- log(spar['A'])
        else if (wtype%in%c(3,'B')) lw <- log(spar['B'])
        else stop('illegal wtype')
        init <- c(pilot$par,'lw'=unname(lw))
        lower <- c(pilot$par-10*spar,'lw'=init['lw']-10)
        upper <- c(pilot$par+10*spar,'lw'=init['lw']+2)
        out <- contingencyfit(init,LL.titterington,lower=lower,
                              upper=upper,XYZ=xyz,wtype=wtype)
        x <- get.titterington.xyz(XYZ=xyz,a=out$par['a'],b=out$par['b'],
                                  A=out$par['A'],B=out$par['B'],
                                  w=exp(out$par['lw']),wtype=wtype)[,'x']
        H <- stats::optimHess(par=c(out$par,x),fn=LL.titterington.abABlwx,
                              XYZ=xyz,wtype=wtype)
        if (invertible(H)) out$cov <- inverthess(H)[1:5,1:5]
        else if (invertible(out$hessian)) out$cov <- inverthess(out$hessian)
        else out$cov <- inverthess(H)[1:5,1:5]
    } else {
        stop('invalid output type for model 3 regression')
    }
    disp <- exp(out$par['lw'])
    sdisp <- disp*sqrt(out$cov['lw','lw'])
    out$disp <- c('w'=unname(disp),'s[w]'=unname(sdisp))
    out
}

LL.york <- function(ablw,XY,wtype='a'){
    x <- get.york.xy(XY=XY,a=ablw['a'],b=ablw['b'],
                     w=exp(ablw['lw']),wtype=wtype)[,'x']
    LL.york.ablwx(c(ablw,x),XY,wtype=wtype)
}
LL.york.ablwx <- function(ablwx,XY,wtype='a'){
    ns <- nrow(XY)
    a <- ablwx[1]
    b <- ablwx[2]
    w <- exp(ablwx[3])
    x <- ablwx[4:(ns+3)]
    DE <- matrix(0,nrow=ns,ncol=5)
    colnames(DE) <- c('X-x','Y-y','vX','vY','sXY')
    DE[,'X-x'] <- XY[,'X']-x
    DE[,'Y-y'] <- XY[,'Y']-a-b*x
    DE[,'vX'] <- XY[,'sX']^2
    if (wtype%in%c('intercept',0,'a')) DE[,'vY'] <- XY[,'sY']^2 + w^2
    else if (wtype%in%c('slope',1,'b')) DE[,'vY'] <- XY[,'sY']^2 + (w*x)^2
    else DE[,'vY'] <- XY[,'sY']^2
    DE[,'sXY'] <- XY[,'rXY']*XY[,'sX']*XY[,'sY']
    detE <- DE[,'vX']*DE[,'vY'] - DE[,'sXY']^2
    O <- invertcovmat(vx=DE[,'vX'],vy=DE[,'vY'],sxy=DE[,'sXY'])
    maha <- (O[,'xx']*DE[,'X-x'] + O[,'xy']*DE[,'Y-y'])*DE[,'X-x'] +
        (O[,'xy']*DE[,'X-x'] + O[,'yy']*DE[,'Y-y'])*DE[,'Y-y']
    sum(log(detE) + maha)/2
}

LL.titterington <- function(abABlw,XYZ,wtype='a'){
    x <- get.titterington.xyz(XYZ=XYZ,a=abABlw['a'],b=abABlw['b'],
                              A=abABlw['A'],B=abABlw['B'],
                              w=exp(abABlw['lw']),wtype=wtype)[,'x']
    LL.titterington.abABlwx(c(abABlw,x),XYZ=XYZ,wtype=wtype)
}
LL.titterington.abABlwx <- function(abABlwx,XYZ,wtype='a'){
    ns <- nrow(XYZ)
    a <- abABlwx['a']
    b <- abABlwx['b']
    A <- abABlwx['A']
    B <- abABlwx['B']
    w <- exp(abABlwx['lw'])
    x <- abABlwx[6:(ns+5)]
    DE <- matrix(0,nrow=ns,ncol=9)
    colnames(DE) <- c('X-x','Y-y','Z-z','vX','vY','vZ','sXY','sXZ','sYZ')
    DE[,'X-x'] <- XYZ[,'X']-x
    DE[,'Y-y'] <- XYZ[,'Y']-a-b*x
    DE[,'Z-z'] <- XYZ[,'Z']-A-B*x
    DE[,'vX'] <- XYZ[,'sX']^2
    DE[,'vY'] <- XYZ[,'sY']^2
    DE[,'vZ'] <- XYZ[,'sZ']^2
    if (wtype=='a'){
        DE[,'vY'] <- XYZ[,'sY']^2 + w^2
    } else if (wtype=='b'){
        DE[,'vY'] <- XYZ[,'sY']^2 + (w*x)^2
    } else if (wtype=='A'){
        DE[,'vZ'] <- XYZ[,'sZ']^2 + w^2
    } else if (wtype=='B'){
        DE[,'vZ'] <- XYZ[,'sZ']^2 + (w*x)^2
    }
    DE[,'sXY'] <- XYZ[,'rXY']*XYZ[,'sX']*XYZ[,'sY']
    DE[,'sXZ'] <- XYZ[,'rXZ']*XYZ[,'sX']*XYZ[,'sZ']
    DE[,'sYZ'] <- XYZ[,'rYZ']*XYZ[,'sY']*XYZ[,'sZ']
    detE <- det3x3(vx=DE[,'vX'],vy=DE[,'vY'],vz=DE[,'vZ'],
                   sxy=DE[,'sXY'],sxz=DE[,'sXZ'],syz=DE[,'sYZ'])
    O <- invertcovmat(vx=DE[,'vX'],vy=DE[,'vY'],vz=DE[,'vZ'],
                      sxy=DE[,'sXY'],sxz=DE[,'sXZ'],syz=DE[,'sYZ'])
    maha <-
        (O[,'xx']*DE[,'X-x']+O[,'xy']*DE[,'Y-y']+O[,'xz']*DE[,'Z-z'])*DE[,'X-x']+
        (O[,'xy']*DE[,'X-x']+O[,'yy']*DE[,'Y-y']+O[,'yz']*DE[,'Z-z'])*DE[,'Y-y']+
        (O[,'xz']*DE[,'X-x']+O[,'yz']*DE[,'Y-y']+O[,'zz']*DE[,'Z-z'])*DE[,'Z-z']
    if (FALSE){
        scatterplot(XYZ[,c('X','sX','Y','sY','rXY')])
        points(x,a+b*x)
        title(sum(log(detE) + maha)/2)
    }
    sum(log(detE) + maha)/2
}
