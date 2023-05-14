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
    if (identical(type,'york')){
        ilw <- init.york.lw(XY=xyz,wtype=wtype,pilot=pilot)$minimum
        init <- c(pilot$a['a'],pilot$b['b'],lw=ilw)
        upper <- init + c(5*pilot$a['s[a]'],5*pilot$b['s[b]'],2)
        lower <- init - c(5*pilot$a['s[a]'],5*pilot$b['s[b]'],2)
        out <- contingencyfit(par=init,fn=LL.york,lower=lower,
                              upper=upper,XY=xyz,wtype=wtype)
        out$a <- c(out$par['a'],'s[a]'=NA)
        out$b <- c(out$par['b'],'s[b]'=NA)
        E <- inverthess(out$hessian)
        out$a['s[a]'] <- unname(sqrt(E['a','a']))
        out$b['s[b]'] <- unname(sqrt(E['b','b']))
        out$cov.ab <- E['a','b']
    } else if (identical(type,'titterington')){
        ilw <- init.titterington.lw(XYZ=xyz,wtype=wtype,pilot=pilot)$minimum
        init <- c(pilot$par,'lw'=unname(ilw))
        upper <- init + c(5*pilot$par[c('a','b','A','B')],2)
        lower <- init - c(5*pilot$par[c('a','b','A','B')],2)
        out <- contingencyfit(par=init,fn=LL.titterington,lower=lower,
                              upper=upper,XYZ=xyz,wtype=wtype)
        out$cov <- E <- inverthess(out$hessian)
    } else {
        stop('invalid output type for model 3 regression')
    }
    disp <- exp(out$par['lw'])
    sdisp <- disp*sqrt(E['lw','lw'])
    out$disp <- c('w'=unname(disp),'s[w]'=unname(sdisp))
    out
}

init.york.lw <- function(XY,wtype='a',pilot){
    err <- ifelse(wtype%in%c('intercept',0,'a'),
                  pilot$a['s[a]'],pilot$b['s[b]'])
    init <- log(sqrt(pilot$mswd)*err)
    stats::optimise(f=LL.york.lw,interval=init+c(-10,5),
                    ab=c(pilot$a['a'],pilot$b['b']),XY=XY,wtype=wtype)
}
LL.york.lw <- function(lw,ab,XY,wtype='a'){
    LL.york(ablw=c(ab,lw=unname(lw)),XY=XY,wtype=wtype)
}
LL.york.ab <- function(ab,lw,XY,wtype='a'){
    LL.york(ablw=c(ab,lw=unname(lw)),XY=XY,wtype=wtype)
}
LL.york <- function(ablw,XY,wtype='a',debug=FALSE){
    if (debug) browser()
    ns <- nrow(XY)
    a <- ablw['a']
    b <- ablw['b']
    w <- exp(ablw['lw'])
    x <- get.york.xy(XY=XY,a=a,b=b,w=w,wtype=wtype)[,'x']
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

init.titterington.lw <- function(XYZ,wtype='a',pilot){
    fact <- max(1,sqrt(pilot$mswd))
    spar <- fact*sqrt(diag(pilot$cov))
    if (wtype%in%c('intercept',0,'a')) init <- log(spar['a'])
    else if (wtype%in%c(1,'b')) init <- log(spar['b'])
    else if (wtype%in%c(2,'A')) init <- log(spar['A'])
    else if (wtype%in%c(3,'B')) init <- log(spar['B'])
    else stop('illegal wtype')
    stats::optimise(f=LL.titterington.lw,interval=init+c(-10,5),
                    abAB=pilot$par,XYZ=XYZ,wtype=wtype)
}
LL.titterington.lw <- function(lw,abAB,XYZ,wtype='a'){
    LL.titterington(abABlw=c(abAB,lw=unname(lw)),XYZ=XYZ,wtype=wtype)
}
LL.titterington <- function(abABlw,XYZ,wtype='a'){
    ns <- nrow(XYZ)
    a <- abABlw['a']
    b <- abABlw['b']
    A <- abABlw['A']
    B <- abABlw['B']
    w <- exp(abABlw['lw'])
    x <- get.titterington.xyz(XYZ=XYZ,a=a,b=b,A=A,B=B,w=w,wtype=wtype)[,'x']
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
