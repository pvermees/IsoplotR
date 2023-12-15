regression <- function(xyz,model=1,type='york',omit=NULL,wtype='a'){
    xyz2calc <- clear(xyz,omit,OGLS=identical(type,'ogls'))
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
    if (identical(type,'ogls')) out$n <- out$n/2
    out$omit <- omit
    out
}

model1regression <- function(xyz,type='york'){
    if (identical(type,'york')){
        out <- york(xyz)
    } else if (identical(type,'titterington')){
        out <- titterington(xyz)
    } else if (identical(type,'ogls')){
        out <- ogls(xyz,random.effects=FALSE)
    } else {
        stop('invalid output type for model 1 regression')
    }
    out
}

model2regression <- function(xyz,type='york'){
    if (identical(type,'york')){
        out <- MLyork(xyz,model=2)
    } else if (identical(type,'titterington')){
        out <- tls(xyz[,c('X','Y','Z')])
        out$df <- 2*nrow(xyz)-4
    } else if (identical(type,'ogls')){
        yd <- data2york(xyz,format=6)
        out <- model2regression(yd,type='york')
    } else {
        stop('invalid output type for model 2 regression')
    }
    out
}

model3regression <- function(xyz,type='york',wtype='a'){
    if (identical(type,'york')){
        return(MLyork(xyz,model=3,wtype=wtype))
    } else if (identical(type,'titterington')){
        pilot <- model1regression(xyz,type=type)
        ilw <- init.titterington.lw(XYZ=xyz,wtype=wtype,pilot=pilot)$minimum
        init <- c(pilot$par,'lw'=unname(ilw))
        upper <- init + c(5*pilot$par[c('a','b','A','B')],2)
        lower <- init - c(5*pilot$par[c('a','b','A','B')],2)
        out <- contingencyfit(par=init,fn=LL.titterington,lower=lower,
                              upper=upper,XYZ=xyz,wtype=wtype)
        out$cov <- E <- inverthess(out$hessian)
    } else if (identical(type,'ogls')){
        out <- ogls(xyz,random.effects=TRUE)
        out$cov <- E <- inverthess(out$hessian)
    } else {
        stop('invalid output type for model 3 regression')
    }
    w <- exp(out$par['lw'])
    sw <- w*sqrt(E['lw','lw'])
    out$w <- c('w'=unname(w),'s[w]'=unname(sw))
    out
}

init.titterington.lw <- function(XYZ,wtype='a',pilot){
    fact <- max(1,sqrt(pilot$mswd))
    spar <- fact*sqrt(diag(pilot$cov))
    if (wtype%in%c('intercept',1,'a')) init <- log(spar['a'])
    else if (wtype%in%c(2,'b')) init <- log(spar['b'])
    else if (wtype%in%c(3,'A')) init <- log(spar['A'])
    else if (wtype%in%c(4,'B')) init <- log(spar['B'])
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
    sum(log(detE) + maha)/2
}
