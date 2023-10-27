# dispatch functions for all the regression algorithms (except Ludwig)

regression <- function(xyz,model=1,type='york',omit=NULL,abwtype=1,abanchor=0){
    xyz2calc <- clear(xyz,omit,OGLS=identical(type,'ogls'))
    if (model==1){
        out <- model1regression(xyz2calc,type=type,abanchor=abanchor)
    } else if (model==2){
        out <- model2regression(xyz2calc,type=type,abanchor=abanchor)
    } else if (model==3){
        out <- model3regression(xyz2calc,type=type,
                                abwtype=abwtype,abanchor=abanchor)
        out$abwtype <- abwtype
    } else if (model==4 && identical(type,'york')){
        out <- irr(xyz2calc)
        out$abwtype <- abwtype
    } else if (model==5 && identical(type,'york')){
        out <- irr(xyz2calc,abwtype=abwtype)
        out$abwtype <- abwtype
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

model1regression <- function(xyz,type='york',abanchor=0){
    if (identical(type,'york')){
        if (abanchor[1]<1 | length(abanchor)<2){
            out <- york(xyz)
            out$par <- c(out$a[1],out$b[1],lw=-Inf)
            out$cov <- rbind(c(out$a[2]^2,out$cov.ab,0),
                             c(out$cov.ab,out$b[2]^2,0),
                             c(0,0,0))
            names(out$par) <- colnames(out$cov) <-
                rownames(out$cov) <- c('a','b','lw')
        } else {
            out <- MLyork(xyz,abanchor=abanchor)
        }
    } else if (identical(type,'titterington')){
        out <- titterington(xyz)
    } else if (identical(type,'ogls')){
        out <- ogls(xyz,random.effects=FALSE)
    } else {
        stop('invalid output type for model 1 regression')
    }
    out
}

model2regression <- function(xyz,type='york',abanchor=0){
    if (identical(type,'york')){
        out <- MLyork(xyz[,c('X','Y')],abanchor=abanchor,model=2)
        out$df <- nrow(xyz)-2
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

model3regression <- function(xyz,type='york',abwtype=1,abanchor=0){
    if (identical(type,'york')){
        return(MLyork(yd=xyz,abanchor=abanchor,model=3,abwtype=abwtype))
    } else if (identical(type,'titterington')){
        pilot <- model1regression(xyz,type=type)
        ilw <- init.titterington.lw(XYZ=xyz,abwtype=abwtype,pilot=pilot)$minimum
        init <- c(pilot$par,'lw'=unname(ilw))
        upper <- init + c(5*pilot$par[c('a','b','A','B')],2)
        lower <- init - c(5*pilot$par[c('a','b','A','B')],2)
        out <- contingencyfit(par=init,fn=LL.titterington,lower=lower,
                              upper=upper,XYZ=xyz,abwtype=abwtype)
        out$cov <- E <- inverthess(out$hessian)
    } else if (identical(type,'titterington')){
        stop('not yet implemented')
    } else if (identical(type,'ogls')){
        out <- ogls(xyz,random.effects=TRUE)
        out$cov <- E <- inverthess(out$hessian)
    } else {
        stop('invalid output type for model 3 regression')
    }
    disp <- exp(out$par['lw'])
    sdisp <- disp*sqrt(E['lw','lw'])
    out$disp <- c('w'=unname(disp),'s[w]'=unname(sdisp))
    out
}

init.titterington.lw <- function(XYZ,abwtype=1,pilot){
    fact <- max(1,sqrt(pilot$mswd))
    spar <- fact*sqrt(diag(pilot$cov))
    if (abwtype%in%c('intercept',0,'a')) init <- log(spar['a'])
    else if (abwtype%in%c(1,'b')) init <- log(spar['b'])
    else if (abwtype%in%c(2,'A')) init <- log(spar['A'])
    else if (abwtype%in%c(3,'B')) init <- log(spar['B'])
    else stop('illegal abwtype')
    stats::optimise(f=LL.titterington.lw,interval=init+c(-10,5),
                    abAB=pilot$par,XYZ=XYZ,abwtype=abwtype)
}
LL.titterington.lw <- function(lw,abAB,XYZ,abwtype=1){
    LL.titterington(abABlw=c(abAB,lw=unname(lw)),XYZ=XYZ,abwtype=abwtype)
}
LL.titterington <- function(abABlw,XYZ,wtype=1){
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
    if (wtype==1){
        DE[,'vY'] <- XYZ[,'sY']^2 + w^2
    } else if (wtype==2){
        DE[,'vY'] <- XYZ[,'sY']^2 + (w*x)^2
    } else if (wtype==3){
        DE[,'vZ'] <- XYZ[,'sZ']^2 + w^2
    } else if (wtype==4){
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
