pkgname <- "IsoplotR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "IsoplotR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('IsoplotR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Pb0corr")
### * Pb0corr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Pb0corr
### Title: Common Pb correction
### Aliases: Pb0corr

### ** Examples

data(examples)
UPb <- Pb0corr(examples$UPb,option=1)
concordia(UPb)
# produces identical results as:
dev.new()
concordia(examples$UPb,common.Pb=1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Pb0corr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("age")
### * age

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: age
### Title: Calculate isotopic ages
### Aliases: age age.default age.UPb age.PbPb age.ArAr age.KCa age.UThHe
###   age.fissiontracks age.ThU age.ReOs age.SmNd age.RbSr age.LuHf

### ** Examples

data(examples)
tUPb <- age(examples$UPb,type=1)
tconc <- age(examples$UPb,type=2)
tdisc <- age(examples$UPb,type=3)
tArAr <- age(examples$ArAr)
tiso <- age(examples$ArAr,isochron=TRUE,i2i=TRUE)
tcentral <- age(examples$FT1,central=TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("age", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("agespectrum")
### * agespectrum

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: agespectrum
### Title: Plot a (40Ar/39Ar) release spectrum
### Aliases: agespectrum agespectrum.default agespectrum.ArAr

### ** Examples

data(examples)
agespectrum(examples$ArAr,ylim=c(0,80))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("agespectrum", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cad")
### * cad

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cad
### Title: Plot continuous data as cumulative age distributions
### Aliases: cad cad.default cad.detritals cad.UPb cad.PbPb cad.ArAr
###   cad.KCa cad.ThU cad.ReOs cad.SmNd cad.RbSr cad.LuHf cad.UThHe
###   cad.fissiontracks

### ** Examples

data(examples)
cad(examples$DZ,verticals=FALSE,pch=20)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cad", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("central")
### * central

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: central
### Title: Calculate U-Th-He and fission track central ages and
###   compositions
### Aliases: central central.default central.UThHe central.fissiontracks

### ** Examples

data(examples)
print(central(examples$UThHe)$age)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("central", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("concordia")
### * concordia

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: concordia
### Title: Concordia diagram
### Aliases: concordia

### ** Examples

data(examples)
concordia(examples$UPb,show.age=2)

dev.new()
concordia(examples$UPb,wetherill=FALSE,
          xlim=c(24.9,25.4),ylim=c(0.0508,0.0518),
          ticks=249:254,exterr=TRUE)

dev.new()
concordia(examples$UPb,wetherill=FALSE,show.age=2,anchor=list(TRUE,0))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("concordia", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("data2york")
### * data2york

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: data2york
### Title: Prepare geochronological data for York regression
### Aliases: data2york data2york.default data2york.UPb data2york.ArAr
###   data2york.KCa data2york.PbPb data2york.PD data2york.UThHe
###   data2york.ThU

### ** Examples

f <- system.file("RbSr1.csv",package="IsoplotR")
dat <- read.csv(f)
yorkdat <- data2york(dat)
fit <- york(yorkdat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("data2york", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ellipse")
### * ellipse

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ellipse
### Title: Get error ellipse coordinates for plotting
### Aliases: ellipse

### ** Examples

x = 99; y = 101;
covmat <- matrix(c(1,0.9,0.9,1),nrow=2)
ell <- ellipse(x,y,covmat)
plot(c(90,110),c(90,110),type='l')
polygon(ell,col=rgb(0,1,0,0.5))
points(x,y,pch=21,bg='black')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ellipse", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("evolution")
### * evolution

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: evolution
### Title: Th-U evolution diagram
### Aliases: evolution

### ** Examples

data(examples)
evolution(examples$ThU)

dev.new()
evolution(examples$ThU,transform=TRUE,
          isochron=TRUE,model=1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("evolution", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("examples")
### * examples

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: examples
### Title: Example datasets for testing 'IsoplotR'
### Aliases: examples

### ** Examples

data(examples)

concordia(examples$UPb)

agespectrum(examples$ArAr)

isochron(examples$ReOs)

radialplot(examples$FT1)

helioplot(examples$UThHe)

evolution(examples$ThU)

kde(examples$DZ)

radialplot(examples$LudwigMixture)

agespectrum(examples$LudwigSpectrum)

weightedmean(examples$LudwigMean)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("examples", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("helioplot")
### * helioplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: helioplot
### Title: Visualise U-Th-He data on a logratio plot or ternary diagram
### Aliases: helioplot

### ** Examples

data(examples)
helioplot(examples$UThHe)
dev.new()
helioplot(examples$UThHe,logratio=FALSE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("helioplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("isochron")
### * isochron

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: isochron
### Title: Calculate and plot isochrons
### Aliases: isochron isochron.default isochron.PbPb isochron.ArAr
###   isochron.KCa isochron.RbSr isochron.ReOs isochron.SmNd isochron.LuHf
###   isochron.ThU isochron.UThHe

### ** Examples

data(examples)
isochron(examples$RbSr)

fit <- isochron(examples$ArAr,inverse=FALSE,plot=FALSE)

dev.new()
isochron(examples$ThU,type=4)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("isochron", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("kde")
### * kde

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: kde
### Title: Create (a) kernel density estimate(s)
### Aliases: kde kde.default kde.UPb kde.detritals kde.PbPb kde.ArAr
###   kde.KCa kde.ThU kde.ReOs kde.SmNd kde.RbSr kde.LuHf kde.UThHe
###   kde.fissiontracks

### ** Examples

kde(examples$UPb)

dev.new()
kde(examples$FT1,log=TRUE)

dev.new()
kde(examples$DZ,from=1,to=3000,kernel="epanechnikov")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("kde", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ludwig")
### * ludwig

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ludwig
### Title: Linear regression of U-Pb data with correlated errors, taking
###   into account decay constant uncertainties.
### Aliases: ludwig ludwig.default ludwig.UPb

### ** Examples

f <- system.file("UPb4.csv",package="IsoplotR")
d <- read.data(f,method="U-Pb",format=4)
fit <- ludwig(d)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ludwig", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mds")
### * mds

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mds
### Title: Multidimensional Scaling
### Aliases: mds mds.default mds.detritals

### ** Examples

data(examples)
mds(examples$DZ,nnlines=TRUE,pch=21,cex=5)
dev.new()
mds(examples$DZ,shepard=TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mds", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("peakfit")
### * peakfit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: peakfit
### Title: Finite mixture modelling of geochronological datasets
### Aliases: peakfit peakfit.default peakfit.fissiontracks peakfit.UPb
###   peakfit.PbPb peakfit.ArAr peakfit.KCa peakfit.ReOs peakfit.SmNd
###   peakfit.RbSr peakfit.LuHf peakfit.ThU peakfit.UThHe

### ** Examples

data(examples)
peakfit(examples$FT1,k=2)

peakfit(examples$LudwigMixture,k='min')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("peakfit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("radialplot")
### * radialplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: radialplot
### Title: Visualise heteroscedastic data on a radial plot
### Aliases: radialplot radialplot.default radialplot.fissiontracks
###   radialplot.UPb radialplot.PbPb radialplot.ArAr radialplot.KCa
###   radialplot.UThHe radialplot.ReOs radialplot.SmNd radialplot.RbSr
###   radialplot.LuHf radialplot.ThU

### ** Examples

data(examples)
radialplot(examples$FT1)

dev.new()
radialplot(examples$LudwigMixture,k='min')




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("radialplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read.data")
### * read.data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read.data
### Title: Read geochronology data
### Aliases: read.data read.data.default read.data.data.frame
###   read.data.matrix

### ** Examples


f1 <- system.file("UPb1.csv",package="IsoplotR")
file.show(f1) # inspect the contents of 'UPb1.csv'
d1 <- read.data(f1,method="U-Pb",format=1)
concordia(d1)

f2 <- system.file("ArAr1.csv",package="IsoplotR")
d2 <- read.data(f2,method="Ar-Ar",format=1)
agespectrum(d2)

f3 <- system.file("ReOs1.csv",package="IsoplotR")
d3 <- read.data(f3,method="Re-Os",format=1)
isochron(d2)

f4 <- system.file("FT1.csv",package="IsoplotR")
d4 <- read.data(f4,method="fissiontracks",format=1)
radialplot(d4)

f5 <- system.file("UThSmHe.csv",package="IsoplotR")
d5 <- read.data(f5,method="U-Th-He")
helioplot(d5)

f6 <- system.file("ThU2.csv",package="IsoplotR")
d6 <- read.data(f6,method="Th-U",format=2)
evolution(d6)

#  one detrital zircon U-Pb file (detritals.csv)
f7 <- system.file("DZ.csv",package="IsoplotR")
d7 <- read.data(f7,method="detritals")
kde(d7)

#  four 'other' files (LudwigMixture.csv, LudwigSpectrum.csv,
#  LudwigMean.csv, LudwigKDE.csv)
f8 <- system.file("LudwigMixture.csv",package="IsoplotR")
d8 <- read.data(f8,method="other")
radialplot(d8)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read.data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("scatterplot")
### * scatterplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: scatterplot
### Title: Create a scatter plot with error ellipses or crosses
### Aliases: scatterplot

### ** Examples

X <- c(1.550,12.395,20.445,20.435,20.610,24.900,
       28.530,50.540,51.595,86.51,106.40,157.35)
Y <- c(.7268,.7809,.8200,.8116,.8160,.8302,
       .8642,.9534,.9617,1.105,1.230,1.440)
sX <- X*0.02
sY <- Y*0.01
dat <- cbind(X,sX,Y,sY)
scatterplot(dat,fit=york(dat),show.ellipses=2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("scatterplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("set.zeta")
### * set.zeta

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: set.zeta
### Title: Calculate the zeta calibration coefficient for fission track
###   dating
### Aliases: set.zeta

### ** Examples

data(examples)
print(examples$FT1$zeta)
FT <- set.zeta(examples$FT1,tst=c(250,5))
print(FT$zeta)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("set.zeta", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("settings")
### * settings

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: settings
### Title: Load settings to and from json
### Aliases: settings

### ** Examples

# load and show the default constants that come with IsoplotR
json <- system.file("constants.json",package="IsoplotR")
settings(fname=json)
print(settings())

# use the decay constant of Kovarik and Adams (1932)
settings('lambda','U238',0.0001537,0.0000068)
print(settings('lambda','U238'))

# returns the 238U/235U ratio of Hiess et al. (2012):
print(settings('iratio','U238U235'))
# use the 238U/235U ratio of Steiger and Jaeger (1977):
settings('iratio','U238U235',138.88,0)
print(settings('iratio','U238U235'))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("settings", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("titterington")
### * titterington

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: titterington
### Title: Linear regression of X,Y,Z-variables with correlated errors
### Aliases: titterington

### ** Examples

d <- matrix(c(0.1677,0.0047,1.105,0.014,0.782,0.015,0.24,0.51,0.33,
              0.2820,0.0064,1.081,0.013,0.798,0.015,0.26,0.63,0.32,
              0.3699,0.0076,1.038,0.011,0.819,0.015,0.27,0.69,0.30,
              0.4473,0.0087,1.051,0.011,0.812,0.015,0.27,0.73,0.30,
              0.5065,0.0095,1.049,0.010,0.842,0.015,0.27,0.76,0.29,
              0.5520,0.0100,1.039,0.010,0.862,0.015,0.27,0.78,0.28),
            nrow=6,ncol=9)
colnames(d) <- c('X','sX','Y','sY','Z','sZ','rXY','rXZ','rYZ')
titterington(d)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("titterington", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("weightedmean")
### * weightedmean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: weightedmean
### Title: Calculate the weighted mean age
### Aliases: weightedmean weightedmean.default weightedmean.UPb
###   weightedmean.PbPb weightedmean.ThU weightedmean.ArAr weightedmean.KCa
###   weightedmean.ReOs weightedmean.SmNd weightedmean.RbSr
###   weightedmean.LuHf weightedmean.UThHe weightedmean.fissiontracks

### ** Examples

ages <- c(251.9,251.59,251.47,251.35,251.1,251.04,250.79,250.73,251.22,228.43)
errs <- c(0.28,0.28,0.63,0.34,0.28,0.63,0.28,0.4,0.28,0.33)
weightedmean(cbind(ages,errs))

data(examples)
weightedmean(examples$LudwigMean)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("weightedmean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("york")
### * york

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: york
### Title: Linear regression of X,Y-variables with correlated errors
### Aliases: york

### ** Examples

X <- c(1.550,12.395,20.445,20.435,20.610,24.900,
       28.530,50.540,51.595,86.51,106.40,157.35)
Y <- c(.7268,.7849,.8200,.8156,.8160,.8322,
       .8642,.9584,.9617,1.135,1.230,1.490)
n <- length(X)
sX <- X*0.01
sY <- Y*0.005
rXY <- rep(0.8,n)
dat <- cbind(X,sX,Y,sY,rXY)
fit <- york(dat)
scatterplot(dat,fit=fit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("york", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
