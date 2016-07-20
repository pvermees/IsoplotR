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
nameEx("age")
### * age

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: age
### Title: Calculate isotopic ages
### Aliases: age age.ArAr age.UPb age.default age.detritals

### ** Examples

data(examples)
print(age(examples$UPb))
print(age(examples$UPb,concordia=1))
print(age(examples$UPb,concordia=2))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("age", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("botev")
### * botev

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: botev
### Title: Compute the optimal kernel bandwidth
### Aliases: botev

### ** Examples

data(examples)
samp <- examples$DZ[['N1']]
bw <- botev(samp)
print(bw)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("botev", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cad")
### * cad

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cad
### Title: Plot continuous data as cumulative age distributions
### Aliases: cad

### ** Examples

data(examples)
cad(examples$DZ)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cad", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
concordia(examples$UPb)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("concordia", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ellipse")
### * ellipse

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ellipse
### Title: Get coordinates of error ellipse for plotting
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
dev.new()
kde(examples$DZ)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("examples", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("iratio")
### * iratio

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: iratio
### Title: Isotopic ratios
### Aliases: iratio

### ** Examples

# returns the 238U/235U ratio of Hiess et al. (2012):
print(iratio('U238U235'))
# use the 238U/235U ratio of Steiger and Jaeger (1977):
iratio('U238U235',138.88,0)
print(iratio('U238U235'))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("iratio", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("isochron")
### * isochron

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: isochron
### Title: Calculate and plot isochrons
### Aliases: isochron isochron.ArAr isochron.default

### ** Examples

data(examples)
isochron(examples$ArAr)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("isochron", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("kde")
### * kde

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: kde
### Title: Create (a) kernel density estimate(s)
### Aliases: kde kde.ArAr kde.UPb kde.default kde.detritals

### ** Examples

data(examples)
kde(examples$DZ[['N1']],kernel="epanechnikov")
kde(examples$DZ,from=0,to=3000)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("kde", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lambda")
### * lambda

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lambda
### Title: Decay constants
### Aliases: lambda

### ** Examples

print(lambda('U238'))
# use the decay constant of Kovarik and Adams (1932)
lambda('U238',0.0001537,0.0000068)
print(lambda('U238'))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lambda", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read.data")
### * read.data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read.data
### Title: Read geochronology data
### Aliases: read.data read.data.default read.data.matrix

### ** Examples

# load one of the built-in .csv files:
data(examples)#fname <- system.file("UPb.csv",package="IsoplotR")
#UPb <- read.data(fname,'U-Pb')
concordia(examples$UPb)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read.data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("settings")
### * settings

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: settings
### Title: Load settings to and from json
### Aliases: settings

### ** Examples

json <- system.file("defaults.json",package="IsoplotR")
settings(json)
print(settings())



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("settings", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("weightedmean")
### * weightedmean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: weightedmean
### Title: Calculate the weighted mean age
### Aliases: weightedmean weightedmean.ArAr weightedmean.UPb
###   weightedmean.default

### ** Examples

data(examples)
weightedmean(examples$ArAr)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("weightedmean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("yorkfit")
### * yorkfit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: yorkfit
### Title: Linear regression of X,Y-variables with correlated errors
### Aliases: yorkfit

### ** Examples

   X <- c(1.550,12.395,20.445,20.435,20.610,24.900,
          28.530,50.540,51.595,86.51,106.40,157.35)
   Y <- c(.7268,.7849,.8200,.8156,.8160,.8322,
          .8642,.9584,.9617,1.135,1.230,1.490)
   n <- length(X)
   sX <- X*0.01
   sY <- Y*0.005
   rXY <- rep(0.8,n)
   fit <- yorkfit(X,Y,sX,sY,rXY)
   covmat <- matrix(0,2,2)
   plot(range(X),fit$a[1]+fit$b[1]*range(X),type='l',ylim=range(Y))
   for (i in 1:n){
       covmat[1,1] <- sX[i]^2
       covmat[2,2] <- sY[i]^2
       covmat[1,2] <- rXY[i]*sX[i]*sY[i]
       covmat[2,1] <- covmat[1,2]
       ell <- ellipse(X[i],Y[i],covmat,alpha=0.05)
       polygon(ell)
   }



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("yorkfit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
