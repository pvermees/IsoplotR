FROM alpine:3.19.1

RUN apk add R R-dev R-doc build-base automake autoconf ttf-freefont

COPY . /isoplotr

RUN Rscript --vanilla -e "install.packages(pkgs='https://cran.r-project.org/package=MASS&version=7.3-60.0.1', repos=NULL)"

RUN R CMD INSTALL /isoplotr