FROM alpine:3.20

RUN apk add R R-dev R-doc build-base automake autoconf ttf-freefont zlib-dev

COPY . /isoplotr

RUN Rscript --vanilla -e "install.packages('MASS', repos='https://cran.r-project.org/')"

RUN R CMD INSTALL /isoplotr