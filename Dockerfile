FROM alpine:3.18.0

RUN apk add R R-dev R-doc build-base automake autoconf ttf-freefont

COPY . /isoplotr

RUN Rscript --vanilla -e \
    "install.packages('MASS', repos='https://cran.r-project.org/')"

RUN R CMD INSTALL /isoplotr