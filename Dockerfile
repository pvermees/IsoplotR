FROM alpine:3.11.6

RUN apk add R R-dev R-doc build-base automake autoconf ttf-freefont

COPY . /app

RUN Rscript --vanilla -e \
    "install.packages('MASS', repos='https://cran.rstudio.com/')"

RUN R CMD INSTALL /app