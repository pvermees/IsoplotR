.onLoad <- function(libname, pkgname){
    settings(fname=system.file("constants.json",package=pkgname))
}
