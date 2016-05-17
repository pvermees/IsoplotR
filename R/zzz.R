.onLoad <- function(libname, pkgname){
    settings(system.file("constants.json",package=pkgname))
}
