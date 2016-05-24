roundit <- function(age,err){
    out <- list()
    out$err <- signif(err,2)
    nd <- log10(trunc(abs(age)/err))+2
    out$x <- signif(age,nd)
    out
}
