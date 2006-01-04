.First.lib <- function(lib, pkg) {
    v <- read.dcf(file=system.file("DESCRIPTION", package="scuba"),
                  fields="Version")
    cat(paste("\nscuba", v, "\n"))
    cat("Type \"help(scuba)\" for an introduction\n")
    cat("Read the warnings in \"help(scuba.disclaimer)\" \n")
}
