.First.lib <- function(lib, pkg)
{
    library.dynam("spc", pkg, lib)
    cat("spc package version 0.4.0\n")  
}

