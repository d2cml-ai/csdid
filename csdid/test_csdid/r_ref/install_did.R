lib <- "C:/Users/h_min/R/R-4.4.2/library"
.libPaths(lib)
# PPM serves prebuilt Windows binaries for R 4.4 (incl. DRDID 1.3.0 and did).
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
options(pkgType = "binary")
options(install.packages.check.source = "no")
install.packages("DRDID")
cat("DRDID:", as.character(packageVersion("DRDID")), "\n")
install.packages("did")
ok <- requireNamespace("did", quietly = TRUE)
cat("DID_INSTALLED:", ok, "\n")
if (ok) cat("did version:", as.character(packageVersion("did")), "\n")