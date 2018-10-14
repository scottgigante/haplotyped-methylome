chooseCRANmirror(ind = 1)
source("https://bioconductor.org/biocLite.R")

install_CRAN <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
}

install_Bioconductor <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) biocLite(pkg, suppressUpdates=TRUE, ask=FALSE)
}

for (pkg in readLines("r_requirements.txt")) {
  install_CRAN(pkg)
}

for (pkg in readLines("bioc_requirements.txt")) {
  install_Bioconductor(pkg)
}
