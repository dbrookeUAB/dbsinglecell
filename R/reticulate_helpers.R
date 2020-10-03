# global reference to scipy (will be initialized in .onLoad)
scipy <- NULL
hdbscan <- NULL
umap <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  scipy <<- reticulate::import("scipy", delay_load = TRUE)
  hdbscan <<- reticulate::import('hdbscan', delay_load = TRUE)
  umap <<- reticulate::import('umap', delay_load = TRUE)
}

install_python_packages <- function(method = "auto", conda = "auto") {
  reticulate::py_install(c("hdscan",'umap'), method = method, conda = conda)
}
