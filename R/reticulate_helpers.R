# global reference to scipy (will be initialized in .onLoad)
scipy <- NULL
hdbscan <- NULL
umap <- NULL
pacmap <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  scipy <<- reticulate::import("scipy", delay_load = TRUE)
  hdbscan <<- reticulate::import('hdbscan', delay_load = TRUE)
  umap <<- reticulate::import('umap', delay_load = TRUE)
  #pacmap <<- reticulate::import('pacmap', delay_load = TRUE)

}

#'  Install python packages to use `HDBSCAN`, `umap-learn`, and `PaCMAP`
#'
#' @param method
#' @param conda
#'
#' @return
#' @export
#'
#' @examples
install_python_packages <- function(method = "auto", conda = "auto") {

  reticulate::conda_install(packages = c("hdbscan",'umap-learn','pacmap'),
                         method = method,
                         conda = conda)
}

