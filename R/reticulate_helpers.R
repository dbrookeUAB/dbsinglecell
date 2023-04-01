# global reference to scipy (will be initialized in .onLoad)
scipy <- NULL
hdbscan <- NULL
umap <- NULL
pacmap <- NULL
#' onLoad
#'
#' @param libname
#' @param pkgname
#'
#' @return
#' @export
#'
#' @examples
.onLoad <- function(libname, pkgname) {
  # reticulate::configure_environment("dbsinglecell", force = TRUE)
  reticulate::configure_environment(pkgname, force = TRUE)
  # reticulate::use_condaenv(condaenv = 'r-reticulate')
  # use superassignment to update global reference to scipy
  # test <<- reticulate::py_module_available('umap')
  # scipy <<- reticulate::import('scipy', delay_load = TRUE)
  hdbscan <<- reticulate::import('hdbscan', delay_load = TRUE)
  umap <<- reticulate::import('umap', delay_load = TRUE)
  pacmap <<- reticulate::import('pacmap', delay_load = TRUE)

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

