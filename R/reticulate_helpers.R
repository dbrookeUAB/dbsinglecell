# global reference to scipy (will be initialized in .onLoad)
scipy <- NULL
numpy <- NULL
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

  reticulate::configure_environment(pkgname, force = TRUE)
  # reticulate::use_condaenv(condaenv = 'r-reticulate')
  # use superassignment to update global reference to scipy
  # test <<- reticulate::py_module_available('umap')
  # scipy <<- reticulate::import('scipy', delay_load = TRUE)

  hdbscan <<- reticulate::import('hdbscan', delay_load = TRUE)
  umap <<- reticulate::import('umap', delay_load = TRUE)
  numpy <<- reticulate::import('numpy', delay_load = TRUE)
  # pacmap <<- reticulate::import('pacmap', delay_load = TRUE)

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
install_python_packages <- function() {

  if(!reticulate::virtualenv_exists('dbsinglecell')){
    if(!dir.exists(reticulate::miniconda_path())){
      reticulate::install_miniconda()
    }
    reticulate::install_python('3.9:latest')
    reticulate::virtualenv_create(envname = 'dbsinglecell',
                                  version = '3.9:latest',
                                  packages = c('numpy','umap-learn','hdbscan'))

    reticulate::use_virtualenv('dbsinglecell')
  } else {
    reticulate::install_python('3.9:latest')
    reticulate::use_virtualenv('dbsinglecell')
  }
  reticulate::configure_environment('dbsinglecell', force = TRUE)
}

