#' Monocle3 workflow for the lazy man
#'
#' @param object Seurat object
#' @param dims dims for preprocess
#' @param verbose Default = TRUE
#'
#' @import monocle3
#'
#' @return
#' @export
#'
#' @examples
monocle_processing <- function(object, dims = 100,verbose = TRUE){
  message_task('Converting Seurat object to cell_data_set')
  cds <- SeuratWrappers::as.cell_data_set(object)

  message_task('preprocess_cds')
  cds <- monocle3::preprocess_cds(cds,
                                  num_dim = dims,
                                  verbose = verbose )

  message_task('align_cds')
  cds <- monocle3::align_cds(cds, verbose = verbose)

  message_task('reduce_dimensions')
  cds <- monocle3::reduce_dimension(cds,reduction_method = 'UMAP', verbose = verbose, umap.min_dist = 0.05, umap.n_neighbors = 10)

  message_task('cluster_cells')
  cds <- monocle3::cluster_cells(cds, reduction_method = 'UMAP', verbose = verbose)

  message_task('learn_graph')
  cds <- monocle3::learn_graph(cds, verbose = verbose)

  message_task('Done!!')
  return(cds)
}
