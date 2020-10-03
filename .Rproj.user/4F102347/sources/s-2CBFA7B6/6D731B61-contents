#' UWOT-UMAP
#'
#' @param object
#' @param reduction
#' @param spread
#' @param n_components
#' @param min_dist
#' @param metric
#' @param n_neighbors
#' @param set_op_mix_ratio
#' @param local_connectivity
#' @param repulsion_strength
#' @param negative_sample_rate
#' @param n_threads
#' @param reduction_name
#' @param return_seurat
#' @param verbose whether to print function messages
#'
#' @return
#' @export
#' @import uwot
#' @import Seurat
#'
#' @examples
visUMAP <- function(object,
                    reduction = 'harmony',
                    spread = 1,
                    n_components = 2,
                    min_dist = 0.3,
                    metric = 'cosine',
                    n_neighbors = 30,
                    set_op_mix_ratio = 1,
                    local_connectivity = 1,
                    repulsion_strength = 1,
                    negative_sample_rate = 5,
                    n_threads =  parallel::detectCores()-1,
                    reduction_name = 'umap',
                    return_seurat = TRUE,
                    verbose = TRUE
){
  embds <- Seurat::Embeddings(object, reduction = reduction)
  umap_res <-   uwot::umap(embds,
                           spread = 1,
                           n_components = n_components,
                           min_dist = min_dist,
                           metric = metric,
                           n_threads = n_threads,
                           n_neighbors = n_neighbors,
                           set_op_mix_ratio = set_op_mix_ratio,
                           local_connectivity = local_connectivity,
                           repulsion_strength = repulsion_strength,
                           negative_sample_rate = negative_sample_rate
  )
  if(return_seurat){
    object[reduction_name] <- Seurat::CreateDimReducObject(embeddings = umap_res, key = 'UMAP_', assay = 'RNA')
    return(object)
  } else {
    return(umap_res)
  }

}


#' UWOT-UAMP: Clustering Specific UMAP
#'
#' @param object
#' @param reduction
#' @param spread
#' @param n_components
#' @param min_dist
#' @param metric
#' @param n_neighbors
#' @param set_op_mix_ratio
#' @param local_connectivity
#' @param repulsion_strength
#' @param negative_sample_rate
#' @param n_threads
#' @param reduction_name
#' @param return_seurat
#' @param verbose   whether to print function messages
#'
#' @return
#' @export
#'
#' @examples
clustUMAP <- function(object,
                      reduction = 'harmony',
                      spread = 1.1,
                      n_components = NULL,
                      min_dist = 0,
                      metric = 'cosine',
                      n_neighbors = 50,
                      set_op_mix_ratio = 1,
                      local_connectivity = 1,
                      repulsion_strength = 1,
                      negative_sample_rate = 5,
                      n_threads = parallel::detectCores()-1,
                      reduction_name = 'umap',
                      return_seurat = TRUE,
                      verbose = TRUE
){
  embds <- Seurat::Embeddings(object, reduction = reduction)
  if(is.null(n_components)){
    n_components <- ncol(embds)
  }
  umap_res <-   uwot::umap(embds,
                           spread = spread,
                           n_components = n_components,
                           min_dist = min_dist,
                           n_threads = n_threads,
                           metric = metric,
                           n_neighbors = n_neighbors,
                           set_op_mix_ratio = set_op_mix_ratio,
                           local_connectivity = local_connectivity,
                           repulsion_strength = repulsion_strength,
                           negative_sample_rate = negative_sample_rate,
                           verbose = verbose
  )

  rownames(umap_res) <- rownames(embds)
  colnames(umap_res) <- paste0('UMAP_', 1:ncol(umap_res))

  if(return_seurat){
    object[reduction_name] <- Seurat::CreateDimReducObject(embeddings = umap_res,
                                                           key = 'clustUMAP_',
                                                           assay = 'RNA')
    return(object)
  } else {
    return(umap_res)
  }

}
