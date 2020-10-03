#' UMAP (Uniform Manifold Approximation and Projection for Dimension Reduction)
#'
#' @param embedding
#' @param a
#' @param angular_rp_forest
#' @param b
#' @param force_approximation_algorithm
#' @param init
#' @param learning_rate
#' @param local_connectivity
#' @param low_memory
#' @param metric
#' @param metric_kwds
#' @param min_dist
#' @param n_components
#' @param n_epochs
#' @param n_neighbors
#' @param negative_sample_rate
#' @param output_metric
#' @param output_metric_kwds
#' @param random_state
#' @param repulsion_strength
#' @param set_op_mix_ratio
#' @param spread
#' @param target_metric
#' @param target_metric_kwds
#' @param target_n_neighbors
#' @param target_weight
#' @param transform_queue_size
#' @param transform_seed
#' @param unique
#' @param verbose
#' @param nThreads  number of parallel threads to be used
#'
#' @return
#' @export
#'
#' @import reticulate
#' @import Seurat
#'
#' @examples
umap <- function(
  embedding,
  a=NULL,
  angular_rp_forest=FALSE,
  b=NULL,
  force_approximation_algorithm=FALSE,
  init='spectral',
  learning_rate=1.0,
  local_connectivity=1.0,
  low_memory=FALSE,
  metric='euclidean',
  metric_kwds=NULL,
  min_dist=0.1,
  n_components=2,
  n_epochs=200,
  n_neighbors=15,
  negative_sample_rate=5,
  output_metric='euclidean',
  output_metric_kwds=NULL,
  random_state=42,
  repulsion_strength=1.0,
  set_op_mix_ratio=1.0,
  spread=1.0,
  target_metric='categorical',
  target_metric_kwds=NULL,
  target_n_neighbors=-1,
  target_weight=0.5,
  transform_queue_size=4.0,
  transform_seed=42,
  unique=FALSE,
  verbose=TRUE,
  nThreads = parallel::detectCores()-1
){
  Sys.setenv(OMP_NUM_THREADS=nThreads)
  umap <- reticulate::import('umap', delay_load = TRUE)
  reducer <- umap$UMAP(
    a=a,
    angular_rp_forest=angular_rp_forest,
    b=b,
    force_approximation_algorithm=force_approximation_algorithm,
    init=init,
    learning_rate=learning_rate,
    local_connectivity=as.intger(local_connectivity),
    low_memory=low_memory,
    metric=metric,
    metric_kwds=metric_kwds,
    min_dist=min_dist,
    n_components=as.integer(n_components),
    n_epochs=as.integer(n_epochs),
    n_neighbors=as.integer(n_neighbors),
    negative_sample_rate=negative_sample_rate,
    output_metric=output_metric,
    output_metric_kwds=output_metric_kwds,
    random_state=as.integer(random_state),
    repulsion_strength=repulsion_strength,
    set_op_mix_ratio=set_op_mix_ratio,
    spread=spread,
    target_metric=target_metric,
    target_metric_kwds=target_metric_kwds,
    target_n_neighbors=as.integer(target_n_neighbors),
    target_weight=target_weight,
    transform_queue_size=transform_queue_size,
    transform_seed=as.integer(transform_seed),
    unique=unique,
    verbose=verbose)

  result <- reducer$fit_transform(embedding)

  return(result)
}

umap.Seurat <- function(
  object,
  reduction = 'pca',
  reduction_name = 'umap',
  dims = NULL,
  a=1.662,
  angular_rp_forest=FALSE,
  b=0.7905,
  force_approximation_algorithm=FALSE,
  init='spectral',
  learning_rate=1.0,
  local_connectivity=1.0,
  low_memory=FALSE,
  metric='euclidean',
  metric_kwds=NULL,
  min_dist=0.1,
  n_components=2,
  n_epochs=100,
  n_neighbors=50,
  negative_sample_rate=5,
  output_metric='euclidean',
  output_metric_kwds=NULL,
  random_state=42,
  repulsion_strength=1.0,
  set_op_mix_ratio=1.0,
  spread=1.0,
  target_metric='categorical',
  target_metric_kwds=NULL,
  target_n_neighbors=-1,
  target_weight=0.5,
  transform_queue_size=4.0,
  transform_seed=42,
  unique=FALSE,
  verbose=TRUE,
  nThreads = parallel::detectCores()-1,
  return_seurat = TRUE
){

  if(is.null(dims)){
    embedding <- Seurat::Embeddings(object, reduction = reduction)
  } else {
    embedding <- Seurat::Embeddings(object, reduction = reduction, dims = dims)
  }

  Sys.setenv(OMP_NUM_THREADS=nThreads)
  umap <- reticulate::import('umap', delay_load = TRUE)
  reducer <- umap$UMAP(
    a=a,
    angular_rp_forest=angular_rp_forest,
    b=b,
    force_approximation_algorithm=force_approximation_algorithm,
    init=init,
    learning_rate=learning_rate,
    local_connectivity=local_connectivity,
    low_memory=low_memory,
    metric=metric,
    metric_kwds=metric_kwds,
    min_dist=min_dist,
    n_components=as.integer(n_components),
    n_epochs=as.integer(n_epochs),
    n_neighbors=as.integer(n_neighbors),
    negative_sample_rate=negative_sample_rate,
    output_metric=output_metric,
    output_metric_kwds=output_metric_kwds,
    random_state=as.integer(random_state),
    repulsion_strength=repulsion_strength,
    set_op_mix_ratio=set_op_mix_ratio,
    spread=spread,
    target_metric=target_metric,
    target_metric_kwds=target_metric_kwds,
    target_n_neighbors=target_n_neighbors,
    target_weight=target_weight,
    transform_queue_size=transform_queue_size,
    transform_seed=as.integer(transform_seed),
    unique=unique,
    verbose=verbose)

  result <- reducer$fit_transform(embedding)

  if(return_seurat){
    object[[reduction_name]] <- Seurat::CreateDimReducObject(embeddings = result, key = 'umap_', assay = 'RNA')
    return(object)
  } else {
    return(result)
  }
}



