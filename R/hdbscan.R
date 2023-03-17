#'  Hierarchical Density-Based Spatial Clustering of Applications with Noise
#'
#' @param x
#' @param algorithm
#' @param alpha
#' @param approx_min_span_tree
#' @param gen_min_span_tree
#' @param leaf_size
#' @param metric
#' @param min_cluster_size
#' @param min_samples
#' @param cluster_selection_epsilon
#' @param cluster_selection_method
#' @param nThreads
#' @param prediction_data  not sure what this is for. Will update later.
#'
#' @return
#' @export
#'
#' @import reticulate
#'
#' @examples
HDBSCAN <- function(object, ...){
  UseMethod(generic = 'HDBSCAN', object = object)
}

#' Hierarchical Density-Based Spatial Clustering of Applications with Noise
#'
#' @param x
#' @param algorithm
#' @param alpha
#' @param approx_min_span_tree
#' @param gen_min_span_tree
#' @param leaf_size
#' @param metric
#' @param prediction_data
#' @param min_cluster_size
#' @param min_samples
#' @param cluster_selection_epsilon
#' @param cluster_selection_method
#' @param nThreads
#'
#' @return
#' @export
#'
#' @examples
HDBSCAN.default <- function(x,
                           algorithm='best',
                           alpha=1.0,
                           approx_min_span_tree = TRUE,
                           gen_min_span_tree=FALSE,
                           leaf_size=30,
                           metric='euclidean',
                           prediction_data=TRUE,
                           min_cluster_size =1000L,
                           min_samples = 1L,
                           cluster_selection_epsilon = 0.1,
                           cluster_selection_method = 'leaf',
                           nThreads = parallel::detectCores()-1,
                           return_full = FALSE
){

  hdbscan <- reticulate::import('hdbscan', delay_load = TRUE)



  clusterer <- hdbscan$HDBSCAN(algorithm = algorithm,
                               allow_single_cluster = FALSE,
                               alpha = alpha,
                               prediction_data = prediction_data,
                               approx_min_span_tree = approx_min_span_tree,
                               gen_min_span_tree = gen_min_span_tree,
                               leaf_size = leaf_size,
                               core_dist_n_jobs = nThreads,
                               metric = metric,
                               min_cluster_size = as.integer(min_cluster_size),
                               min_samples = as.integer(min_samples),
                               cluster_selection_epsilon =  cluster_selection_epsilon,
                               cluster_selection_method = cluster_selection_method
  )



  reticulate::py_main_thread_func(clusterer$fit(x))


  result <- list(
      labels = clusterer$labels_,
      probabilities = clusterer$probabilities_,
      cluster_persistance = clusterer$cluster_persistence_,
      exemplars = clusterer$exemplars_,
      outlier_scores = clusterer$outlier_scores_)



  if(return_full==TRUE){
    return(clusterer)
  } else {
    return(result)
  }

}

#' Hierarchical Density-Based Spatial Clustering of Applications with Noise
#'
#' @param x
#' @param algorithm
#' @param alpha
#' @param approx_min_span_tree
#' @param gen_min_span_tree
#' @param leaf_size
#' @param metric
#' @param prediction_data
#' @param min_cluster_size
#' @param min_samples
#' @param cluster_selection_epsilon
#' @param cluster_selection_method
#' @param nThreads
#'
#' @return
#' @export
#'
#' @examples
HDBSCAN.matrix <- function(x,
                    algorithm='best',
                    alpha=1.0,
                    approx_min_span_tree = TRUE,
                    gen_min_span_tree=FALSE,
                    leaf_size=30,
                    metric='euclidean',
                    prediction_data=TRUE,
                    min_cluster_size =1000L,
                    min_samples = 1L,
                    cluster_selection_epsilon = 0.1,
                    cluster_selection_method = 'leaf',
                    nThreads = parallel::detectCores()-1,
                    return_full = FALSE
){

  hdbscan <- reticulate::import('hdbscan', delay_load = TRUE)



  clusterer <- hdbscan$HDBSCAN(algorithm = algorithm,
                               allow_single_cluster = FALSE,
                               alpha = alpha,
                               prediction_data = prediction_data,
                               approx_min_span_tree = approx_min_span_tree,
                               gen_min_span_tree = gen_min_span_tree,
                               leaf_size = leaf_size,
                               core_dist_n_jobs = nThreads,
                               metric = metric,
                               min_cluster_size = as.integer(min_cluster_size),
                               min_samples = as.integer(min_samples),
                               cluster_selection_epsilon =  cluster_selection_epsilon,
                               cluster_selection_method = cluster_selection_method
  )



  reticulate::py_main_thread_func(clusterer$fit(x))


  result <- list(
    labels = clusterer$labels_,
    probabilities = clusterer$probabilities_,
    cluster_persistance = clusterer$cluster_persistence_,
    exemplars = clusterer$exemplars_,
    outlier_scores = clusterer$outlier_scores_)

  if(return_full==TRUE){
    return(clusterer)
  } else {
    return(result)
  }

  reticulate::py_de
  reticulate::py_run_string('del clusterer')
  reticulate::py_gc <- import("gc")
  reticulate::py_gc$collect()
}



#'  Hierarchical Density-Based Spatial Clustering of Applications with Noise
#'
#' @param object
#' @param reduction
#' @param dims
#' @param algorithm
#' @param alpha
#' @param approx_min_span_tree
#' @param gen_min_span_tree
#' @param leaf_size
#' @param metric
#' @param min_cluster_size
#' @param min_samples
#' @param cluster_selection_epsilon
#' @param cluster_selection_method
#' @param nThreads
#' @param return_seurat  logical to return the result within the orignal object or as the raw HDBSCAN result
#' @param prediction_data not sure what this is for. Will update later.
#'
#' @return
#' @export
#'
#' @examples
HDBSCAN.Seurat <- function(object,
                           reduction = 'umap',
                           dims = NULL,
                           algorithm='best',
                           alpha=1.0,
                           prediction_data = TRUE,
                           approx_min_span_tree = TRUE,
                           gen_min_span_tree=FALSE,
                           leaf_size=30,
                           metric='euclidean',
                           min_cluster_size =1000L,
                           min_samples = 1L,
                           cluster_selection_epsilon = 0.1,
                           cluster_selection_method = 'leaf',
                           nThreads = parallel::detectCores()-1,
                           return_seurat = TRUE,
                           return_full = FALSE
){

  if(is.null(dims)){
    x <- Seurat::Embeddings(object, reduction = reduction)
  } else {
    x <- Seurat::Embeddings(object, reduction = reduction)[,dims]
  }

  hdbscan <- reticulate::import('hdbscan', delay_load = TRUE)



  clusterer <- hdbscan$HDBSCAN(algorithm=algorithm,
                               alpha = alpha,
                               prediction_data = prediction_data,
                               approx_min_span_tree = approx_min_span_tree,
                               gen_min_span_tree = gen_min_span_tree,
                               leaf_size = leaf_size,
                               core_dist_n_jobs = nThreads,
                               metric = metric,
                               min_cluster_size = as.integer(min_cluster_size),
                               min_samples = as.integer(min_samples),
                               cluster_selection_epsilon =  cluster_selection_epsilon,
                               cluster_selection_method = cluster_selection_method
  )
  reticulate::py_main_thread_func(clusterer$fit(x))



  result <- list(
    labels = factor(clusterer$labels_),
    probabilities = clusterer$probabilities_,
    cluster_persistance = clusterer$cluster_persistence_,
    exemplars = clusterer$exemplars_,
    outlier_scores = clusterer$outlier_scores_)

  if(return_seurat){
    if(return_full==TRUE){
      object@misc$hdbscan <- clusterer
    } else {
      object@misc$hdbscan <- result
    }

    object$hdbscan_clusters <- result$labels
    object$outlier_scores <- result$outlier_scores
    return(object)
  } else {

    if(return_full==TRUE){
      return(clusterer)
    } else {
      return(result)
    }



  }
reticulate::py_de
reticulate::py_run_string('del clusterer')
reticulate::py_gc <- import("gc")
reticulate::py_gc$collect()
}
