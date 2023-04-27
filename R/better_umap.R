#' umap-learn (python implementation of umap)
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
#' @param nThreads
#'
#' @return
#' @export
#'
#' @import reticulate
#' @import Seurat
#'
#' @examples

umapl <- function(object, ...) {
  UseMethod(generic = 'umapl', object = object)
}

#' umap-learn (python implementation of umap)
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
#' @param nThreads
#'
#' @return
#' @export
#'
#' @examples
umapl.matrix <- function(embedding,
                         a = NULL,
                         angular_rp_forest = FALSE,
                         b = NULL,
                         densmap = FALSE,
                         dens_frac = 0.3,
                         dens_var_shift = 0.1,
                         output_dens = FALSE,
                         disconnection_distance = NULL,
                         dens_lambda = 2.0,
                         force_approximation_algorithm = FALSE,
                         init = 'spectral',
                         learning_rate = 1.0,
                         local_connectivity = 1.0,
                         low_memory = FALSE,
                         metric = 'euclidean',
                         metric_kwds = NULL,
                         min_dist = 0.1,
                         n_components = 2,
                         n_epochs = 200,
                         n_neighbors = 15,
                         negative_sample_rate = 5,
                         output_metric = 'euclidean',
                         output_metric_kwds = NULL,
                         random_state = 42,
                         repulsion_strength = 1.0,
                         set_op_mix_ratio = 1.0,
                         spread = 1.0,
                         target_metric = 'categorical',
                         target_metric_kwds = NULL,
                         target_n_neighbors = -1,
                         target_weight = 0.5,
                         transform_queue_size = 4.0,
                         transform_seed = 42,
                         unique = FALSE,
                         verbose = TRUE,
                         nThreads = 3) {


  UMAP <- reticulate::import('umap', delay_load = TRUE)
  Sys.setenv(OMP_NUM_THREADS = nThreads)
  embds <- embedding
  reticulate::py_capture_output({
    reducer <- UMAP$UMAP(
      a = a,
      angular_rp_forest = angular_rp_forest,
      b = b,
      densmap = densmap,
      dens_frac = dens_frac,
      dens_var_shift = dens_var_shift,
      output_dens = output_dens,
      disconnection_distance = disconnection_distance,
      force_approximation_algorithm = force_approximation_algorithm,
      init = init,
      learning_rate = learning_rate,
      local_connectivity = local_connectivity,
      low_memory = low_memory,
      metric = metric,
      metric_kwds = metric_kwds,
      min_dist = min_dist,
      n_components = as.integer(n_components),
      n_epochs = as.integer(n_epochs),
      n_neighbors = as.integer(n_neighbors),
      negative_sample_rate = negative_sample_rate,
      output_metric = output_metric,
      output_metric_kwds = output_metric_kwds,
      random_state = as.integer(random_state),
      repulsion_strength = repulsion_strength,
      set_op_mix_ratio = set_op_mix_ratio,
      spread = spread,
      target_metric = target_metric,
      target_metric_kwds = target_metric_kwds,
      target_n_neighbors = target_n_neighbors,
      target_weight = target_weight,
      transform_queue_size = transform_queue_size,
      transform_seed = transform_seed,
      unique = unique,
      verbose = verbose
    )
  })


  result <- reducer$fit_transform({
    embedding
  })

  colnames(result) <- paste0('umap_', 1:ncol(result))
  rownames(result) <- rownames(embedding)
  return(result)


}

#' umap-learn (python implementation of umap)
#'
#' @param object
#' @param reduction
#' @param reduction_name
#' @param dims
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
#' @param nThreads
#' @param return_seurat
#'
#' @return
#' @export
#'
#' @examples
umapl.Seurat <- function(object,
                         reduction = 'pca',
                         reduction_name = 'umap',
                         clustering = FALSE,
                         dims = NULL,
                         a = 1.662,
                         angular_rp_forest = FALSE,
                         b = 0.7905,
                         densmap = FALSE,
                         dens_frac = 0.3,
                         dens_var_shift = 0.1,
                         output_dens = FALSE,
                         disconnection_distance = NULL,
                         dens_lambda = 2.0,
                         force_approximation_algorithm = FALSE,
                         init = 'spectral',
                         learning_rate = 1.0,
                         local_connectivity = 1.0,
                         low_memory = FALSE,
                         metric = 'euclidean',
                         metric_kwds = NULL,
                         min_dist = 0.1,
                         n_components = 2,
                         n_epochs = 100,
                         n_neighbors = 50,
                         negative_sample_rate = 5,
                         output_metric = 'euclidean',
                         output_metric_kwds = NULL,
                         random_state = 42,
                         repulsion_strength = 1.0,
                         set_op_mix_ratio = 1.0,
                         spread = 1.0,
                         target_metric = 'categorical',
                         target_metric_kwds = NULL,
                         target_n_neighbors = -1,
                         target_weight = 0.5,
                         transform_queue_size = 4.0,
                         transform_seed = 42,
                         unique = FALSE,
                         verbose = TRUE,
                         nThreads = parallel::detectCores() - 1,
                         return_seurat = TRUE) {
  message_task('importing umap')
  UMAP <- reticulate::import('umap', delay_load = TRUE)
  message_task('loading embeddings')
  if (is.null(dims)) {
    embedding <- Seurat::Embeddings(object, reduction = reduction)
  } else {
    embedding <-
      Seurat::Embeddings(object, reduction = reduction)[,dims]
  }

  message_task('defining parameters')
  Sys.setenv(OMP_NUM_THREADS = nThreads)
  reticulate::py_capture_output({
    reducer <- UMAP$UMAP(
      a = a,
      angular_rp_forest = angular_rp_forest,
      b = b,
      densmap = densmap,
      dens_frac = dens_frac,
      dens_var_shift = dens_var_shift,
      output_dens = output_dens,
      disconnection_distance = disconnection_distance,
      dens_lambda = dens_lambda,
      force_approximation_algorithm = force_approximation_algorithm,
      init = init,
      learning_rate = learning_rate,
      local_connectivity = local_connectivity,
      low_memory = low_memory,
      metric = metric,
      metric_kwds = metric_kwds,
      min_dist = min_dist,
      n_components = as.integer(n_components),
      n_epochs = as.integer(n_epochs),
      n_neighbors = as.integer(n_neighbors),
      negative_sample_rate = negative_sample_rate,
      output_metric = output_metric,
      output_metric_kwds = output_metric_kwds,
      random_state = as.integer(random_state),
      repulsion_strength = repulsion_strength,
      set_op_mix_ratio = set_op_mix_ratio,
      spread = spread,
      target_metric = target_metric,
      target_metric_kwds = target_metric_kwds,
      target_n_neighbors = target_n_neighbors,
      target_weight = target_weight,
      transform_queue_size = transform_queue_size,
      transform_seed = as.integer(transform_seed),
      unique = unique,
      verbose = verbose
    )
  }, type = 'stdout')



  message_task('performing umap')
  result <- reducer$fit_transform(embedding)
  message_task('done')
  colnames(result) <- paste0('UMAP_', 1:ncol(result))
  rownames(result) <- Seurat::Cells(object)

  if(clustering&reduction_name!='umap'){
    reduction_name <- 'clustering'
  }

  if (return_seurat) {
    object[[reduction_name]] <-
      Seurat::CreateDimReducObject(
        embeddings = result,
        key = 'umap_',
        assay =  Seurat::DefaultAssay(object)
      )
    return(object)
  } else {
    return(result)
  }
}

# plot.umap <- function(x){
#   plot(umap_2~umap_1, x,  pch = 21, bg = '#77777732', col = '#22222233', cex = 0.75)
# }


# # global reference to scipy (will be initialized in .onLoad)
# HDBSCAN <- NULL
# UMAP <- NULL
#
# .onLoad <- function(libname, pkgname) {
#   # use superassignment to update global reference to scipy
#   HDBSCAN <<- reticulate::import(c("hdbscan"), delay_load = TRUE)
#   UMAP <<- reticulate::import(c('umap-learn'), delay_load = TRUE)
# }
#




#' UWOT-UAMP: multithreading capable
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#' @import uwot
#'
#' @examples
uwot <- function(object, ...) {
  UseMethod(generic = 'uwot', object = object)
}



#' UWOT-UAMP: multithreading capable
#'
#' @param x
#' @param n_neighbors
#' @param n_components
#' @param metric
#' @param n_epochs
#' @param learning_rate
#' @param scale
#' @param init
#' @param init_sdev
#' @param spread
#' @param min_dist
#' @param set_op_mix_ratio
#' @param local_connectivity
#' @param bandwidth
#' @param repulsion_strength
#' @param negative_sample_rate
#' @param a
#' @param b
#' @param nn_method
#' @param n_trees
#' @param search_k
#' @param approx_pow
#' @param y
#' @param target_n_neighbors
#' @param target_metric
#' @param target_weight
#' @param pca
#' @param pca_center
#' @param pcg_rand
#' @param fast_sgd
#' @param ret_model
#' @param ret_nn
#' @param ret_extra
#' @param n_threads
#' @param n_sgd_threads
#' @param grain_size
#' @param tmpdir
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
uwot.matrix <- function(x,
                        n_neighbors = 15,
                        n_components = 2,
                        metric = "euclidean",
                        n_epochs = NULL,
                        learning_rate = 1,
                        scale = FALSE,
                        init = "spectral",
                        init_sdev = NULL,
                        spread = 1,
                        min_dist = 0.01,
                        set_op_mix_ratio = 1,
                        local_connectivity = 1,
                        bandwidth = 1,
                        repulsion_strength = 1,
                        negative_sample_rate = 5,
                        a = NULL,
                        b = NULL,
                        nn_method = NULL,
                        n_trees = 50,
                        search_k = 2 * n_neighbors * n_trees,
                        approx_pow = FALSE,
                        y = NULL,
                        target_n_neighbors = n_neighbors,
                        target_metric = "euclidean",
                        target_weight = 0.5,
                        pca = NULL,
                        pca_center = TRUE,
                        pcg_rand = TRUE,
                        fast_sgd = FALSE,
                        ret_model = FALSE,
                        ret_nn = FALSE,
                        ret_extra = c(),
                        n_threads = parallel::detectCores()-1,
                        n_sgd_threads = 1,
                        grain_size = 1,
                        tmpdir = tempdir(),
                        verbose = getOption("verbose", TRUE)) {
  embds <- x
  if (is.null(n_components)) {
    n_components <- ncol(embds)
  }



  umap_res <-   uwot::umap(
    embds,
    n_neighbors = n_neighbors,
    n_components = n_components,
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    scale = scale,
    init = init,
    init_sdev = init_sdev,
    spread = spread,
    min_dist = min_dist,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    bandwidth = bandwidth,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    nn_method = nn_method,
    n_trees = n_trees,
    search_k = search_k,
    approx_pow = approx_pow,
    y = y,
    target_n_neighbors = target_n_neighbors,
    target_metric = target_metric,
    target_weight = target_weight,
    pca = pca,
    pca_center = pca_center,
    pcg_rand = pcg_rand,
    fast_sgd = fast_sgd,
    ret_model = ret_model,
    ret_nn = ret_nn,
    ret_extra = ret_extra,
    n_threads = n_threads,
    n_sgd_threads = n_sgd_threads,
    grain_size = grain_size,
    tmpdir = tmpdir,
    verbose = verbose
  )

  rownames(umap_res) <- rownames(embds)
  colnames(umap_res) <- paste0('UMAP_', 1:ncol(umap_res))

  return(umap_res)


}



#' UWOT-UAMP: multithreading capable
#'
#' @param object
#' @param reduction
#' @param n_neighbors
#' @param n_components
#' @param metric
#' @param n_epochs
#' @param learning_rate
#' @param scale
#' @param init
#' @param init_sdev
#' @param spread
#' @param min_dist
#' @param set_op_mix_ratio
#' @param local_connectivity
#' @param bandwidth
#' @param repulsion_strength
#' @param negative_sample_rate
#' @param a
#' @param b
#' @param nn_method
#' @param n_trees
#' @param search_k
#' @param approx_pow
#' @param y
#' @param target_n_neighbors
#' @param target_metric
#' @param target_weight
#' @param pca
#' @param pca_center
#' @param pcg_rand
#' @param fast_sgd
#' @param ret_model
#' @param ret_nn
#' @param ret_extra
#' @param n_threads
#' @param n_sgd_threads
#' @param grain_size
#' @param tmpdir
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
uwot.Seurat <- function(object,
                        reduction = 'harmony',
                        clustering = FALSE,
                        dims = NULL,
                        n_neighbors = 15,
                        n_components = 2,
                        metric = "euclidean",
                        n_epochs = NULL,
                        learning_rate = 1,
                        scale = FALSE,
                        init = "spectral",
                        init_sdev = NULL,
                        spread = 1,
                        min_dist = 0.01,
                        set_op_mix_ratio = 1,
                        local_connectivity = 1,
                        bandwidth = 1,
                        repulsion_strength = 1,
                        negative_sample_rate = 5,
                        a = NULL,
                        b = NULL,
                        nn_method = NULL,
                        n_trees = 50,
                        search_k = 2 * n_neighbors * n_trees,
                        approx_pow = FALSE,
                        y = NULL,
                        target_n_neighbors = n_neighbors,
                        target_metric = "euclidean",
                        target_weight = 0.5,
                        pca = NULL,
                        pca_center = TRUE,
                        pcg_rand = TRUE,
                        fast_sgd = FALSE,
                        ret_model = FALSE,
                        ret_nn = FALSE,
                        ret_extra = c(),
                        n_threads = parallel::detectCores()-1,
                        n_sgd_threads = 1,
                        grain_size = 1,
                        tmpdir = tempdir(),
                        verbose = getOption("verbose", TRUE),
                        return_seurat = TRUE,
                        reduction_name = 'umap'
                        ) {
  reds <- Seurat::Reductions(object)



  if (!reduction %in% reds) {
    pca_pos <- grepl(pattern = 'PCA',
                     ignore.case = TRUE,
                     x = reds)

    if (any(pca_pos)) {
      if (length(which(pca_pos)) > 1) {
        stop(
          paste(
            reduction,
            'is not found as a reduction and more than 1 reduction PCA detected.'
          )
        )
      }
      reduction <- reds[which(pca_pos)]
    } else {
      stop(paste(
        reduction,
        'is not found as a reduction and no reduction named PCA was found.'
      ))
    }
  }

  embds <- Seurat::Embeddings(object, reduction = reduction)

  # reducing embedding if dims is set
  if(!is.null(dims)){

    # checks if you just want the first 10 dims or a specific set
    if(length(dims==1)){
      embds <- embds[,1:dims]
    } else {
      embds <- embds[,dims]
    }

  }

  # checks for the number of components to return
  if (is.null(n_components)) {
    n_components <- ncol(embds)
  }

  umap_res <-   uwot::umap(
    embds,
    n_neighbors = n_neighbors,
    n_components = n_components,
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    scale = scale,
    init = init,
    init_sdev = init_sdev,
    spread = spread,
    min_dist = min_dist,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    bandwidth = bandwidth,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    nn_method = nn_method,
    n_trees = n_trees,
    search_k = search_k,
    approx_pow = approx_pow,
    y = y,
    target_n_neighbors = target_n_neighbors,
    target_metric = target_metric,
    target_weight = target_weight,
    pca = pca,
    pca_center = pca_center,
    pcg_rand = pcg_rand,
    fast_sgd = fast_sgd,
    ret_model = ret_model,
    ret_nn = ret_nn,
    ret_extra = ret_extra,
    n_threads = n_threads,
    n_sgd_threads = n_sgd_threads,
    grain_size = grain_size,
    tmpdir = tmpdir,
    verbose = verbose
  )

  rownames(umap_res) <- rownames(embds)
  colnames(umap_res) <- paste0('UMAP_', 1:ncol(umap_res))

  if(clustering&reduction_name!='umap'){
    reduction_name <- 'clustering'
  }

  if (return_seurat) {
    object[[reduction_name]] <-
      Seurat::CreateDimReducObject(embeddings = umap_res,
                                   key = 'clustering_',
                                   assay = 'RNA')
    return(object)
  } else {
    return(umap_res)
  }

}
