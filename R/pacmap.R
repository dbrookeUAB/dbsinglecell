#' Pairwise Controlled Manifold Approximation (R wrapper)
#'
#' @param x
#' @param n_components
#' @param n_neighbors
#' @param MN_ratio
#' @param FP_ratio
#' @param init
#' @param num_iters
#' @param pair_neighbors
#' @param pair_MN
#' @param pair_FP
#' @param lr
#' @param apply_pca
#' @param intermediate
#' @param nThreads
#'
#' @return
#' @export
#'
#' @examples
PaCMAP <- function(object, ...){
  UseMethod(generic = 'PaCMAP', object = object)
}


#' Pairwise Controlled Manifold Approximation (R wrapper)
#'
#' @param x
#' @param n_components
#' @param n_neighbors
#' @param MN_ratio
#' @param FP_ratio
#' @param num_iters
#' @param pair_neighbors
#' @param pair_MN
#' @param pair_FP
#' @param lr
#' @param apply_pca
#' @param intermediate
#' @param nThreads
#'
#' @return
#' @export
#'
#' @examples
PaCMAP.default <- function(x,
                           n_components = 2,
                           n_neighbors = 10,
                           MN_ratio = 0.5,
                           FP_ratio = 2,
                           num_iters = 450,
                           pair_neighbors = NULL,
                           pair_MN = NULL,
                           pair_FP = NULL,
                           lr = 1,
                           apply_pca = TRUE,
                           intermediate = FALSE,
                           nThreads = parallel::detectCores()-1
){

  pacmap <- reticulate::import('pacmap', delay_load = FALSE)

  if(nrow(x) > 10000){
    n_neighbors <- round(10 + 15 * (log10(nrow(x)) - 4))
  }

  EMBEDDING <- pacmap$PaCMAP(n_components = as.integer(n_components),
                             n_neighbors = as.integer(n_neighbors),
                             MN_ratio = MN_ratio,
                             FP_ratio = FP_ratio,
                             num_iters = as.integer(num_iters),
                             pair_neighbors = pair_neighbors,
                             pair_MN = pair_MN,
                             pair_FP = pair_FP,
                             lr = as.integer(lr),
                             apply_pca = apply_pca,
                             intermediate = intermediate
  )

  result <- EMBEDDING$fit_transform(X = x)
  colnames(result) <- paste0('PaCMAP_',1:ncol(result))
  rownames(result) <- rownames(x)
  # reticulate::py_main_thread_func(clusterer$fit(x))

  return(result)
}

#' Pairwise Controlled Manifold Approximation (R wrapper)
#'
#' @param x
#' @param n_components
#' @param n_neighbors
#' @param MN_ratio
#' @param FP_ratio
#' @param num_iters
#' @param pair_neighbors
#' @param pair_MN
#' @param pair_FP
#' @param lr
#' @param apply_pca
#' @param intermediate
#' @param nThreads
#'
#' @return
#' @export
#'
#' @examples
PaCMAP.matrix <- function(x,
                          n_components = 2,
                          n_neighbors = 10,
                          MN_ratio = 0.5,
                          FP_ratio = 2,
                          num_iters = 450,
                          pair_neighbors = NULL,
                          pair_MN = NULL,
                          pair_FP = NULL,
                          lr = 1,
                          apply_pca = TRUE,
                          intermediate = FALSE,
                          nThreads = parallel::detectCores()-1
){

  pacmap <- reticulate::import('pacmap', delay_load = FALSE)

  if(nrow(x) > 10000){
    n_neighbors <- round(10 + 15 * (log10(nrow(x)) - 4))
  }

  EMBEDDING <- pacmap$PaCMAP(n_components = as.integer(n_components),
                             n_neighbors = as.integer(n_neighbors),
                             MN_ratio = MN_ratio,
                             FP_ratio = FP_ratio,
                             num_iters = as.integer(num_iters),
                             pair_neighbors = pair_neighbors,
                             pair_MN = pair_MN,
                             pair_FP = pair_FP,
                             lr = as.integer(lr),
                             apply_pca = apply_pca,
                             intermediate = intermediate
  )

  result <- EMBEDDING$fit_transform(X = x)
  colnames(result) <- paste0('PaCMAP_',1:ncol(result))
  rownames(result) <- rownames(x)
  # reticulate::py_main_thread_func(clusterer$fit(x))

  return(result)

}



#' Pairwise Controlled Manifold Approximation (R wrapper)
#'
#' @param n_components
#' @param n_neighbors
#' @param MN_ratio
#' @param FP_ratio
#' @param num_iters
#' @param pair_neighbors
#' @param pair_MN
#' @param pair_FP
#' @param lr
#' @param apply_pca
#' @param intermediate
#' @param nThreads
#' @param object
#' @param reduction
#' @param reduction_name
#' @param dims
#' @param return_seurat
#'
#' @return
#' @export
#'
#' @examples
PaCMAP.Seurat <- function(object,
                          reduction = 'pca',
                          reduction_name = 'pacmap',
                           dims = NULL,
                          n_components = 2,
                          n_neighbors = 10,
                          MN_ratio = 0.5,
                          FP_ratio = 2,
                          num_iters = 450,
                          pair_neighbors = NULL,
                          pair_MN = NULL,
                          pair_FP = NULL,
                          lr = 1,
                          apply_pca = TRUE,
                          intermediate = FALSE,
                           return_seurat = TRUE,
                          nThreads = parallel::detectCores() - 1
){

  if(is.null(dims)){
    x <- Seurat::Embeddings(object, reduction = reduction)
  } else {
    x <- Seurat::Embeddings(object, reduction = reduction)[,dims]
  }
  Sys.setenv(OMP_NUM_THREADS = nThreads)
  pacmap <- reticulate::import('pacmap', delay_load = FALSE)

  if(nrow(x) > 10000){
    n_neighbors <- round(10 + 15 * (log10(nrow(x)) - 4))
  }

  EMBEDDING <- pacmap$PaCMAP(n_components = as.integer(n_components),
                             n_neighbors = as.integer(n_neighbors),
                             MN_ratio = MN_ratio,
                             FP_ratio = FP_ratio,
                             num_iters = as.integer(num_iters),
                             pair_neighbors = pair_neighbors,
                             pair_MN = pair_MN,
                             pair_FP = pair_FP,
                             lr = as.integer(lr),
                             apply_pca = apply_pca,
                             intermediate = intermediate
  )

  result <- EMBEDDING$fit_transform(X = x)
  colnames(result) <- paste0('PaCMAP_',1:ncol(result))
  rownames(result) <- Seurat::Cells(object)
  # reticulate::py_main_thread_func(clusterer$fit(x))

  if (return_seurat==TRUE){
    object[[reduction_name]] <-
      Seurat::CreateDimReducObject(
        embeddings = result,
        key = 'pacmap_',
        assay =  Seurat::DefaultAssay(object)
      )
    return(object)
  } else {
    return(result)
  }
  }



