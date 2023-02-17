#' Easy Add Meta data to Seurat Object
#'
#' @param object
#' @param meta
#' @param col.name   name of the column for the new meta data
#'
#' @return
#' @export
#' @import Seurat
#'
#' @examples
NewMeta <- function(object, meta, col.name){
  test <- meta[as.character(Seurat::Idents(object))]
  names(test) <- colnames(object)
  result <- Seurat::AddMetaData(object, test, col.name)
  return(result)
}

#' Easily add new embeddings to Seurat Object
#'
#' @param object
#' @param embedding a matrix containing the embedding.
#' @param reduction_name name to save the reduction as
#' @param type embedding type (e.g. umap (default), 'harmony','hdbscan','pca')
#'
#' @return
#' @export
#'
#' @examples
add_seurat_embedding <- function(object, embedding,reduction_name, type = 'umap'){
  key_name <- paste0(type,'_')
  colnames(embedding) <- paste0(key_name, 1:ncol(embedding))
  rownames(embedding) <- colnames(object)
  object[[reduction_name]] <- Seurat::CreateDimReducObject(embeddings = embedding, key = 'umap_', assay =  Seurat::DefaultAssay(object))
  return(object)
}

#' Simple method for creating Seurat Objects
#'
#' @param filepath
#' @param sample  sample name to use
#'
#' @return
#' @export
#' @import Seurat
#'
#' @examples
create_seurat <- function(filepath, sample = NULL ){
  if(is.null(sample)){
    sample <- basename(filepath)
  }

  # read in 10X data
  x <- Seurat::Read10X(data.dir = filepath)

  # create unique cell ids
  cell_ids <- paste0(sample, '_', colnames(x))
  colnames(x) <-cell_ids

  # create Seurat Object and include meta data
  suppressWarnings({
    res <- Seurat::CreateSeuratObject(x, meta.data = meta, project = sample)
  })

  return(res)
}

#' Seurat Preprocessing
#'
#' @param object
#' @param species
#' @param nfeatures
#' @param npcs  number of principle component dimensions to calculate
#'
#' @return
#' @export
#' @import Seurat
#' @importFrom stringr str_to_title
#' @import crayon
#'
#' @examples
pre_processing <- function(object, species = 'Homo sapiens', nfeatures = 3000, npcs = 50, filters  = list(nFeature = c(200,3000), percent.mt = 20), verbose = TRUE){
  message_section('Start')
  if(species == 'Homo sapiens'){
    mt_pattern <- '^MT-'
  } else {
    mt_pattern <- '^mt-'
  }

  object <- Seurat::PercentageFeatureSet(object,
                                         pattern = mt_pattern,
                                         col.name = "percent.mt")

  message_task('Filtering out low quality cells and doublets (1/5)')

  # Removing low quality cells and doublets
  object <- subset(object,percent.mt < filters[['percent.mt']] & nFeature_RNA > filters[['nFeature']][1] & nFeature_RNA <  filters[['nFeature']][2] )

  message_task('Normalizing data (2/5)')
  # Normalization
  object<- Seurat::NormalizeData(object, verbose = verbose)
  # Variable Features

  message_task(paste('Finding',nfeatures,'most variable fatures (3/5)'))

  object<- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures, verbose = verbose)
  data('cc.genes.updated.2019', package = 'Seurat')
  if(species =='Mus musculus'){ # change gene name format to title capitalization
    ccss <- stringr::str_to_title(cc.genes.updated.2019$s.genes)
    ccg2m <- stringr::str_to_title(cc.genes.updated.2019$g2m.genes)
  } else { # use build in gene names
    ccss <- cc.genes.updated.2019$s.genes
    ccg2m <- cc.genes.updated.2019$g2m.genes
  }

  # scoring function
  object<- Seurat::CellCycleScoring(object,s.features = ccss, g2m.features = ccg2m, verbose = verbose)

  # difference between s and g2m scores
  object$CC.Difference <-object$S.Score -object$G2M.Score

  message_task('Scaling data (4/5)')
  # Scaling Data ----
  object<- Seurat::ScaleData(object,vars.to.regress = c('CC.Difference','percent.mt'), verbose = verbose)

  message_task('Performing PCA (5/5)')
  message_append(paste('using npcs =',npcs))
  # PCA ----
  object<- Seurat::RunPCA(
    object,
    pc.genes = Seurat::VariableFeatures(object),
    npcs = npcs,
    verbose = verbose)
  return(object)
  message_section('Finished')

}


