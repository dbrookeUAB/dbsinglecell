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
pre_processing <- function(object, species = 'Homo sapiens', nfeatures = 3000, npcs = 50){
  if(species == 'Homo sapiens'){
    mt_pattern <- '^MT-'
  } else {
    mt_pattern <- '^mt-'
  }

  object <- Seurat::PercentageFeatureSet(object,
                                 pattern = mt_pattern,
                                 col.name = "percent.mt")

  message_section('Filtering out low quality cells and doublets')

  # Removing low quality cells and doublets
  object <- subset(object,percent.mt < 20 &nFeature_RNA >500 & nFeature_RNA < 4100)

  message_section('Normalizing data')
  # Normalization
  object<- Seurat::NormalizeData(object, verbose = TRUE)
  # Variable Features

  message_section(paste('Finding',nfeatures,'most variable fatures'))
  object<- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures)

  if(species =='Mus musculus'){ # change gene name format to title capitalization
    ccss <- stringr::str_to_title(cc.genes.updated.2019$s.genes)
    ccg2m <- stringr::str_to_title(cc.genes.updated.2019$g2m.genes)
  } else { # use build in gene names
    ccss <- cc.genes.updated.2019$s.genes
    ccg2m <- cc.genes.updated.2019$g2m.genes
  }

  # scoring function
  object<- Seurat::CellCycleScoring(object,s.features = ccss, g2m.features = ccg2m)

  # difference between s and g2m scores
  object$CC.Difference <-object$S.Score -object$G2M.Score

  message_section('Scaling data')
  # Scaling Data ----
  object<- Seurat::ScaleData(object,vars.to.regress = c('CC.Difference','percent.mt'))

  message_section('Performing PCA')
  message_append(paste('using npcs =',npcs))
  # PCA ----
  object<- Seurat::RunPCA(
    object,
    pc.genes =object@var.genes,
    npcs = npcs)
  return(object)

}

message_section <- function(text){
  n <- ceiling(options()$width*0.75)
  if(n >120){
    n <- 120
  }
  cat("\n",rep('-',n), "\n",sep = '')
  cat(crayon::bold(crayon::yellow(paste0('[',Sys.time(),']'))), crayon::bold(crayon::green(text)),'\n')
}

message_task <- function(text){
  n <- ceiling(options()$width*0.75)
  if(n >120){
    n <- 120
  }

  if(nchar(text) > n -22){
    cat(crayon::yellow(paste0('[',Sys.time(),']')),'\n')
  } else {
    cat(crayon::yellow(paste0('[',Sys.time(),']')), text,'\n')
  }

}

message_append <- function(text){
  n <- ceiling(options()$width*0.75)
  if(n >120){
    n <- 120
  }

  if(nchar(text) > n - 22){
    invisible()
  } else {
    cat(rep(' ',23),crayon::silver('- '),crayon::silver(text),'\n', sep = '')
  }
}
