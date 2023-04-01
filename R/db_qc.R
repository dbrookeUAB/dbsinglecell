#' Add QC metrics to Seurat Object
#'
#' @param object
#' @param species
#' @param subset
#' @param mito
#' @param ribo
#' @param Complexity
#' @param min.features
#'
#' @return
#' @export
#'
#' @examples
db_qc <- function(object,
                  species = 'Mus musculus',
                  subset = FALSE, mito = 10,ribo = 1,Complexity = 0.8, min.features = 500, cell_cycle = FALSE){
  if(species=='Mus musculus'){
    mt_genes <- '^mt-'
    ribo_genes <- '^Rp[sl]'

  } else if(species=='Homo sapiens'){
    mt_genes <- '^MT-'
    ribo_genes <- '^RP[SL]'

  } else {
    stop('incorrect species')
  }

  object <- Seurat::PercentageFeatureSet(Seurat::PercentageFeatureSet(object, pattern = mt_genes,col.name = 'percent.mt'), pattern = ribo_genes, col.name = 'percent.ribo')
  object[['complexity']] <- log10(object$nFeature_RNA)/log10(object$nCount_RNA)

  if(cell_cycle==TRUE){
    if(!exists('cc.genes')){
      data('cc.genes',package = 'dbsinglecell')
    }

    object <- Seurat::CellCycleScoring(object = object,
                                       s.features = cc.genes[[species]]$s.genes,
                                       g2m.features = cc.genes[[species]]$g2m.genes)
    object$CC.Difference <-object$S.Score -object$G2M.Score
  }

  if(subset==TRUE){
    object <- subset(object, percent.mt < mito & percent.ribo > ribo & complexity > Complexity & nFeature_RNA > min.features & percent.mt > 0)
  }
  return(object)
}
