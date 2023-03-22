#' Quality Control
#'
#' @param object
#'
#' @return
#' @export
#' @importFrom
#'
#' @examples
db_count_genes <- function(object){

# counting the total number of genes with counts per barcode ----------------------------------
  result <- object[,.(genes = .N),barcodes]

# setting keys to allowing joining by barcodes ------------------------------------------------
  data.table::setkey(result, barcodes)
  return(result)
}

#' Quality Control
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
db_count_reads <- function(object){

# counting the total number of reads per barcode ----------------------------------------------
  result <- object[,.(reads = sum(count)),barcodes]

# setting keys to allowing joining by barcodes ------------------------------------------------

  data.table::setkey(result, barcodes)
  return(result)
}


#' Quality Control
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
db_qc_cell <- function(object){
  result <- object[,.(genes = .N, reads = sum(count)),barcodes]
  result[,genes.out:=scater::isOutlier(genes)][]
  result[,reads.out:=scater::isOutlier(reads)][]
  setkey(result,barcodes)
  return(result)
}


Outliers <- function(object){
  t1 <- scuttle::isOutlier(object$nFeature_RNA/log(object$nCount_RNA), log = T)
  t2 <- scuttle::isOutlier(log(object$nCount_RNA)/object$nFeature_RNA,log = T)
  t3 <- scuttle::isOutlier(object$nFeature_RNA)
  t4 <- scuttle::isOutlier(log(object$nCount_RNA),log = F)

  res <- t1|t2|t3|t4
  return(res)
}
