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
