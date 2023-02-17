#' Remove empty columns from sparse matrix
#'
#' @param M
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
db_mtx_remove_empty <- function(M, ...){
  M[,Matrix::colSums(M != 0) != 0]
}



#' Convert a sparse matrix to a data.table
#'
#' @param matrix  sparse matrix to be used
#' @param with.names  return the row and column names from the sparse matrix [default:FALSE]
#' @param single.cell returns the row names as Genes and column names as cells if from a single cell dataset [default: FALSE]
#' @param key logical for whether the data.table key should be set [default:FALSE]
#'
#' @return
#' @export
#'
#' @import data.table
#' @import Seurat
#' @import Matrix
#'
#' @examples
#'
#'
sparse2DT <- function(matrix, with.names = FALSE,single.cell = FALSE, key = TRUE){

    i <- as.integer(matrix@i+1)
    j  <- as.integer(findInterval(seq_along(matrix@x)-1,matrix@p[-1])+1)
    x  <- matrix@x

    if(single.cell==TRUE){
      with.names <- TRUE
    }


    # creating workable dataset of count data
    if(with.names==FALSE){
      result <- data.table::data.table(i,j,x)
    } else {
      rn <- matrix@Dimnames[[1]][i]
      cn <- matrix@Dimnames[[2]][j]

      if(single.cell==TRUE){
        result <- data.table::data.table(Genes = rn,
                                         Cell = cn,
                                         i,
                                         j,
                                         Value = x)
      } else {
        result <- data.table::data.table(rn,
                                         cn,
                                         i,
                                         j,
                                         Value = x)
      }
    }

if(key==TRUE){
  if(with.names==FALSE){
    data.table::setkey(result, i,j)
  } else if(single.cell==TRUE){
    data.table::setkey(result, Genes,Cell)
  } else {
    setkeyv(result, rn, cn)
  }

}

  return(result)

}


#' Convert a sparse matrix to a data.table
#'
#' @param object Seurat object
#'
#' @return
#' @export
#' @import data.table
#' @import Seurat
#' @importMatrix
#'
#' @examples
sparse2DT.Seurat <- function(object){

  # exporting normalized data
  mat <- object@assays$RNA@data[, Cells(object)]

  result <- sparse2DT(mat)
  return(result)
}
