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
