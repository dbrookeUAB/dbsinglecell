#' Fast Import of 10X Data
#'
#' @param return Set the format for the imported data. By default, 10X data is imported as a sparse matrix (`return='sparse'`), yet it can also return it as either a `SingleCellExperiment` (`return='sce'`) or  `SeuratObject` (`return='Seurat`).
#' @param gene_col The column position to
#' @param barcode_col
#' @param seurat_project
#' @param QC
#' @param path
#'
#' @return
#' @export
#' @import data.table
#' @import Matrix
#' @import SingleCellExperiment
#'
#' @examples
db_read10x <- function(path,
                       return = 'Seurat',
                       gene_col = 2,
                       barcode_col = 1,
                       seurat_project = basename(path),
                       QC = TRUE
){

  fl <- dir(path)

  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(
    file.path(path, fl[grepl('^matrix.mtx', fl)]),
    skip = 3,
    col.names = c('i', 'j', 'value'),
    key = c('i', 'j'),
    colClasses = c('integer', 'integer', 'integer'),
    header = FALSE,
    integer64 = 'numeric'
  )

  # imports barcode ---------------------------------------------------------
  barcode <- data.table::fread(file.path(path, fl[grepl('^barcodes.tsv', fl)]),
                               header = FALSE,
                               colClasses = 'character')[[barcode_col]]

  # imports gene ------------------------------------------------------------
  GENE <- data.table::fread(file.path(path, fl[grepl('^[gf][e][na][te].+.tsv', fl)]),
                            header = FALSE)

  gene <- GENE[[gene_col]]

  # duplicate gene names for row names --------------------------------------

  if (any(duplicated(gene) == TRUE)) {
    dg <- data.table(position = which(duplicated(gene)),
                     name = gene[duplicated(gene)])[, N := .N, name][]
    dg[, new.name := paste0(name, '.', 1:.N), name]
    gene[dg$position] <- dg$new.name
  }


  # creates sparse matrix ---------------------------------------------------

  COUNTS <-  Matrix::sparseMatrix(
    i =  mat$i,
    j = mat$j,
    x = mat$value,
    dimnames = list(gene[1:(max(mat$i))], barcode[1:(max(mat$j))])
  )

  # object creation ---------------------------------------------------------

  #  Seurat
  if (return == 'Seurat') {
    res <- Seurat::RenameCells(Seurat::CreateSeuratObject(counts = COUNTS, project = seurat_project), add.cell.id = basename(path))

    #  SingleCellExperiment
  } else if (return == 'sce') {
    res <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = COUNTS))
  } else if (return == 'sparse') {
    res <- COUNTS
  }

  return(res)
}


#' Read Count data from Salmon Alevin
#'
#' @param path
#' @param return
#'
#' @return
#' @export
#'
#' @examples
db_read_alevin  <- function(path = '.', return = 'Seurat'){
  DIR <- path
  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(file.path(DIR,'quants_mat.mtx'),
                           skip = 3,
                           col.names = c('i','j','value'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE)
  # data.table::setcolorder(mat, 'i')
  # imports barcode ---------------------------------------------------------
  barcode <- readLines(file.path(DIR, 'quants_mat_rows.txt'))

  # imports gene ------------------------------------------------------------
  gene<- readLines(file.path(DIR, 'quants_mat_cols.txt'))

  # duplicate gene names for row names --------------------------------------

  if(!all(duplicated(gene)==FALSE)){
    dg <- data.table(
      position = which(duplicated(gene)),
      name = gene[duplicated(gene)])[,N:=.N,name][]
    dg[,new.name:=paste0(name,'.',1:.N), name]
    gene[dg$position] <- dg$new.name
  }
  max_i <- max(mat$j)
  res <-  (Matrix::sparseMatrix(
    i =  mat$j,
    j = mat$i,
    x = mat$value,
    dimnames = list(gene[1:max_i],barcode),
    repr = 'C'
  ))

  if(return=='SCE'){
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = res))
  } else if(return=='Seurat') {
    Seurat::CreateSeuratObject(counts = res)
  } else {
    return(res)
  }

}


#' Import `.mtx` files
#'
#' @param mtx
#' @param features
#' @param cells
#' @param return
#' @param transpose_mtx
#'
#' @return
#' @export
#'
#' @examples
db_read_mtx  <- function(mtx, features, cells, return = 'matrix', transpose_mtx = FALSE){

  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(mtx,
                           skip = 3,
                           col.names = c('i','j','value'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE)
  # data.table::setcolorder(mat, 'i')
  # imports barcode ---------------------------------------------------------
  barcode <- readLines(features)

  # imports gene ------------------------------------------------------------
  gene<- readLines(cells)

  # duplicate gene names for row names --------------------------------------

  if(!all(duplicated(gene)==FALSE)){
    dg <- data.table(
      position = which(duplicated(gene)),
      name = gene[duplicated(gene)])[,N:=.N,name][]
    dg[,new.name:=paste0(name,'.',1:.N), name]
    gene[dg$position] <- dg$new.name
  }
  max_i <- max(mat$i)
  res <-  Matrix::sparseMatrix(
    i =  mat$i,
    j = mat$j,
    x = mat$value,
    dimnames = list(gene[1:max_i],barcode),
    repr = 'C')

  if(transpose_mtx==TRUE){
    res <- t(res)
  }

  if(return=='SCE'){
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = res))
  } else if(return=='Seurat') {
    Seurat::CreateSeuratObject(counts = res)
  } else {
    return(res)
  }

}


