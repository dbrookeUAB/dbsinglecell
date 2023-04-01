#' Fast Import of 10X Data
#'
#' @param return Set the format for the imported data. By default, 10X data is imported as a sparse matrix (`return='sparse'`), yet it can also return it as either a `SingleCellExperiment` (`return='sce'`) or  `SeuratObject` (`return='Seurat`).
#' @param gene_col The column position to
#' @param barcode_col
#' @param seurat_project
#' @param QC
#' @param path
#' @param trim_fat
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
                       QC = TRUE,
                       trim_fat = FALSE,
                       verbose = TRUE
){

    THREADS <- parallel::detectCores()
    data.table::setDTthreads(threads =  THREADS)
  message_section('Starting db_read10x', verbose = verbose)
  fl <- dir(path)
  message_task('Importing files', verbose = verbose)
  message_append('matrix.mtx', verbose = verbose)
  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(
    file.path(path, fl[grepl('matrix.mtx$', fl)|grepl('matrix.mtx.gz$', fl)]),
    skip = 3,
    col.names = c('i', 'j', 'value'),
    key = c('i', 'j'),
    colClasses = c('integer', 'integer', 'integer'),
    header = FALSE,
    integer64 = 'numeric',
    nThread = THREADS,
    showProgress = TRUE
  )

  message_append('barcodes.tsv', verbose = verbose)
  # imports barcode ---------------------------------------------------------
  barcode <- data.table::fread(file.path(path, fl[grepl('barcodes.tsv$', fl)|grepl('barcodes.tsv.gz$', fl)]),
                               nThread = THREADS,
                               header = FALSE,
                               colClasses = 'character')[[barcode_col]]
  message_append('features.tsv', verbose = verbose)
  # imports gene ------------------------------------------------------------
  GENE <- data.table::fread(file.path(path, fl[grepl('[gf][e][na][te].+.tsv$', fl)|grepl('[gf][e][na][te].+.tsv.gz$', fl)]),
                            nThread = THREADS,
                            header = FALSE)

  gene <- GENE[[gene_col]]

  message_append('done', verbose = verbose)
  # duplicate gene names for row names --------------------------------------

  if (any(duplicated(gene) == TRUE)) {
    dg <- data.table::data.table(position = which(duplicated(gene)),
                     name = gene[duplicated(gene)])[, N := .N, name][]
    dg[, new.name := paste0(name, '.', 1:.N), name]
    gene[dg$position] <- dg$new.name
  }

  message_task('Assembling sparse matrix', verbose = verbose)
  # creates sparse matrix ---------------------------------------------------

  COUNTS <-  Matrix::sparseMatrix(
    i =  mat$i,
    j = mat$j,
    x = mat$value,
    dimnames = list(gene[1:(max(mat$i))], barcode[1:(max(mat$j))])
  )
  if(trim_fat==TRUE){
    COUNTS <- COUNTS[Matrix::rowSums(COUNTS)!=0,]
  }


  # object creation ---------------------------------------------------------

  #  Seurat
  if (return == 'Seurat') {
    message_task('Creating Seurat Object', verbose = verbose)
    res <- Seurat::RenameCells(Seurat::CreateSeuratObject(counts = COUNTS, project = seurat_project), add.cell.id = basename(path))

    #  SingleCellExperiment
  } else if (return == 'sce') {
    message_task('Creating SingleCellExperiemtn Object', verbose = verbose)
    res <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = COUNTS))
  } else if (return == 'sparse') {
    res <- COUNTS
  } else if(return =='dbsc'){
    res <- data.table::data.table(gene_name = gene[mat$i],
                                  barcodes = barcode[mat$j],
                                  count = mat$value, key = c('gene_name','barcodes'))
  }

  message_task('done', verbose = verbose)
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
                           showProgress = TRUE,
                           col.names = c('i','j','value'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE,
                           nThread = parallel::detectCores()
                           )
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
                           showProgress = TRUE,
                           skip = 3,
                           col.names = c('i','j','value'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE,
                           nThread = parallel::detectCores() )

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

#' Fast Import of 10X Data (experimental)
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
#' db_read10x_experimental  <- function(path){
#'   fl <- dir(path)
#'
#'   # reads in matrix file ----------------------------------------------------
#'   mat <- data.table::fread(file.path(path,fl[grepl('^matrix.mtx',fl)]),
#'                            skip = 3,
#'                            col.names = c('i','j','value'),
#'                            key = c('i','j'),
#'                            colClasses = c('integer','integer','integer'),
#'                            header = FALSE)
#'
#'   # imports barcode ---------------------------------------------------------
#'   barcode <- data.table::fread(file.path(path, fl[grepl('^barcodes.tsv',fl)]),
#'                                header = FALSE,
#'                                colClasses = 'character')$V1
#'
#'   # imports gene ------------------------------------------------------------
#'   gene<- data.table::fread(
#'     file.path(path,fl[grepl('^[gf][e][na][te].+.tsv',fl)]),
#'     header = FALSE)
#'   colnames(gene) <- c('gene_id','gene_name','gene_type')[1:ncol(gene)]
#'
#'
#'   # duplicate gene names for row names --------------------------------------
#'
#'   # if(!all(duplicated(gene[['gene_id']])==FALSE)){
#'   #   dg <- data.table(
#'   #     position = which(duplicated(gene)),
#'   #     name = gene[['gene_id']][duplicated(gene[['gene_id']])])[,N:=.N,name][]
#'   #   dg[,new.name:=paste0(name,'.',1:.N), name]
#'   #   gene[['gene_id']][dg$position] <- dg$new.name
#'   # }
#'
#'   res <- data.table::data.table(gene_id = gene[['gene_id']][mat$i],
#'                                 gene_name = gene[['gene_name']][mat$i],
#'                                 barcodes = barcode[mat$j],
#'                                 count = mat$value, key = c('gene_id','barcodes'))
#'   res$gene_id <- as.factor(res$gene_id)
#'   res$gene_name <- as.factor(res$gene_name)
#'   res$barcodes <- as.factor(res$barcodes)
#'   data.table::setkey(res, gene_name, barcodes)
#'
#'   return(res)
#'
#' }
#'
#'
#'
#'
#'
#' #' Merge dbsc objects
#' #'
#' #' @param ...
#' #' @param return.DT
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' merge_dbsc <- function(..., return.DT = FALSE){
#'   obj_list <- list(...)
#'   if(length(obj_list)==1 && is.list(obj_list[[1]])){
#'     result <- data.table::rbindlist(obj_list[[1]])
#'   } else {
#'     result <-   data.table::rbindlist(obj_list)
#'   }
#'   data.table::setkey(result, gene_name, barcodes)
#'   if(return.DT==FALSE){
#'     res <- Matrix::sparseMatrix(i = as.integer(result$gene_name), j = as.integer(result$barcodes), x = result$count, repr = 'C')
#'     rownames(res) <- levels(result$gene_name)
#'     colnames(res) <- levels(result$barcodes)
#'     return(res)
#'   } else {
#'     return(result)
#'   }
#' }




