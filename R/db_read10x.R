#' Read10x Atlas
#'
#' @param filepaths
#' @param project
#' @param meta  meta data to include with the various datasets
#'
#' @return
#' @export
#' @import doParallel
#' @import foreach
#' @import doSNOW
#' @import snow
#' @import progress
#'
#' @examples
read10x_atlas <- function(filepaths, project = 'scRNAseq', meta = NULL){
  int_list <- 1:length(filepaths)

# checking meta data ------------------------------------------------------
  #  if(is.null(meta)){
  #   meta = list()
  # } else if(nrow(meta)!=length(filepaths)){
  #   stop('meta data needs to be the same length as filepaths')
  # } else {
  #   meta <- as.list(meta)
  # }

# setting project vector --------------------------------------------------
  # if(length(project)!=length(filepaths)){
  #   if( length(project) == 1){
  #     project <- rep(project, times = length(filepaths))
  #   } else {
  #     stop('supply either one project or a vector the same length as filepaths')
  #   }
  # }

# creating cluster and registering doSNOW ---------------------------------
  numCores <- parallel::detectCores() -1
  cl <- snow::makeCluster(numCores)
  doSNOW::registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl))
  e <- simpleError("error occured")

# progress bar ------------------------------------------------------------
  iterations <- length(int_list)                               # used for the foreach loop

  pb <- progress::progress_bar$new(
    format = ":percent item = :item [:bar] :elapsed | eta: :eta",
    total = iterations,
    width = floor(options()$width*0.9),
    clear = TRUE
  )

  # allowing progress bar to be used in foreach -----------------------------

  progress <- function(n) {
    pb$tick(tokens = list(item = int_list[n]))     # report the int_list item
  }

  opts <- list(progress = progress)  # used in the the foreach loop

    result <- foreach::foreach( i = 1:iterations,
                       .options.snow = opts,
                       .export = 'db_read10x',
                       .combine = 'cbind',
                       .packages = c('data.table','SingleCellExperiment','Matrix')) %dopar% {
                         db_read10x(path = filepaths[i])
                       }


  return(result)
}


#' Read10x (v1)
#'
#' @param path
#' @param return.sce  return result as SingleCellExperiment object
#'
#' @return
#' @export
#' @import data.table
#' @import Matrix
#' @import SingleCellExperiment
#'
#' @examples
db_read10x  <- function(path, return.sce = TRUE){
  fl <- dir(path)

  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(file.path(path,fl[grepl('^matrix.mtx',fl)]),
                           skip = 3,
                           col.names = c('i','j','value'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE)

  # imports barcode ---------------------------------------------------------
  barcode <- data.table::fread(file.path(path, fl[grepl('^barcodes.tsv',fl)]),
                               header = FALSE,
                               colClasses = 'character')$V1

  # imports gene ------------------------------------------------------------
  gene<- data.table::fread(
    file.path(path,fl[grepl('^[gf][e][na][te].+.tsv',fl)]),
    header = FALSE)$V1


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
    dimnames = list(gene[1:max_i],barcode))

  if(return.sce){
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = res))
  } else {
    return(gene)
  }

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



#' Read 10X Genomics Datq (V2)
#'
#' @param path
#' @param return
#'
#' @return
#' @export
#'
#' @examples
db_read10x_v2 <- function(path, return = 'sce'){
  fl <- dir(path)

  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(file.path(path,fl[grepl('^matrix.mtx',fl)]),
                           skip = 3,
                           col.names = c('i','j','value'),
                           key = c('i','j'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE)

  # imports barcode ---------------------------------------------------------
  barcode <- data.table::fread(file.path(path, fl[grepl('^barcodes.tsv',fl)]),
                               header = FALSE,
                               colClasses = 'character')$V1

  # imports gene ------------------------------------------------------------
  gene<- data.table::fread(
    file.path(path,fl[grepl('^[gf][e][na][te].+.tsv',fl)]),
    header = FALSE)$V1


  # duplicate gene names for row names --------------------------------------

  if(!all(duplicated(gene)==FALSE)){
    dg <- data.table(
      position = which(duplicated(gene)),
      name = gene[duplicated(gene)])[,N:=.N,name][]
    dg[,new.name:=paste0(name,'.',1:.N), name]
    gene[dg$position] <- dg$new.name
  }

  if(return %in% c('sparse','sce')){
    max_i <- max(mat$i)
    res <-  Matrix::sparseMatrix(
      i =  mat$i,
      j = mat$j,
      x = mat$value,
      dimnames = list(gene[1:max_i],barcode))

    if(return=='sce'){
      res <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = res))
    }
  } else if(return=='DT'){
    res <- data.table::data.table(genes = gene[mat$i],
                                  barcodes = barcode[mat$j],
                                  count = mat$value, key = c('genes','barcodes'))
  } else if(return=='raw'){
    res <- list(data = mat, genes = gene, barcodes = barcode)
  }
  return(res)

}



#' Read10X (v3)
#'
#' @param path
#' @param return
#'
#' @return
#' @export
#'
#' @examples
db_read10x_v3  <- function(path, return = 'sce'){
  fl <- dir(path)

  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(file.path(path,fl[grepl('^matrix.mtx',fl)]),
                           skip = 3,
                           col.names = c('i','j','value'),
                           key = c('i','j'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE)

  # imports barcode ---------------------------------------------------------
  barcode_info <- data.table::fread(file.path(path, fl[grepl('^barcodes.tsv',fl)]),
                               header = FALSE,
                               colClasses = 'character',
                               col.names = c('barcodes')
                               )

  # imports gene ------------------------------------------------------------
  gene_info <- data.table::fread(
    file.path(path,fl[grepl('^[gf][e][na][te].+.tsv',fl)]),
    header = FALSE)

  colnames(gene_info) <- c('gene_id','gene_name','gene_type','gene_type')[1:ncol(gene_info)]


  # duplicate gene names for row names --------------------------------------

  if(!all(duplicated(gene_info[['gene_id']])==FALSE)){
    dg <- data.table::data.table(
      position = which(duplicated(gene_info[['gene_id']])),
      name = gene_info[['gene_id']][duplicated(genegene_info[['gene_id']])])[,N:=.N,name][]
    dg[,new.name:=paste0(name,'.',1:.N), name]
    gene_info[['gene_id']][dg$position] <- dg$new.name
  }

  if(return %in% c('sparse','sce')){
    max_i <- max(mat$i)
    res <-  Matrix::sparseMatrix(
      i =  mat$i,
      j = mat$j,
      x = mat$value,
      dimnames = list(gene_info[['gene_id']][1:max_i],barcode_info[['barcodes']]))

    if(return=='sce'){
      res <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = res),
        rowData = gene_info,
        colData = barcode_info
        )
    }
  } else if(return=='DT'){
    res <- data.table::data.table(gene_id = gene_info[['gene_id']][mat$i],
                                  barcodes = barcode_info[['barcodes']][mat$j],
                                  count = mat$value, key = c('gene_id','barcodes'))
  } else if(return=='raw'){
    res <- dbsc(data = mat, genes = gene_info, barcodes = barcode_info)
  }
  return(res)

}



#' Read10x (v4)
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
db_read10x_v4  <- function(path){
  fl <- dir(path)

  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(file.path(path,fl[grepl('^matrix.mtx',fl)]),
                           skip = 3,
                           col.names = c('i','j','value'),
                           key = c('i','j'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE)

  # imports barcode ---------------------------------------------------------
  barcode <- data.table::fread(file.path(path, fl[grepl('^barcodes.tsv',fl)]),
                               header = FALSE,
                               colClasses = 'character')$V1

  # imports gene ------------------------------------------------------------
  gene<- data.table::fread(
    file.path(path,fl[grepl('^[gf][e][na][te].+.tsv',fl)]),
    header = FALSE)$V1


  # duplicate gene names for row names --------------------------------------

  if(!all(duplicated(gene)==FALSE)){
    dg <- data.table(
      position = which(duplicated(gene)),
      name = gene[duplicated(gene)])[,N:=.N,name][]
    dg[,new.name:=paste0(name,'.',1:.N), name]
    gene[dg$position] <- dg$new.name
  }

  res <- data.table::data.table(genes = gene[mat$i],
                                barcodes = barcode[mat$j],
                                count = mat$value, key = c('genes','barcodes'))
  return(res)

}



#' Read10x (v6)
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
db_read10x_v6  <- function(path){
  fl <- dir(path)

  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(file.path(path,fl[grepl('^matrix.mtx',fl)]),
                           skip = 3,
                           col.names = c('i','j','value'),
                           key = c('i','j'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE)

  # imports barcode ---------------------------------------------------------
  barcode <- data.table::fread(file.path(path, fl[grepl('^barcodes.tsv',fl)]),
                               header = FALSE,
                               colClasses = 'character')$V1

  # imports gene ------------------------------------------------------------
  gene<- data.table::fread(
    file.path(path,fl[grepl('^[gf][e][na][te].+.tsv',fl)]),
    header = FALSE)
  colnames(gene) <- c('gene_id','gene_name','gene_type')[1:ncol(gene)]


  # duplicate gene names for row names --------------------------------------

  # if(!all(duplicated(gene[['gene_id']])==FALSE)){
  #   dg <- data.table(
  #     position = which(duplicated(gene)),
  #     name = gene[['gene_id']][duplicated(gene[['gene_id']])])[,N:=.N,name][]
  #   dg[,new.name:=paste0(name,'.',1:.N), name]
  #   gene[['gene_id']][dg$position] <- dg$new.name
  # }

  res <- data.table::data.table(gene_id = gene[['gene_id']][mat$i],
                                gene_name = gene[['gene_name']][mat$i],
                                barcodes = barcode[mat$j],
                                count = mat$value, key = c('gene_id','barcodes'))
  return(res)

}




#' Read10X (v5)
#'
#' @param path
#' @param return
#'
#' @return
#' @export
#'
#' @examples
db_read10x_v5  <- function(path, return = 'sce'){
  fl <- dir(path)

  # reads in matrix file ----------------------------------------------------
  mat <- data.table::fread(file.path(path,fl[grepl('^matrix.mtx',fl)]),
                           skip = 3,
                           col.names = c('i','j','value'),
                           key = c('i','j'),
                           colClasses = c('integer','integer','integer'),
                           header = FALSE)

  # imports barcode ---------------------------------------------------------
  barcode_info <- data.table::fread(file.path(path, fl[grepl('^barcodes.tsv',fl)]),
                                    header = FALSE,
                                    colClasses = 'character',
                                    col.names = c('barcodes')
  )

  # imports gene ------------------------------------------------------------
  gene_info <- data.table::fread(
    file.path(path,fl[grepl('^[gf][e][na][te].+.tsv',fl)]),
    header = FALSE)

  colnames(gene_info) <- c('gene_id','gene_name','gene_type','gene_type')[1:ncol(gene_info)]


  # duplicate gene names for row names --------------------------------------

  if(!all(duplicated(gene_info[['gene_id']])==FALSE)){
    dg <- data.table::data.table(
      position = which(duplicated(gene_info[['gene_id']])),
      name = gene_info[['gene_id']][duplicated(genegene_info[['gene_id']])])[,N:=.N,name][]
    dg[,new.name:=paste0(name,'.',1:.N), name]
    gene_info[['gene_id']][dg$position] <- dg$new.name
  }

  res <- dbsc2(data = mat, genes = gene_info, barcodes = barcode_info)


  return(res)

}


#' `dbsc` data class
#'
#' @slot data data.table.
#' @slot genes character.
#' @slot barcodes character.
#'
#' @return
#' @export
#'
#' @examples
dbsc <- setClass(Class = "dbsc",
                 slots = c(data = 'data.table',
                           genes = 'data.table',
                           barcodes = 'data.table'
                 ),
                 prototype = list(data = data.table::data.table(i = integer(),
                                                                j = integer(),
                                                                value= numeric()),
                                  genes = data.table::data.table(
                                    gene_id = character(),
                                    gene_name = character()
                                  ),
                                  barcodes = data.table::data.table(
                                    barcodes = character()
                                  )
                 ))

setMethod(f = "show",
          signature = "dbsc",
          definition = function(object){

            print(
              slot(object,'data')[ , .(genes =  slot(object,'genes')[i,gene_id],
                                       barcodes = slot(object,'barcodes')[j,barcodes],
                                       counts = value)]
            )
          })


#' `dbsc2` data class
#'
#' @slot data data.table.
#' @slot genes data.table.
#' @slot barcodes data.table.
#'
#' @return
#' @export
#'
#' @examples
dbsc2 <- setClass(Class = "dbsc2",
                 slots = c(data = 'data.table',
                           genes = 'data.table',
                           barcodes = 'data.table'
                 ),
                 prototype = list(data = data.table::data.table(gene_id = character(),
                                                                barcodes = character(),
                                                                value = numeric(),
                                                                key = c('gene_id',
                                                                        'barcodes')),
                                  genes = data.table::data.table(
                                    gene_id = character(),
                                    gene_name = character(),
                                    key = 'gene_id'
                                  ),
                                  barcodes = data.table::data.table(
                                    barcodes = character(),
                                    key = 'barcodes'
                                  )
                 ))



# setMethod(f = "show",
#           signature = "dbsc2",
#           definition = function(object){
#
#             print(
#               slot(object,'data')[ , .(genes =  slot(object,'genes')[i,gene_id],
#                                        barcodes = slot(object,'barcodes')[j,barcodes],
#                                        counts = value)]
#             )
#           })

