#' Read10x v1
#'
#' @param path
#' @param return.sce return result as SingleCellExperiment object
#'
#' @return
#' @export
#' @import data.table
#' @import Matrix
#' @import SingleCellExperiment
#'
#' @examples
read10x  <- function(path, return.sce = TRUE){
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
    header = FALSE)$V2


# duplicate gene names for row names --------------------------------------

  if(!all(duplicated(gene)==FALSE)){
    dg <- data.table::data.table(
      position = which(duplicated(gene)),
      name = gene[duplicated(gene)])[,N:=.N,name][]
    dg[,new.name:=paste0(name,'.',1:.N), name]
    gene[dg$position] <- dg$new.name
  }
 res <-  Matrix::sparseMatrix(
    i =  mat$i,
    j = mat$j,
    x = mat$value,
    dimnames = list(gene,barcode))

 if(return.sce){
   SingleCellExperiment::SingleCellExperiment(list(counts = res), meta = meta)
 } else {
   return(res)
 }

}

#' Read10x v2
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


#' Read10x v3
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
    SingleCellExperiment::SingleCellExperiment(list(counts = res))
  } else {
    return(gene)
  }

}
