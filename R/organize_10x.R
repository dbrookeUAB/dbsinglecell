#' Reorganize another person's mess into a usable 10X dataset
#'
#' @param x  path containing the unorganized disaster
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples
organize_10x <- function(x ){
  path_main <- x
  file_list <- dir(path = x, full.names = T)
  file_list <- file_list[grepl('tsv.gz$',file_list)|grepl('mtx.gz$',file_list)]

  res <- data.table::data.table(strcapture('(GSM\\d+)_.+([fmbg][ae][tarn]\\w+.\\w{3}.gz)', x = basename(file_list),
                               proto = data.table::data.table(accession_id = character(),
                                                  file_type = character())))
  res$old_path <- file_list
  res$old_name<- basename(file_list)
  res[res$file_type=='genes.tsv.gz']$file_type<-'features.tsv.gz'
  res$new_folder <- file.path(path_main,paste0(res$accession_id))
  res$new_path <- file.path(res$new_folder, res$file_type)

  new_dirs <- unique(res$new_folder)

  length(file_list)
  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent eta: :eta",
    clear = FALSE, total = length(file_list), width = 80)


  for(i in new_dirs){

    if(!dir.exists(i)){
      dir.create(i)
    }

    f2m <- res[new_folder==i]

    for(j in 1:nrow(f2m)){
    pb$tick()
      file.copy(f2m[j,old_path ],f2m[j,new_path] )
    }

  }

for(i in file_list){
  file.remove(i)
}

}


#' Convert a flat text file into a 10X-formated folder
#'
#' @param file
#' @param new_folder
#'
#' @return
#' @export
#' @import Matrix
#' @import data.table
#'
#' @examples
flatMatrix_2_10x <- function(file, new_folder = NULL, extract = NULL){

  counts <- as(as.matrix(suppressWarnings( data.table::fread(file)),rownames = 1),'dgCMatrix')

  if(!is.null(extract)){
    new_folder <- file.path(dirname(file),stringr::str_extract(basename(file),'GS[EM]\\d+'))
  } else if(is.null(new_folder)&is.null(extract)){
    new_folder <- basename(dirname(file))
  }

  DropletUtils::write10xCounts(path = new_folder, x = counts, overwrite = TRUE,version = '3')
}

flatMatrix_2_10x_flip <- function(file, new_folder = NULL, extract = NULL){
  message_task('Importing Count Matrix')
  counts <- data.table::fread(file)
  message_append('converting to matrix')
  counts <- as.matrix(counts, rownames = 1)
  message_append('tansposing')
  counts <- t(counts)
  message_append('converting to dgCMatrix')
  counts <- as(counts, 'dgCMatrix')

  # counts <- as((t(as.matrix(suppressWarnings( data.table::fread(file)),rownames = 1))),'dgCMatrix')

  if(!is.null(extract)){
    new_folder <- file.path(dirname(file),stringr::str_extract(basename(file),'GS[EM]\\d+'))
  } else if(is.null(new_folder)&is.null(extract)){
    new_folder <- basename(dirname(file))
  }

  message_task('Writing 10X files')
  DropletUtils::write10xCounts(path = new_folder, x = counts, overwrite = TRUE,version = '3')
}


#' Write 10X Files from Sparse Matrix
#'
#' @param matrix
#' @param path
#' @param gene_id
#' @param gene_symbol
#'
#' @return
#' @export
#'
#' @examples
db_write10X <- function(matrix, path = 'output', gene_id = NULL, gene_symbol = NULL){

  i <- matrix@i+1L
  j <- findInterval(seq_along(matrix@x)-1L,matrix@p[-1])+1L
  x  <- matrix@x

  feats <- matrix@Dimnames[[1]]
  features <- data.table::data.table(gene_id = feats)

  if(is.null(gene_id)){
    features <- data.table::data.table(gene_id = feats)
  } else {
    features <- data.table::data.table(gene_id = gene_id)
  }


  if(is.null(gene_symbol)){
    features$gene_symbol <- feats
  } else {
    features$gene_symbol <- gene_symbol
  }

  features$type <- 'Gene Expression'

  if(is.null(gene_symbol))

    barcodes <- matrix@Dimnames[[2]]



  result <- data.table::data.table(i,j,x)

  if(!dir.exists(path)){
    dir.create(path)
  }

  matrix_file <- file.path(path,'matrix.mtx.gz')
  features_file <- file.path(path, 'features.tsv.gz')
  barcodes_file <- file.path(path, 'barcodes.tsv.gz')

  con <- gzfile(matrix_file)
  writeLines(text = c( "%%MatrixMarket matrix coordinate real general" ,paste(matrix@Dim[1], matrix@Dim[2], length(x))), con =con)
  data.table::fwrite(result, file = matrix_file, append = TRUE,sep = ' ',nThread =  parallel::detectCores()-1)
  close.connection(con)

  con <- gzfile(features_file)

  data.table::fwrite(features, file =features_file, append = TRUE, sep = '\t', col.names = FALSE)

  close.connection(con)

  con <- gzfile(barcodes_file)

  writeLines(barcodes, con = con)
  close.connection(con)
}


#' Create Counts objet
#'
#' @param matrix
#'
#' @return
#' @export
#'
#' @examples
db_counts <- function(matrix){
  result <- data.table::data.table(i =matrix@i+1L,
                                   j =findInterval(seq_along(matrix@x)-1L,matrix@p[-1])+1L,
                                   x = matrix@x)
  class(result) <- c('db_assay','counts',class(result))
  return(result)
}



# db_assays <- setClass('db_assays',
                      # slots = list(counts = 'data.table'),
                      # prototype = prototype(list(counts = data.table::data.table(i = integer(), j = integer(),x = integer()))), contains = 'list')




db_sce <- setClass('db_sce',slots = c(metaData = 'data.table', geneData = 'data.table', assays = 'list'),

                   prototype  = prototype(metaData = data.table::data.table(id = character(), nFeatures = integer(), nCounts = integer()),
                                          geneData = data.table::data.table(gene_id = character(), gene_symbol = character()),
                                          assays = list()
                   )
)


#' Print
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.db_sce <- function(x,...){
  tmp <- c(Class = 'class: db_sce',
           Dim = paste('dim:',nrow(x@geneData),nrow(x@metaData)),
           Assays = paste('assays: ', paste(names(x@assays), collapse = ', ')),
           gene_ids = paste0('gene_ids: ', paste(x@geneData$gene_id[c(1,2)],x@geneData$gene_id[c(ncol(x@geneData)-c(1,0))],collapse = ' .... ')),
           MetaData = paste0('metaData: ',paste(colnames(x@metaData)[-c(1:2)], collapse = ', ')),
           cell_id = paste0('cell_id: ',paste(x@metaData$cell_id[c(1,2)],x@metaData$cell_id[c(ncol(x@geneData)-c(1,0))],collapse = ' .... '))
  )

  cat(paste(tmp, collapse = '\n'))


}

#' Show
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
show.db_sce <-   function(x){
  tmp <- c(Class = 'class: db_sce',
           Dim = paste('dim:',nrow(x@geneData),nrow(x@metaData)),
           Assays = paste('assays: ', paste(names(x@assays), collapse = ', ')),
           gene_ids = paste0('gene_ids: ', paste(x@geneData$gene_id[c(1,2)],x@geneData$gene_id[c(ncol(x@geneData)-c(1,0))],collapse = ' .... ')),
           MetaData = paste0('metaData: ',paste(colnames(x@metaData)[-c(1:2)], collapse = ', ')),
           cell_id = paste0('cell_id: ',paste(x@metaData$cell_id[c(1,2)],x@metaData$cell_id[c(ncol(x@geneData)-c(1,0))],collapse = ' .... '))
  )
  cat(paste(tmp, collapse = '\n'))

}


#' Read10X as db_sce object
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
db_read10X <- function(path){
  file_info <- as.integer(strsplit( readLines(file.path(path,'matrix.mtx.gz'),n = 2)[2],'\\s')[[1]])

  names(file_info) <- c('features','barcodes','length')

  assays <- list(counts = data.table::fread(file.path(path,'matrix.mtx.gz'),
                                                 header = FALSE,
                                                 colClasses = c('integer','integer','integer'),
                                                 skip = 2,
                                                 col.names = c('i','j','counts'),
                                                 nrows = file_info['length']+2,
                                                 nThread = parallel::detectCores()-1 )
  )

  result <- db_sce( assays = assays,
                    metaData = data.table::data.table(cell_id = readLines(file.path(path,'barcodes.tsv.gz')),
                                                      assays$counts[,.(nFeatures = .N, nCounts = sum(counts)) ,j]),

                    geneData = data.table::fread(file.path(path,'features.tsv.gz'),
                                                 header = FALSE,
                                                 col.names = c('gene_id','gene_symbol'),
                                                 select = c('V1','V2')))


  return(result)

}




