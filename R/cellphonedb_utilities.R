#' CellPhoneDB Summary File
#'
#' @param path the directory containing the CellPhoneDB Output
#' @param pvalue setting this will return results less than it
#'
#' @return
#' @export
#'
#' @examples
#' @import data.table
#' @import Matrix
#' @import crayon
#'
cellphonedb_summary <- function(path, pvalue = 'all'){
  means <- data.table::fread(file.path(path,'means.txt'))
  pvalues <- data.table::fread(file.path(path, 'pvalues.txt'))
  id.vars <- colnames(means)[1:11]

  cat(crayon::green(paste0('\n',Sys.time(),'|',crayon::white(' Reading CellPhoneDB files'))))
  means <- data.table::melt(means, id.vars = id.vars, variable.name = 'cell_pair', value.name = 'mean')
  pvalues <- data.table::melt(pvalues, id.vars = id.vars, variable.name = 'cell_pair', value.name = 'pvalue')

  data.table::setkeyv(means, c('cell_pair',id.vars))
  data.table::setkeyv(pvalues, c('cell_pair',id.vars))

  cat(crayon::green(paste0('\n',Sys.time(),'|',crayon::white(' Merging datasets'))))
  result <- data.table::merge.data.table(means, pvalues)
  result <- as.data.table(result)

  cat(crayon::green(paste0('\n',Sys.time(),'|',crayon::white(' Capturing Gene Pairs'))))
 int_pairs <- strcapture('(.+)\\_(.+)',result$interacting_pair,
            data.table::data.table(gA = character(),
                                   gB = character()))

  cat(crayon::green(paste0('\n',Sys.time(),'|',crayon::white(' Capturing Cell Pairs'))))
 cell_pair <- strcapture('(.+)\\|(.+)',
                         result$cell_pair,
                         data.table::data.table(cell_a = character(),
                                                      cell_b = character()))

 result <- data.table(cell_pair, int_pairs, result)
 if(pvalue=='significant'){
   result <- result[pvalue<0.05]
 }
  cat(crayon::green(paste0('\n',Sys.time(),crayon::yellow('| Finished'))))

  return(result)
}


# prep_cellphonedb <- function(rds, meta_column, path){
#   require(data.table)
#   require(Seurat)
#
#   object <- readRDS(rds)
#
#   res <- sparse2DT.Seurat(object)
#
#   new.meta <- object@meta.data[,meta_column]
#   names(new.meta) <- rownames(object@meta.data)
#
#   # add cell_types to res
#   res[,cell_subset:=new.meta[res$Cell]]
#   data.table::setkey(res, Genes, cell_subset)
#
#   # generate summary information to be used for filtering uninformative genes
#   test <- res[,.(disp = var(Count)/mean(Count), N = .N), c('Genes',meta_column)]
#   test[,total:=sum(N),Genes]
#   test <- test[total>500&!grepl('^mt-',Genes)&!is.na(disp)]
#
#   # create vector with leftover genes
#   gl <- unique(test$Genes)
#
#   # subset count dataset
#   res <- res[Genes %in% gl]
#
#   # create counts file
#   counts <- dcast(res, Genes~Cell, value.var = 'Count', fill = 0)
#   colnames(counts)[1] <- 'Gene'
#   setkey(counts, Gene)
#
#   m2h <- fread('/data/user/dbrooke/db/CellPhoneDB/data/mouse2human.csv', key = 'mouse')
#   mz_genes <- m2h$Ensembl_gene_id
#   names(mz_genes) <- m2h$mouse
#   new_genes <- mz_genes[counts$Gene]
#   names(new_genes) <- counts$Gene
#   new_genes <- new_genes[!is.na(new_genes)]
#
#   dim(counts)
#   counts <- counts[Gene  %in% names(new_genes)]
#   counts[,Gene:=new_genes[Gene]]
#
#   # create meta file
#   meta <- data.table(Cell = colnames(counts)[-1],cell_type =  new.meta[colnames(counts)[-1]])
#   meta <- meta[Cell %in% colnames(counts)[-1]]
#   fwrite(meta, 'PerNiche_int/meta.csv', quote = FALSE)
#
#   fwrite(counts,'PerNiche_int/counts.csv', nThread = 20, showProgress = TRUE)
# }

#' Convert a sparse matrix to a data.table
#'
#' @param matrix  sparse matrix to be used
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
sparse2DT <- function(matrix){

  # creating i,j,x format
  mm.sum <- Matrix::summary(matrix)

  # creating workable dataset of count data
  result <- data.table::data.table(Genes = rownames(matrix)[mm.sum$i], Cell = colnames(matrix)[mm.sum$j], Count = mm.sum$x)
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
