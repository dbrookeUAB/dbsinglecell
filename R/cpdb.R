#' Simplified importing CellphoneDB Files
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
import_cpdb_results <- function(path){
  message_section('Importing CellphoneDB Files')
  fl <- list.files(path, full.names = TRUE, pattern = '.txt',include.dirs = FALSE)
  means_path <- fl[grepl('means.txt$',basename(fl))&!grepl('significant',basename(fl))]
  pvalue_path <- fl[grepl('pvalues.txt$',basename(fl))]
  message_task('importing means.txt')
  means <- data.table::fread(means_path)
  message_task('importing pvalues.txt')
  pvalues <- data.table::fread(pvalue_path)

  message_task('converting data to long form')
  means <- data.table::melt(means, colnames(means)[1:11], variable.name = 'pair', value.name = 'mean' )
  pvalues <- data.table::melt(pvalues, colnames(pvalues)[1:11], variable.name = 'pair', value.name = 'pvalue' )
  data.table::setkeyv(means, colnames(means)[1:12])
  data.table::setkeyv(pvalues, colnames(pvalues)[1:12])

  message_task('combining means.txt with pvalues.txt')
  result <- means[pvalues]

  message_task('creating new columns: cell_a, cell_b')
  result[,c('cell_a','cell_b'):=strcapture('(^\\w.+)\\|(\\w.+$)',pair, data.frame(cell_a = character(), cell_b = character()))][]
  message_task('creating new columns: protein_a, protein_b')
  result[,c('protein_a','protein_b'):=strcapture('(^\\w.+)\\_(\\w.+$)',interacting_pair, data.frame(protein_a = character(), protein_b = character()))][]

  message_task('returning result')
  return(result)

}
