f_list <- gl3


pb <- progress::progress_bar$new(
  format = " :what [:bar] :elapsed :percent eta: :eta",
  clear = FALSE, total = length(f_list), width = 50)

f <- function(x){
  pb$tick(tokens = list(what = x))
  Sys.sleep(2 / 100)
  result <- miknn::mi_knn(blah2[gene_name==x], 'cell_type', 'Count', warnings = FALSE,nthreads = 3, global = FALSE, k = 5)
  result <- dcast(result, .~cell_type, value.var = 'I')
  result[1,1] <- x
  colnames(result)[1] <- 'gene_name'
  return(result)
}


 res <- do.call(rbind,lapply(f_list, f))



