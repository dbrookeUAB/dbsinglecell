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


