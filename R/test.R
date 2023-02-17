#
#
# dbsc.obj <- setClass(Class = "dbsc.obj",
#                   slots = c(data = 'data.table',
#                             gene_info = 'data.table',
#                             barcode_info = 'data.table'
#                   ),
#                   prototype = list(data = data.table::data.table(gene_id = character(),
#                                                                  barcodes = character(),
#                                                                  count = numeric(),
#                                                                  key = c('gene_id',
#                                                                          'barcodes')),
#                                    gene_info = data.table::data.table(
#                                      gene_id = character(),
#                                      gene_name = character(),
#                                      key = 'gene_id'
#                                    ),
#                                    barcode_info = data.table::data.table(
#                                      barcodes = character(),
#                                      key = 'barcodes'
#                                    )
#                   ))
#
# dbsc.obj2 <- setClass(Class = "dbsc.obj2",
#                      slots = c(data = 'data.table',
#                                gene_info = 'data.table',
#                                barcode_info = 'data.table'
#                      ),
#                      prototype = list(data = data.table::data.table(id = character(),
#                                                                     gene_id = character(),
#                                                                     barcodes = character(),
#                                                                     count = numeric(),
#                                                                     key = c('gene_id',
#                                                                             'barcodes')),
#                                       gene_info = data.table::data.table(
#                                         gene_id = character(),
#                                         gene_name = character(),
#                                         key = 'gene_id'
#                                       ),
#                                       barcode_info = data.table::data.table(
#                                         id = character(),
#                                         barcodes = character(),
#                                         key = 'barcodes'
#                                       )
#                      ))
#
# db_create_object <- function(counts, gene_info, barcode_info){
#
# # general tests -------------------------------------------------------------------------------
# test_dt <- sapply(list(counts = counts,
#                        barcode_info = barcode_info,
#                        gene_info = gene_info),
#                   data.table::is.data.table)
#
# if(!all(test_dt)){
#   stop(paste0('objects must be data.tables: ',
#          paste(names(which(!test_dt)), collapse = ', ')))
# }
#
# if(!identical(colnames(counts),c('gene_id','barcodes','count'))){
#   stop('counts object must have three columns named: gene_id, barcodes, count')
# }
#
# # gene_info tests -----------------------------------------------------------------------------
#
#   if(!'gene_id' %in% colnames(gene_info)){
#     stop("gene_id column is not found in gene_info")
#   }
#
#
#
# # barcode_info tests --------------------------------------------------------------------------
# if(!'barcodes' %in% colnames(barcode_info)){
#   stop("barcodes column is not found in barcode_info")
# }
#
# # counts tests --------------------------------------------------------------------------------
#
#
#   result <- dbsc.obj(data = counts,
#                      gene_info = gene_info,
#                      barcode_info = barcode_info
#                      )
#   data.table::setkey(result@data, gene_id, barcodes)
#   data.table::setkey(result@gene_info, gene_id)
#   data.table::setkey(result@barcode_info,barcodes)
#   return(result)
# }
#
# setMethod(f = "show",
#           signature = "dbsc.obj",
#           definition = function(object){
#
#             print(
#               slot(object,'data')
#             )})
#
# setMethod(f = "show",
#           signature = "dbsc.obj2",
#           definition = function(object){
#
#             print(
#               slot(object,'data')
#             )})
#
# DuplicateBarcodes <- function(x,y){
#   any(duplicated(c(slot(x,'barcode_info')[['barcodes']],slot(y,'barcode_info')[['barcodes']])))
# }
#
# MergeBarcodes <- function(x,y){
#  tmp <- list(slot(x,'barcode_info'),slot(y,'barcode_info'))
#  names(tmp) <- c(as.character(substitute(x)),as.character(substitute(y)))
#
#  BARCODE_INFO <- data.table::rbindlist(
#    tmp,
#    use.names = TRUE, fill = TRUE,
#    idcol = 'id'
#  )
#  return(BARCODE_INFO)
# }
#
# MergeData <- function(x,y){
#   tmp <- list(slot(x,'data'),slot(y,'data'))
#   names(tmp) <- c(as.character(substitute(x)),as.character(substitute(y)))
#
#   COUNTS <- data.table::rbindlist(
#     tmp,
#     use.names = TRUE, fill = TRUE,
#     idcol = 'id'
#   )
#   return(COUNTS)
# }
#
#
# db_merge_objects <- function(object_list){
#   object_names <- as.character(substitute(object_list))[-1]
#   names(object_list) <- object_names
#   dbsc.obj2(data = data.table::data.table(data.table::rbindlist(
#     lapply(object_list, function(x) slot(x,'data')),
#     use.names = TRUE,
#     fill = TRUE,
#     idcol = 'id'
#     ),key = c('id','gene_id','barcodes')),
#     gene_info = data.table::data.table(unique(data.table::rbindlist(
#       lapply(object_list, function(x) slot(x,'gene_info')),
#       use.names = TRUE,
#       fill = TRUE
#     )),
#     key = 'gene_id'),
#     barcode_info = data.table::data.table(data.table::rbindlist(
#       lapply(object_list, function(x) slot(x,'barcode_info')),
#       use.names = TRUE,
#       fill = TRUE,
#       idcol = 'id'
#     ),
#     key = c('id','barcodes'))
#     )
# }
#
#
# fp <- system.file('extdata/GSE174241/GSM5289774',package = 'dbsinglecell')
# fp2 <- system.file('extdata/GSE174241/GSM5289775',package = 'dbsinglecell')
#
#
# tmp <- dbsinglecell::db_read10x_v6(fp)
# db1 <- db_create_object(counts = tmp[,.(gene_id,barcodes,count)],
#                         gene_info = unique(tmp[,.(gene_id, gene_name, gene_type = 'Gene Expression')]),
#                         barcode_info = unique(tmp[,.(barcodes)]))
#
# tmp2 <- dbsinglecell::db_read10x_v6(fp2)
# db2 <- db_create_object(counts = tmp2[,.(gene_id,barcodes,count)],
#                         gene_info = unique(tmp2[,.(gene_id, gene_name, gene_type = 'Gene Expression')]),
#                         barcode_info = unique(tmp2[,.(barcodes)]))
#
# tmp3 <- dbsinglecell::db_read10x_v6('~/rdx/db/scRNA/mouse/datasets/GSE199563/GSM5975235')
# db3 <- db_create_object(counts = tmp3[,.(gene_id,barcodes,count)],
#                         gene_info = unique(tmp3[,.(gene_id, gene_name, gene_type = 'Gene Expression')]),
#                         barcode_info = unique(tmp3[,.(barcodes)]))
#
#
#
#
#
#
#
#
#
#
#
