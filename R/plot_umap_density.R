plot_umap_density <- function(object, reduction = 'umap',dims = 1:2, label = TRUE, cols = NULL, group.by = NULL, bins = 256){
  Label <- logical(length = 1)
  Label <- as.logical(label)
  # plotting embedding ------------------------------------------------------
  emb <-  data.table::as.data.table(Seurat::Embeddings(object, reduction = reduction)[,dims], keep.rownames = 'id')

  # cluster assignment ------------------------------------------------------
  if(is.null(group.by)){
    emb$cluster <- as.character(Seurat::Idents(object)[emb$id])
  } else {
    emb$cluster <- object@meta.data[ names(Seurat::Idents(object)),group.by]
  }

  # plot parameters and labels ----------------------------------------------
  f <- function(z){
    res <- density(z)
    res$x[which.max(res$y)]
    }

  xy <- colnames(emb[,c((1:2)+1)])
  label <- emb[,.(x =f(get(xy[1])), y = f(get(xy[2]))),cluster]

  # plotting ----------------------------------------------------------------
  suppressWarnings(
    p1 <- ggplot2::ggplot(data = emb,
                          mapping = ggplot2::aes_string(x = xy[1], y = xy[2], fill = 'cluster'))+
      ggplot2::stat_density_2d(geom = "polygon",
                               ggplot2::aes(alpha = ..level.., fill = cluster),
                               contour_var = 'count',n = 10000,
                               bins = bins) +
      ggplot2::guides(alpha = FALSE)+
      ggplot2::theme_linedraw()+
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank())
  )

  if(Label){
    p1 <- p1+ ggplot2::geom_text(data = label, mapping = ggplot2::aes(x, y, label = cluster))+
      ggplot2::theme(legend.position = 'none')
  }


  return(p1)

}




