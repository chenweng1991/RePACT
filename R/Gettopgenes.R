#' Gettopgenes
#'
#' This function is based on seurat function to generate marker genes for each cluster
#' @param object  The seurat object that has been analyzed.
#' @param number  number of genes fro each cluster
#' @return  generate a pdf,  if  doreturn=T then return the figures
#' @export
#' @examples
#' Gettopgenes(XX.ob,number)

Gettopgenes<-function(object,number)
{
  pbmc.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  pbmc.markers %>% group_by(cluster) %>% top_n(number, avg_diff) -> topgene
  return(topgene)
}
