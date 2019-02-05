#' Fullplot_v2
#'
#' This function plot basic clustering result for seurat object
#' @param object  The seurat object that has been analyzed.
#' @param name  the name for pdf file
#' @param topgene An optional input of cluster marker genes if "Gettopgenes" has been run outside. This can help to save some time
#' @param resolusion This is clustering resolution used in clustering, which can be found in object@data.info
#' @param signiture  a vector of gene names that are intersting to be checked on tsne by color intensity
#' @param p.tsnecluster if true,  plot the tsne figure stained by clusters
#' @param Pheatmap  if true,  plot the marker gene heatmap
#' @param p.tsnesample  if true,  plot the tsne figure stained by samples
#' @param p.pcacluster  if true,  plot the pca12 figure stained by clusters
#' @param p.pcasample  if true,  plot the pca12 figure stained by sample
#' @param darwPCdrive  if true,  plot the heatmap showing the pca driving genes
#' @param cell.use  the cell number used to draw PCA driving gene plot.  Default is 500
#' @param color  The color palette used for tsne and pca, default is "Paired"
#' @param PC34    If ture, also plot PCA34 additionally, default is FALSE
#' @param dotsize   The dot size for tsne and pca
#' @param doreturn   If ture, retun plot to a variable
#' @param pdfwidth  The width of pdf file.  Default is 7
#' @param pdfheight  The height of pdf file.  Default is 7
#' @return  generate a pdf,  if  doreturn=T then return the figures
#' @export
#' @examples
#'Fullplot_v2(S7rock_1.ob,"S7rock_1.ob.pdf",topgene=NULL,resolusion="res.0.6",signiture=c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A","DNAJB1","GHRL"))
Fullplot_v2<-function(object,name,topgene=NULL,resolusion="res.0.6",signiture=c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A","DNAJB1","GHRL"),p.tsnecluster=T,Pheatmap=T,p.tsnesample=T,p.pcacluster=T,p.pcasample=T,darwPCdrive=T,cell.use=500,color="Paired",PC34=F,dotsize=0.1,heatmapannosize=0.5,doreturn=F,pdfwidth=7,pdfheight=7){
require(RColorBrewer)
require(Seurat)
	titlename<-paste(unique(object@data.info$Sample),collapse="_")
	if(Pheatmap==T)
	{
		if(length(unique(object@data.info[,resolusion]))==1)
		{
			nocluster=TRUE
		}else
		{
			nocluster=FALSE
			if(is.null(topgene))
			{
				topgene<-Gettopgenes(object,10)
			}
		}
	}
	cellPCAstate_ALL<-GetPCAcelldata_v2(object)
	#pbmc.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
	#pbmc.markers %>% group_by(cluster) %>% top_n(5, avg_diff) -> topgene
	getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
	data.infoandtSNE<-merge(object@tsne.rot,object@data.info,by="row.names")
	data.infoandtSNE$Sample<-factor(data.infoandtSNE$Sample)
	if(darwPCdrive)
	{
		pc1.driver<-GetPCheat(object,1,cell=cell.use)+ggtitle("PC1")
		pc2.driver<-GetPCheat(object,2,cell=cell.use)+ggtitle("PC2")
		pc3.driver<-GetPCheat(object,3,cell=cell.use)+ggtitle("PC3")
		pc4.driver<-GetPCheat(object,4,cell=cell.use)+ggtitle("PC4")
	}
	fig1<-ggplot(data.infoandtSNE)+aes(tSNE_1,tSNE_2,color=eval(parse(text=resolusion)))+geom_point(size=dotsize)+theme_bw()+coord_fixed()+scale_color_brewer(palette=color,guide_legend(title="Cell clustering"))+guides(color=guide_legend(override.aes=list(size=5)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename)
	fig2<-ggplot(data.infoandtSNE)+aes(tSNE_1,tSNE_2,color=Sample)+geom_point(size=dotsize)+theme_bw()+coord_fixed()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename)   # page 3
	fig3<-ggplot(cellPCAstate_ALL)+aes(PC1,PC2,color=eval(parse(text=resolusion)))+geom_point(size=dotsize)+theme_bw()+theme_bw()+coord_fixed()+scale_color_brewer(palette=color,guide_legend(title="Cell clustering"))+guides(color=guide_legend(override.aes=list(size=5)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename)
	fig4<-ggplot(cellPCAstate_ALL)+aes(PC1,PC2,color=Sample)+geom_point(size=dotsize)+theme_bw()+theme_bw()+coord_fixed()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename)
	fig5<-ggplot(cellPCAstate_ALL)+aes(PC3,PC4,color=eval(parse(text=resolusion)))+geom_point(size=dotsize)+theme_bw()+theme_bw()+coord_fixed()+scale_color_brewer(palette=color,guide_legend(title="Cell clustering"))+guides(color=guide_legend(override.aes=list(size=5)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename)
	fig6<-ggplot(cellPCAstate_ALL)+aes(PC3,PC4,color=Sample)+geom_point(size=dotsize)+theme_bw()+theme_bw()+coord_fixed()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename)
	fig.list<-list()
	pdf(name,width=pdfwidth,height=pdfheight)
	if(p.tsnecluster)
	{
		print(fig1)
		fig.list<-c(fig.list,list(fig1))
	}
	if(Pheatmap)
	{
		if (!nocluster)
		{
			DoHeatmap(object, genes.use = topgene$gene, order.by.ident = TRUE, col.use = rev(getPalette(10)),slim.col.label = TRUE, remove.key = TRUE, cexRow=heatmapannosize,cex.col=heatmapannosize)  # page2
		}
	}
	if(p.tsnesample)
	{
		print(fig2)
		fig.list<-c(fig.list,list(fig2))
	}
	if(p.pcacluster)
	{
		print(fig3)
		fig.list<-c(fig.list,list(fig3))
	}
	if(p.pcasample)
	{
		print(fig4)
		fig.list<-c(fig.list,list(fig4))
	}
	if(PC34)
	{
		print(fig5)
		print(fig6)
		fig.list<-c(fig.list,list(fig5))
		fig.list<-c(fig.list,list(fig6))
	}

	if(darwPCdrive)
	{
		grid.arrange(pc1.driver,pc2.driver,pc3.driver,pc4.driver,ncol=2)      #page5
		fig.list<-c(fig.list,list(pc1.driver,pc2.driver,pc3.driver,pc4.driver))
	}

	if (!is.null(signiture))
	{
		sig.plot<-GettsnesignatureSuper(object,object,signiture=signiture,toreturn=T)
    cycles<-as.integer(length(sig.plot)/6)
    for (i in 0:(cycles-1)){
      grid.arrange(grobs=sig.plot[(6*i+1):(6*i+6)])
    }
    grid.arrange(grobs=sig.plot[(6*cycles+1):length(sig.plot)])
	}
	dev.off()
	if(doreturn)
	{
		return(fig.list)
	}
}


#' GetPCheat
#'
#' This function plot heatmap for PCA driving genes
#' @param object  The seurat object that has been analyzed.
#' @param pc  the PC#  to plot.  input could be 1 or 2 or 3 ...
#' @param cell  A number of cells to be used for plotting
#' @return  Return the figure
#' @export
#' @examples
#'pc1.driver<-GetPCheat(object,1,cell=cell.use)+ggtitle("PC1")

GetPCheat<-function(object,pc,cell=500){
PCheatdata<-t(PCHeatmap(object, pc.use = pc, cells.use =cell , do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, do.return=TRUE))
PCheatdata.m<-reshape2::melt(PCheatdata)
#PCheatdata.m<-PCheatdata.m[order(PCheatdata.m$value,decreasing=F),]
names(PCheatdata.m)[c(1,2)]<-c("Cells","Genes")
p<-ggplot(PCheatdata.m)+aes(Cells,Genes,fill=value)+geom_tile()+scale_fill_gradient2(low="darkgreen",high="red")+theme(axis.text.x=element_blank(),axis.text.y=element_text(size=6))
return(p)
}
