#' docluster
#'
#' This function is to do clustering analysis for one dataset, given a clustering resolusion
#' @param Istopcells  The formatted dge from dgepreprocess
#' @param mylist the informative gene list from GetinformativeGene
#' @param name a customized name for this clustering object
#' @param reso the clustering resolusion, default=0.6
#' @return this will return the clean and square plot
#' @import Seurat ggplot2 Matrix RColorBrewer
#' @export
#' @examples
#' S7rock_1.ob<-docluster(dgepreprocess(s7.RockII_1.dge,500,norowname=T),GetinformativeGene(dgepreprocess(s7.RockII_1.dge,500,norowname=T),500),"s7.RockII_1",reso=0.6)

docluster<-function(Istopcells,mylist,name,reso=0.6)
{
require(Seurat)
require(ggplot2)
require(Matrix)
require(RColorBrewer)
data<-Matrix(as.matrix(Istopcells))
pbmc <- new("seurat", raw.data = data)
pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, do.logNormalize = T, total.expr = 1e4, project = "10X_PBMC")
pbmc@var.genes<-row.names(x=mylist)
mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
if(length(mito.genes)!=0){
	percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))
}else{
	percent.mito<-rep(0,ncol(pbmc@data))
	names(percent.mito)<-colnames(pbmc@data)
}
pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
basic<-VlnPlot(pbmc, c("nGene", "nUMI", "percent.mito"),size.use=0.3,do.ret=TRUE)
pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = 2500)
pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.05)
pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))
pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
pbmc <- ProjectPCA(pbmc)
pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = reso, print.output = 0, save.SNN = T)
pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = (ncol(Istopcells)>200))
pbmc@data.info$Sample<-name
return(pbmc)
}




#' dgepreprocess
#'
#' This function is to preprocess raw dge into a format for clustering analysis, including add rownames and cutoff low expressing cells
#' @param isletdge  The raw dge data
#' @param Txcutoff The UMI cutoff as for defined "cell"
#' @param norowname The input dge file with or without row names? if true, input already has rowname(usually when data taken from other analyzed Seurat object), if FALSE, the input without rowname and ususally is the very raw data
#' @return this will return the formated dge for clustering
#' @export
#' @examples
#' dgepreprocess(s7.RockII_1.dge,500,norowname=T)

dgepreprocess<-function(isletdge,Txcutoff=1000,norowname=T){
if(norowname)
{
	rownames(isletdge)<-isletdge[,1]
	isletdge<-isletdge[,-1]
}else
{
	isletdge<-as.data.frame(isletdge)
}
isletdge<-rbind(isletdge,Sum=apply(isletdge,2,sum))
Istopcells<-isletdge[,which(isletdge["Sum",]>Txcutoff)]
Istopcells<-Istopcells[1:(nrow(Istopcells)-2),]
Istopcells<-cbind(Istopcells,genesum=apply(Istopcells,1,sum))
Istopcells<-subset(Istopcells,genesum!=0)
Istopcells<-Istopcells[,1:ncol(Istopcells)-1]
return(Istopcells)
}





#' GetinformativeGene
#'
#' This function is to identify the most informative genes for clustering analysis
#' @param Istopcells  The formatted dge from dgepreprocess
#' @param Number number of genes for out put
#' @param cutoff the percentile cutoff for defining a variable gene, default=0.95
#' @return this will return a vector of informative genes for clustering
#' @export
#' @examples
#' GetinformativeGene(dgepreprocess(s7.RockII_1.dge,500,norowname=T),500)
GetinformativeGene<-function(Istopcells,Number,cutoff=0.95){
require(dplyr)
#require(raster)
data<-as.data.frame(Istopcells)
expvalue<-function(vector){
x<-5000*vector/sum(vector[1:length(vector)-1])
return(x)
}
CV<-function(vector){
x<-(sd(vector))^2/mean(vector)
return(x)
}
GeneMean<-apply(data,1,mean)
Genecv<-apply(data,1,cv)
GeneCV<-apply(data,1,CV)
Genesd<-apply(data,1,sd)
data.meran_cv_sd<-cbind(data,GeneMean,Genecv,GeneCV,Genesd)
data.meran_cv_sd<-data.meran_cv_sd[order(data.meran_cv_sd$GeneMean,decreasing=T),]
decision<-rep(F,10)
if(dim(Istopcells)[1]>=1000)
{
	top1000<-data.meran_cv_sd[1:1000,]
	list1<-top1000[which(top1000$GeneCV>quantile(top1000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list1)[1]>1)
	{
	decision[1]<-T
	}
}
if(dim(Istopcells)[1]>=2000)
{
	top2000<-data.meran_cv_sd[1001:2000,]
	list2<-top2000[which(top2000$GeneCV>quantile(top2000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list2)[1]>1)
	{
	decision[2]<-T
	}
}
if(dim(Istopcells)[1]>=3000)
{
	top3000<-data.meran_cv_sd[2001:3000,]
	list3<-top3000[which(top3000$GeneCV>quantile(top3000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list3)[1]>1)
	{
	decision[3]<-T
	}
}
if(dim(Istopcells)[1]>=4000)
{
	top4000<-data.meran_cv_sd[3001:4000,]
	list4<-top4000[which(top4000$GeneCV>quantile(top4000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list4)[1]>1)
	{
	decision[4]<-T
	}
}
if(dim(Istopcells)[1]>=5000)
{
	top5000<-data.meran_cv_sd[4001:5000,]
	list5<-top5000[which(top5000$GeneCV>quantile(top5000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list5)[1]>1)
	{
	decision[5]<-T
	}
}
if(dim(Istopcells)[1]>=6000)
{
	top6000<-data.meran_cv_sd[5001:6000,]
	list6<-top6000[which(top6000$GeneCV>quantile(top6000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list6)[1]>1)
	{
	decision[6]<-T
	}
}
if(dim(Istopcells)[1]>=7000)
{
	top7000<-data.meran_cv_sd[6001:7000,]
	list7<-top7000[which(top7000$GeneCV>quantile(top7000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list7)[1]>1)
	{
	decision[7]<-T
	}
}
if(dim(Istopcells)[1]>=8000)
{
	top8000<-data.meran_cv_sd[7001:8000,]
	list8<-top8000[which(top8000$GeneCV>quantile(top8000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list8)[1]>1)
	{
	decision[8]<-T
	}
}
if(dim(Istopcells)[1]>=9000)
{
	top9000<-data.meran_cv_sd[8001:9000,]
	list9<-top9000[which(top9000$GeneCV>quantile(top9000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list9)[1]>1)
	{
	decision[9]<-T
	}
}
if(dim(Istopcells)[1]>=10000)
{
	top10000<-data.meran_cv_sd[9001:10000,]
	list10<-top10000[which(top10000$GeneCV>quantile(top10000$GeneCV,cutoff)),1:(ncol(data.meran_cv_sd)-4)]
	if(dim(list10)[1]>1)
	{
	decision[10]<-T
	}
}


if(all(decision[1:10]))
{
mylist.un.500<-rbind(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10)
}
if(all(decision[1:8]))
{
mylist.un.400<-rbind(list1,list2,list3,list4,list5,list6,list7,list8)
}
if(all(decision[1:6]))
{
mylist.un.300<-rbind(list1,list2,list3,list4,list5,list6)
}
if(all(decision[1:4]))
{
mylist.un.200<-rbind(list1,list2,list3,list4)
}
if(all(decision[1:2]))
{
mylist.un.100<-rbind(list1,list2)
}

if (length(which(decision))==10)
{
return(mylist.un.500)
}else if(length(which(decision))>=8)
{
return(mylist.un.400)
}else if(length(which(decision))>=6)
{
return(mylist.un.300)
}else if(length(which(decision))>=4)
{
return(mylist.un.200)
}else if(length(which(decision))>=4)
{
return(mylist.un.100)
}
}


#' docluster.multi
#'
#' This function is to do clustering analysis for a number of sets combined
#' @param Number  the number of informative genes for clustering, default is 500
#' @param txcutoff the informative gene list from GetinformativeGene, default is 500
#' @param sets a list of dge filesw
#' @param nms  a vector of names for each sample
#' @param selected  NULL
#' @param filterstuff NULL
#' @param reso  The resolusion for clustering,  default is 0.6
#' @param israw,  the dge in sets is from raw raw data that has no rowname if TRUE, or if FALSE it is extracted raw matrix from Seurat object
#' @param dofindcluster if TRUE, do find clustering, default is true
#' @param dotsne  if TRUE, do tsne, default is true
#' @param acphi  accepted hiest gene number per cell, default is 2500(Should increase when processing fludigam data)
#' @return this will return the clean and square plot
#' @export
#' @examples
#' S7rock_1.S7<-docluster.multi(500,txcutoff=500,sets=list(s7.RockII_1=s7.RockII_1.dge,s7.B=s7.B.dge),nms=c("s7.RockII_1","s7.B"),israw=T)

docluster.multi<-function(Number,txcutoff=500,sets,nms,selected=NULL,filterstuff=NULL,reso=0.6,israw=F,dofindcluster=T,dotsne=T,acphi=2500)
{
	require(Seurat)
	require(ggplot2)
	require(Matrix)
	require(RColorBrewer)
	if(israw)
	{
		for (i in 1: length(sets))
		{
			row.names(sets[[i]])<-sets[[i]][,1]
			sets[[i]]<-sets[[i]][,-1]
		}
	}
	Mergedset<-sets[[1]][,colSums(sets[[1]])>txcutoff]
	if(length(sets)<2)
	{
	Mergedset<-as.matrix(Mergedset)
	}else
	{
		for (i in 2: length(sets))
		{
			sets[[i]]<-sets[[i]][,colSums(sets[[i]])>txcutoff]
			colnames(sets[[i]])<-paste(colnames(sets[[i]]),letters[i],sep="_")
			Mergedset<-Tomerge(Mergedset,sets[[i]])
		}
	}
	if (!is.null(selected))
	{
		Mergedset<-Mergedset[,colnames(Mergedset) %in% selected]
	}
	if (!is.null(filterstuff))
	{
#		filterstuff<-filterstuff[which(nchar(filterstuff)==12)]
		Mergedset<-Mergedset[,filterstuff[filterstuff %in% colnames(Mergedset)]]
	}
	Mergedset[is.na(Mergedset)]<-0
	data<-Matrix(as.matrix(Mergedset))
	mylist<-GetinformativeGene(Mergedset,Number)
	pbmc <- new("seurat", raw.data = data)
	pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, do.logNormalize = T, total.expr = 1e4, project = "10X_PBMC")
	pbmc@var.genes<-row.names(x=mylist)
	mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
	if (length(mito.genes)!=0)
	{
	percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))
	}else{
	percent.mito<-rep(0,length(colnames(pbmc@data)))
	names(percent.mito)<-colnames(pbmc@data)
  }

	pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
	donor<-c()
	for (i in 1: length(sets))
	{
	donor[match(colnames(sets[[i]]),row.names(pbmc@data.info))]<-nms[i]
	}
	pbmc@data.info$Sample<-donor
	pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = acphi)
	pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.05)
	pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))
	pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
	pbmc <- ProjectPCA(pbmc)
	if(dofindcluster)
	{
		pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = reso, print.output = 0, save.SNN = T)
	}
	if(dotsne)
	{
	pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = T)
	data.infoandtSNE<-merge(pbmc@tsne.rot,pbmc@data.info,by="row.names")
	data.infoandtSNE$Sample<-factor(data.infoandtSNE$Sample)
	}
	return(pbmc)
}


#' docluster.single
#'
#' This function a simple function to do clustering for one dge file
#' @param Number number of the most informative genes for clustering, default is 500
#' @param Mergedset  The dge matrix
#' @param nm1 the name for this sample
#' @param geneminThe the minimum gene number as hard filtering for a valid cell
#' @param cellmin the minimum cell number that express a specific gene as  hard filtering for an expressed gene
#' @return this will return a serurat object for other analysis
#' @export
#' @examples
#' GetinformativeGene(dgepreprocess(s7.RockII_1.dge,500,norowname=T),500)

docluster.single<-function(Number,Mergedset,nm1="",dict=NULL,reso=0.6,genemin=200,cellmin=3)
{
require(Seurat)
require(ggplot2)
require(Matrix)
require(RColorBrewer)

Mergedset[is.na(Mergedset)]<-0
if(length(which(rowSums(Mergedset)==0)))
{
Mergedset<-Mergedset[-which(rowSums(Mergedset)==0),]
}
data<-Matrix(as.matrix(Mergedset))
mylist<-GetinformativeGene(Mergedset,Number)
pbmc <- new("seurat", raw.data = data)
pbmc <- Setup(pbmc, min.cells = cellmin, min.genes = genemin, do.logNormalize = T, total.expr = 1e4, project = "10X_PBMC")
pbmc@var.genes<-row.names(x=mylist)
mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))
pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
if (!is.null(dict))
  {pbmc@data.info<-Tomerge_v2(pbmc@data.info,dict)
  } else
  {pbmc@data.info<-cbind(pbmc@data.info,Sample=nm1)
  }
pbmc@data.info$nGene<-as.numeric(as.character(pbmc@data.info$nGene))
pbmc@data.info$nUMI<-as.numeric(as.character(pbmc@data.info$nUMI))
pbmc@data.info$percent.mito<-as.numeric(as.character(pbmc@data.info$percent.mito))

pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = 2500)
pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.05)
pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))
pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
pbmc <- ProjectPCA(pbmc)
pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = reso, print.output = 0, save.SNN = T)
pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = T)
data.infoandtSNE<-merge(pbmc@tsne.rot,pbmc@data.info,by="row.names")
data.infoandtSNE$Sample<-factor(data.infoandtSNE$Sample)
return(pbmc)
}
