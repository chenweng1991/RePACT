#' GettsnesignatureSuper
#'
#' This function plot tsne or PCA stained by gene relative expression
#' @param object  The seurat object that has been analyzed.
#' @param object.all   The total object for normalization.  Usually it is the same as   object
#' @param signiture  a vector of gene names that are intersting to be checked on tsne by color intensity
#' @param doPCA  if ture then stain on PCA12 plot.  default is F
#' @param dotSNE   if ture then stain on tsne plot.  default is T
#' @param usePC34 if true,  plot the PC34 istead of PC12. Default is F
#' @param extratitle  A string that can be added to the title of figures. default is ""
#' @param dotsize   The dot size for the figures,  default is 0.3
#' @param buttomgrey  if true,  The cells with no expression set color as grey
#' @param nolegend  if true,  no color intensity legend is put
#' @param highcolor  The color for highest expression, A color spectrum is formed between "grey" to the color set. Default is red
#' @return  draw figures. if toreturn=T, then return a list of figures
#' @export
#' @examples
#'Fullplot_v2(S7rock_1.ob,"S7rock_1.ob.pdf",topgene=NULL,resolusion="res.0.6",signiture=c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A","DNAJB1","GHRL"))

GettsnesignatureSuper<-function(object,object.all,signiture=c("INS","RBP4","FFAR4","ID1","ID2","ID3","DNAJB1"),doPCA=F,dotSNE=T,usePC34=F,extratitle="",toreturn=F,dotsize=0.3,buttomgrey=T,nolegend=T,highcolor="red")
{
PC1used<-"PC1"
PC2used<-"PC2"
if(usePC34)
{
	PC1used<-"PC3"
	PC2used<-"PC4"
}

	geneNA<-signiture[which(!signiture %in% row.names(object.all@scale.data))]
	signiture.expre<-setdiff(signiture,geneNA)
	tSNE_signature<-as.data.frame(Tomerge(object@tsne.rot,t(object.all@scale.data[signiture.expre,])))
	for (gene in geneNA)
	{
		tSNE_signature<-cbind(gene=0,tSNE_signature)
		names(tSNE_signature)[1]<-gene
	}

	toPlottSNEandPCA<-Tomerge(tSNE_signature,GetPCAcelldata_v2(object)[,c(1,2,3,4)])
	all.list<-list()
	if(dotSNE)
	{
		for (i in 1:length(signiture))
		{
			genename<-signiture[i]
			if(genename %in% geneNA)
			{
				lowlim<-0
				highlim<-0
			}else
			{
				lowlim<-min(t(object.all@scale.data[signiture.expre,])[,genename])
				highlim<-max(t(object.all@scale.data[signiture.expre,])[,genename])
			}
			if(buttomgrey)
			{
				if(genename %in% geneNA)
				{
				p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)+scale_color_gradient(low="grey",high="grey",limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)
				}else
				p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)
			}else
			{
				p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)
			}
			p<-format2(p,toPlottSNEandPCA,"tSNE_1","tSNE_2",nolegend=nolegend)
			all.list[[i]]<-p
			#print(p)
		}
	}
	if(doPCA)
	{
		for (i in 1:length(signiture))
		{
			genename<-signiture[i]
			if(genename %in% geneNA)
			{
				lowlim<-0
				highlim<-0
			}else
			{
				lowlim<-min(t(object.all@scale.data[signiture.expre,])[,genename])
				highlim<-max(t(object.all@scale.data[signiture.expre,])[,genename])
			}
			if(buttomgrey)
			{
				if(genename %in% geneNA)
				{
				p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string(PC1used,PC2used,color=genename),size=dotsize)+scale_color_gradient(low="grey",high="grey",limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"PC_space",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string(PC1used,PC2used,color=genename),size=dotsize)
				}else
				p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string(PC1used,PC2used,color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"tSNE_space",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string(PC1used,PC2used,color=genename),size=dotsize)
			}else
			{
				p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string(PC1used,PC2used,color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"tSNE_space",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string(PC1used,PC2used,color=genename),size=dotsize)
			}
			p<-format2(p,toPlottSNEandPCA,PC1used,PC2used)
			#print(p)
			all.list<-c(all.list,list(p))
		}
	}

	if(toreturn)
	{
		return(all.list)
	}
}
