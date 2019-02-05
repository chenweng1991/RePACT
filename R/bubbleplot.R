#' Generalbubbleplot
#'
#' This function is to make bubbleplot out of a data pair.  Notice,  genelist is a dataframe that contain a column that is gene name, where the column name is "gene",  if fa facet plot is wanted , an extra column named "tag" is required
#' @param pair.data the data pair that was prepred by datapair.mk
#' @param cpcol The column name for comparison, the default is name
#' @param toskip A vector of variables from col(eg,  a cell type), that I don't want to show in bubble plot
#' @param genelist A vector of genes that are interesting
#' @param hlen  The high color span length, default is 30
#' @param llen The low color span length, default is 25
#' @param midlen  The middle color span length, default is 0
#' @param limlen   The extreme high or low color span length
#' @param showcmcol   If true, print out a couple of lines of compared column
#' @param titlename   The figure title
#' @param angl  The angle of x axis that I want it to rotate
#' @param donormalscale  If true,  further normalize to zero for all genes
#' @param doreturn  If true,  reurn a list including bubble data and bubble plot
#' @param usefacet  If true,  plot with facet, an extra column in genelist is required
#' @return  return data pair that can be used for DE, bubble plot and others.
#' @export
#' @examples
#' Generalbubbleplot(ROCKvsnorock.non.paired,cpcol="name",genelist=c("TP53","TNFRSF1A","BAK1","CASP1"),donormalscale=F)

Generalbubbleplot<-function(pair.data=NULL,cpcol="celltype2",toskip=NULL,genelist=NULL,hlen=30,llen=25,midlen=0,limlen=30,showcmcol=F,titlename="",angl=0,donormalscale=F,doreturn=F,usefacet=F)
#									#this is the column to do compare; toskip is used if one or more levels in cpcol is not gonna used for compare in bubble
{
	require(RColorBrewer)
	#  A color pallet prepare
	Getpallet<-function(highlen,lowlen,wN,n,midcolor="white",killwhitelow=0,killwhitehigh=0)
	{
	myred<-colorRampPalette(brewer.pal(9,"Reds"))(highlen)[killwhitehigh:length(colorRampPalette(brewer.pal(9,"Reds"))(highlen))]
	myblue<-rev(colorRampPalette(brewer.pal(9,"Blues"))(lowlen))[1:(length(rev(colorRampPalette(brewer.pal(9,"Blues"))(lowlen)))-killwhitelow)]
	mypanel<-c(rep(myblue[1],n),myblue,rep(midcolor,wN),myred,rep(myred[length(myred)],n))
	return(mypanel)
	}

	if(is.null(pair.data))
	{
		print("parameter for this function is likebelow  :pair.data,cpcol='celltype2',toskip=NULL,genelist,hlen=30,llen=25,midlen=0,limlen=30,showcmcol=F")
		print("please enter the XXX.pair from DE analysis for 'pair.data' ")
	}
	if(is.null(genelist))
	{
		print("please enter a vector of genes you are interested for 'genelist'")
	}
	if(showcmcol==T)
	{
		print(head(pair.data$info))
		print(levels(pair.data$info[,ncol(pair.data$info)]))
		stop("Please choose from above and set showcmcol as F")
	}else{
	newlevel<-setdiff(levels(pair.data$info[,cpcol]),toskip)
	pair.data$info<-pair.data$info[which(!(pair.data$info[,cpcol] %in% toskip)),]
	pair.data$info[,cpcol]<-factor(pair.data$info[,cpcol],levels=newlevel)
	bubbledata<-subset(bubblePrep(pair.data,cpcol), gene %in% row.names(genelist))
	bubbledata$groupname<-factor(bubbledata$groupname,levels=newlevel)
		if(donormalscale)
		{
			scaled.result=c()
			n=1
				for(curgene in unique(bubbledata$gene))
				{
					scaled.result<-rbind(scaled.result,data.frame(gene=subset(result,gene==curgene)[,1],groupname=subset(result,gene==curgene)[,2],nonzeroratio=subset(result,gene==curgene) %>% .[,4] %>% scale, nTrance=subset(result,gene==curgene) %>% .[,4] %>% scale))
					print(n)
					n=n+1
				}
				bubbledata<-scaled.result
		}
		bubbledata$gene<-factor(bubbledata$gene,levels=row.names(genelist))
		if(usefacet)
		{
		bubbledata<-merge(bubbledata,genelist,by="gene")
		p<-ggplot(bubbledata)
		p<-p+aes(groupname,gene,color=log10(nTrance),size=nonzeroratio)+geom_point()+scale_color_gradientn(colors=Getpallet(hlen,llen,midlen,limlen))+ggtitle(titlename)
		p<-p+theme(axis.text.x = element_text(angle = angl, hjust = 1))+facet_grid("tag~.",scales="free",space="free")
		print(p)
	}else
		p<-ggplot(bubbledata)
		p<-p+aes(groupname,gene,color=log10(nTrance),size=nonzeroratio)+geom_point()+scale_color_gradientn(colors=Getpallet(hlen,llen,midlen,limlen))+ggtitle(titlename)
		p<-p+theme(axis.text.x = element_text(angle = angl, hjust = 1))
		print(p)
	}
	if(doreturn)
	{
		return(list(bubbledata,p))

	}
}


#' bubblePrep
#'
#' This function is the internal function to generate bubble data for plotting
#' @param pair.data the data pair that was prepred by datapair.mk
#' @param cpcol The column name for comparison, the default is name
#' @return  return bubble plot for plotting
#' @export



bubblePrep<-function(datapair,cpcol)
{
	Nonzeroratio.cal<-function(vector)
	{
		nonzero.n<-length(which(vector>0))
		nonzero.ratio<-nonzero.n/length(vector)
	}
		datapair$data[is.na(datapair$data)]<-0
		result<-data.frame()
		for (groupname in unique(datapair$info[,cpcol]))
		{
			tmpdata<-datapair$data[,row.names(subset(datapair$info,get(cpcol)==groupname))]
			tmpresult<-data.frame(gene=row.names(tmpdata),groupname=groupname,nonzeroratio=apply(tmpdata,1,Nonzeroratio.cal)/100*mean(apply(tmpdata,1,Nonzeroratio.cal)),nTrance=apply(tmpdata,1,mean)/(mean(apply(tmpdata,1,mean))*10),row.names=NULL)
			result<-rbind(result,tmpresult)
		}
	return(result)
}
