#' Prepareforpseudoregress.g
#'
#' This function is to perform the initial regression to prepare the trajectory study
#' @param object The seurat object for Repact study
#' @param PCrange The PC dimensions used for modeling
#' @param phenodic.use If there are extra data to be added, then this is a dataframe taht with one column named "Sample" for merge, also, the object@data.info has to have a column named "Sample". Otherwise, if the column for regression is already in object@data.info and phenodic.use=NULL
#' @param pheno  The column name used for modeling
#' '@param linear If true, perform a linear regression, else, perform a logistic regression
#' @return  The output include: PCvariance this summrize the percentage of variance chosen PC can explain; PCanfpheno this is a data frame including PC information, phenotype information as well as pseudo.index and residues after regression; $object.raw.withinfo: this is a datafrom of raw data with pseudo.index and residues; $model: this is the model. reg.plot.2d this is the 2d regression plot. model.para: this is the regression parameter.
#' @export
#' @examples
#' BMI.tjct.ob<-Prepareforpseudoregress.g(Beta.HSnegonly.ob,PCrange=1:10,phenodic.use=phenodic,pheno="BMI",linear=T)

Prepareforpseudoregress.g<-function(object=NULL,PCrange=1:10,phenodic.use=NULL,pheno,linear=T)
{   # if downsample==T, I'll do downsampling for donor3'
	require(ggplot2)
	require(RColorBrewer)
	require(gridExtra)

	if (is.null(object))
	{
		print("The parameter to enter is like : object=NULL,PCrange=1:40,phenodic.use=phenodic,pheno,linear=T,fam='binomial'")
		print("object(This is a Seutrat object),PCrange=1:40(decide which PCs to use for regression),phenodic.use=phenodic(This is a phenotype table loaded from outside),pheno(this is the column name that is used for regresion analysis),linear=T(if it is F, a non-linear relationship such as T2D/nonT2D that could be logistic using binomial),fam='binomial'(only used when linear=F)")
	}else
	{
		PCvariance<-summary(object@pca.obj[[1]])$importance[,40]
		KeyPCinformation<-data.frame(name=row.names(GetPCAcelldata_v2(object)),GetPCAcelldata_v2(object)[,c(PCrange,(ncol(GetPCAcelldata_v2(object))-1),(ncol(GetPCAcelldata_v2(object))-3))])
		if(!is.null(phenodic.use))
		{
		PCandPheno<-merge(KeyPCinformation,phenodic.use,by="Sample")
	}else
	{
		PCandPheno<-KeyPCinformation
	}
		PCnames<-paste("PCandPheno$PC",PCrange,sep="")
		form<-formula(paste("PCandPheno[,pheno]",paste(PCnames,collapse="+"),sep="~"))
		if(linear==T)
		{
			model<-lm(form)
			PCandPheno[,pheno]<-as.factor(PCandPheno[,pheno])
			p<-ggplot(PCandPheno)+aes_string("PC1","PC2",color=pheno)+geom_point()+scale_color_brewer(palette="YlOrRd")+geom_abline(slope=model$coefficients[3]/model$coefficients[2])
			PCandPheno<-cbind(PCandPheno,pseudo.index=model$fitted.values,residues=model$residuals)
			PCandPheno[,pheno]<-as.numeric(as.character(PCandPheno[,pheno]))
			model.para<-summary(model)
			model.para<-list(modelsummary=model.para)
		}else
		{
			PCandPheno[, pheno]<-factor(PCandPheno[, pheno])
			model<-glm(form,family="binomial")
			p<-ggplot(PCandPheno)+aes_string("PC1","PC2",color=pheno)+geom_point()+scale_color_brewer(palette="Set2")+geom_abline(slope=model$coefficients[3]/model$coefficients[2])
			PCandPheno<-cbind(PCandPheno,pseudo.index=model$linear.predictors,residues=model$residuals)
			model.para<-summary(model)
			lrt.pvalue<-1-pchisq(summary(model)$null.deviance-summary(model)$deviance,summary(model)$df.null-summary(model)$df.residual)
			pR2<-pR2(model)[4]
			model.para<-list(modelsummary=model.para,lrt.pvalue=lrt.pvalue,pR2=pR2)
		}

		row.names(PCandPheno)<-PCandPheno$name
		PCandPheno.names<-row.names(PCandPheno)
		object.raw<-t(as.matrix(object@raw.data))
		object.raw<-object.raw[PCandPheno.names,]
		object.raw<-cbind(object.raw,totalreads=apply(object.raw,1,sum))
		object.raw.withinfo<-Tomerge_v2(object.raw,PCandPheno[,c(ncol(PCandPheno)-1,ncol(PCandPheno))])
		return(list(PCvariance=PCvariance,PCanfpheno=PCandPheno,object.raw.withinfo=object.raw.withinfo,model=model,reg.plot.2d=p,model.para=model.para))
	}
}


#' Toplot3Dtjct
#'
#' This function is to make 3D example plot for regression analysis
#' @param trajectory.ob  This is the object generated from Prepareforpseudoregress.g
#' @param PCrange The PC dimensions used for 3D plotting,  usually pick top 3 significant pCs
#' @param pheno the column name for regression
#' @param linear If true, perform a linear regression, else, perform a logistic regression
#' '@param fam   For logistic regression, "binomial" is used for glm
#' '@param enlag    this is to adjust the length of regression line, default is 10
#' '@param decided   If TRUE, then only plot one based on angles by theta, and phi, otherwise, plot a series of figures in one big PDF.
#' '@param theta   A number indicating angle, I will work if decided
#' '@param phi    A another number indicating angle, I will work if decided
#' '@param singeplotname   a pdf file name if I have a decided single pdf
#' '@param multiplotname  a pdf file name if generating a series of figures
#' '@param titlename  This is the name in figure title
#' '@param theta   A number indicating angle, I will work if decided

#' @return  The output include: PCvariance this summrize the percentage of variance chosen PC can explain; PCanfpheno this is a data frame including PC information, phenotype information as well as pseudo.index and residues after regression; $object.raw.withinfo: this is a datafrom of raw data with pseudo.index and residues; $model: this is the model. reg.plot.2d this is the 2d regression plot. model.para: this is the regression parameter.
#' @export
#' @examples
#' Toplot3Dtjct(T2D.tjct.ob,PCrange=c(1,3,4),pheno="disease",linear=F,multiplotname="test.pdf",titlename="tittle")

Toplot3Dtjct<-function(trajectory.ob=tjct.ob,PCrange=1:3,pheno,linear=T, fam="binomial",enlag=10,decided=F,theta=60,phi=180,singeplotname=NULL,multiplotname,titlename="PC1-3 BMI trajectory regression\n")
{
	Getpallet2<-function(wN,topn,highlen,lowcolor="white")
	{
	myred<-colorRampPalette(brewer.pal(9,"Reds"))(highlen)
	#myblue<-rev(colorRampPalette(brewer.pal(9,"Blues"))(lowlen))
	mypanel<-c(rep(lowcolor,wN),myred,rep(myred[length(myred)],topn))
	return(mypanel)
	}
	PCandPheno<-trajectory.ob$PCanfpheno
	PCnames<-paste("PCandPheno$PC",PCrange,sep="")
	PCshornames<-paste("PC",PCrange,sep="")
	form<-formula(paste("PCandPheno[,pheno]",paste(PCnames,collapse="+"),sep="~"))
	if(linear==T)
		{
			model.eg<-lm(form)
			intercept<-model.eg$coefficients[1]
			a<-model.eg$coefficients[2]
			b<-model.eg$coefficients[3]
			c<-model.eg$coefficients[4]
			adjustrange=seq(0, 180, length.out = 13)
			if(decided)
			{
				pdf(singeplotname)
				title<-paste(titlename,"theta=",theta,"phi=",phi)
				scatter3D(PCandPheno[,PCshornames[1]],PCandPheno[,PCshornames[2]],PCandPheno[,PCshornames[3]],ticktype = "detailed",pch=20,theta=theta,phi=phi,colvar=as.integer(PCandPheno[,pheno]),bty="b2",cex=0.6,col=alpha.col(col=c(Getpallet2(0,0,11)[3:(length(Getpallet2(0,0,11))-2)]),0.4),main=title,xlab="PC1",ylab="PC2",zlab="PC3")+scatter3D(x=c(a*enlag,-a*enlag),y=c(b*enlag,-b*enlag),z=c(-c*enlag,c*enlag),type="l",ticktype = "detailed",add=T,colkey=F)
				dev.off()
			}else
			{
				pdf(multiplotname)
				for (i in adjustrange)
				{
					for(j in adjustrange)
					{
						title<-paste(titlename,"theta=",i,"phi=",j)
						scatter3D(PCandPheno[,PCshornames[1]],PCandPheno[,PCshornames[2]],PCandPheno[,PCshornames[3]],ticktype = "detailed",pch=20,theta=i,phi=j,colvar=as.integer(PCandPheno[,pheno]),bty="b2",cex=0.6,col=alpha.col(col=c(Getpallet2(0,0,11)[3:(length(Getpallet2(0,0,11))-2)]),0.4),main=title,xlab="PC1",ylab="PC2",zlab="PC3")+scatter3D(x=c(a*enlag,-a*enlag),y=c(b*enlag,-b*enlag),z=c(-c*enlag,c*enlag),type="l",ticktype = "detailed",add=T,colkey=F)
					}
				}
				dev.off()
			}
		}else
		{
			model.eg<-glm(form,family=fam)
			intercept<-model.eg$coefficients[1]
			a<-model.eg$coefficients[2]
			b<-model.eg$coefficients[3]
			c<-model.eg$coefficients[4]
			adjustrange=seq(0, 180, length.out = 13)
			if (decided)
			{
				pdf(singeplotname)
				title<-paste(titlename,"theta=",theta,"phi=",phi)
				scatter3D(PCandPheno[,PCshornames[1]],PCandPheno[,PCshornames[2]],PCandPheno[,PCshornames[3]],ticktype = "detailed",pch=20,theta=theta,phi=phi,colvar=as.integer(PCandPheno[,pheno]),bty="b2",cex=0.6,col=alpha.col(col=c("blue","red"),0.4),main=title,xlab="PC1",ylab="PC2",zlab="PC3",colkey = list(at = c(1, 2), side = 4, addlines = TRUE, length = 0.5, width = 0.5,labels = levels(PCandPheno[,pheno])))+scatter3D(x=c(a*enlag,-a*enlag),y=c(-b*enlag,b*enlag),z=c(c*enlag,-c*enlag),type="l",ticktype = "detailed",add=T,colkey=F)
				dev.off()
			}else
			{
				pdf(multiplotname)
				for (i in adjustrange)
				{
					for(j in adjustrange)
					{
						title<-paste(titlename,"theta=",i,"phi=",j)
						scatter3D(PCandPheno[,PCshornames[1]],PCandPheno[,PCshornames[2]],PCandPheno[,PCshornames[3]],ticktype = "detailed",pch=20,theta=i,phi=j,colvar=as.integer(PCandPheno[,pheno]),bty="b2",cex=0.6,col=alpha.col(col=c("blue","red"),0.4),main=title,xlab="PC1",ylab="PC2",zlab="PC3",colkey = list(at = c(1, 2), side = 4, addlines = TRUE, length = 0.5, width = 0.5,labels = levels(PCandPheno[,pheno])))+scatter3D(x=c(a*enlag,-a*enlag),y=c(-b*enlag,b*enlag),z=c(c*enlag,-c*enlag),type="l",ticktype = "detailed",add=T,colkey=F)
					}
				}
			dev.off()
			}
		}

}




#' Tjct.core.gen
#'
#' This function is to compute significant trajectory genes by linear regression
#' @param object The seurat object for Repact study
#' @param binnumber  Number of bins to divide the whole trajectory, default is 20.  This caqn vary upon total cell numbers available
#' @param qcut  q value cutoff to be used ti call a trajectory gene hits.
#' @return  A list of upregulated and downregulated trajectory genes.
#' @export
#' @examples
#' T2D.tjct.2nd.ob<-Tjct.core.gen(T2D.tjct.ob)

Tjct.core.gen<-function(object=NULL,binnumber=20,qcut=0.05)
{
	mydplyr<-function(df,by="tag",thefunction=mean)
	{
		newdf<-c()
		for (i in as.character(unique(df[,by])))
		{
		dftm<-df[which(as.character(df[,by])==i),]
		newdf<-rbind(newdf,apply(dftm,2,thefunction))
	 	}
	return(newdf)
	}
	DopseuBINregress<-function(rawdata,regressby="pseudo.index",start=1,end=ncol(rawdata)-4){
	require(qvalue)
	  pseudoregress.all<-c()
	  for (i in start:end)     #Loop gene by gene
	  {
	   gmodel1<-glm(rawdata[,i]~rawdata[,regressby],offset=log10(totalreads),data=rawdata,family=gaussian)
	   cur.slope1<-summary(gmodel1)$coefficients[,1][2]   # to get slope of pseudoBMIindex
	   cur.pvalues1<-summary(gmodel1)$coefficients[,4][2]  # to get p value for pseudoBMIindex
	   names(cur.slope1)<-paste("slope.",regressby,sep="")
	   names(cur.pvalues1)<-paste("p.",regressby,sep="")
	#   names(meanUMI)<-"meanUMI"
	   pseudoregress.all<-rbind(pseudoregress.all,c(cur.slope1,cur.pvalues1))
	  }
	  pseudoregress.all<-as.data.frame(pseudoregress.all)
	  row.names(pseudoregress.all)<-names(rawdata)[start:end]
	return(pseudoregress.all)
	}
	givebintag<-function(df,ordername="pseudoBMIindex",bin=10,genenumbercut=50){
	leftendpoint<-quantile(df[,ordername],genenumbercut/nrow(df))
	rightendpoint<-quantile(df[,ordername],1-genenumbercut/nrow(df))
	binunit<-(rightendpoint-leftendpoint)/(bin-2)
	cutoff1<-leftendpoint
	cutoff2<-cutoff1+binunit
	roadmark<-c(cutoff1,cutoff2)
	newdf<-c()
	newdf<-rbind(newdf,cbind(df[which(df[,ordername]<=leftendpoint),],tag=1))
	newdf<-rbind(newdf,cbind(df[which(df[,ordername]>rightendpoint),],tag=bin))
	for (i in 2:(bin-1)){
	newdf<-rbind(newdf,cbind(df[which(df[,ordername]>cutoff1 & df[,ordername]<=cutoff2),],tag=i))
	cutoff1<-cutoff2
	cutoff2<-cutoff2+binunit
	roadmark<-c(roadmark,cutoff2)
	}
	roadmark<-roadmark[-length(roadmark)]
	return(list(newdf,roadmark))
	}
	GetBirank.g<-function(df,cutoff=0.05,Tscut=100,rankby="slope"){  # rankby == "slope" or "qvalue"
	require(qvalue)
	tellquality<-function(vector){   # T2D  pos=3,  BMI pos=4
	if (vector["qvalue"]>cutoff){
	  return("NS")
	  }
	else if (vector["slope.pseudo.index"]>0){
	  return("UP")
	 }

	else if (vector["slope.pseudo.index"]<0){
	  return("DOWN")
	}
	}
	df2<-cbind(df,qvalue=qvalue(df[,"p.pseudo.index"])$qvalues)
	df3<-subset(df2,Transcripts>Tscut & qvalue<cutoff)
	df4<-cbind(df3,feature=apply(df3,1,tellquality))
	UP<-subset(df4,feature=="UP")
	DOWN<-subset(df4,feature=="DOWN")
	NS<-subset(df4,feature=="NS")
	if (rankby=="slope")
	{
		UP.rank<-data.frame(rank=rank(-abs(UP$slope.pseudo.index)),row.names=row.names(UP))
		DOWN.rank<-data.frame(rank=rank(-abs(DOWN$slope.pseudo.index)),row.names=row.names(DOWN))
	}else if (rankby=="qvalue")
	{
		UP.rank<-data.frame(rank=rank(-abs(UP$qvalue)),row.names=row.names(UP))
		DOWN.rank<-data.frame(rank=rank(-abs(DOWN$qvalue)),row.names=row.names(DOWN))
	}

	UP.rank<-UP.rank/max(UP.rank)
	DOWN.rank<-DOWN.rank/max(DOWN.rank)
	rank.dict<-rbind(UP.rank,DOWN.rank)
	df4<-Tomerge_v2(df4,rank.dict)
	UP<-subset(df4,feature=="UP")
	UP<-UP[order(UP$rank),]

	DOWN<-subset(df4,feature=="DOWN")
	DOWN<-DOWN[order(DOWN$rank),]

	return(list(UP=UP,DOWN=DOWN))

	}
#---------------
	geneUMI.dic<- data.frame(Transcripts=apply(object$object.raw.withinfo,2,sum)[-((ncol(object$object.raw.withinfo)-2):ncol(object$object.raw.withinfo))])
	# To add Bin tag very everey single cell   #  The result is list with two elements., First is the raw data with info;  Second is the position of making bin
	raw.bin<-givebintag(object$object.raw.withinfo,ordername="pseudo.index",genenumbercut=round(nrow(object$object.raw.withinfo)/binnumber),bin=binnumber)
	print("rawBin has been made")
	# To calculate mean number of every bin for every gene
	bin.data<-as.data.frame(mydplyr(raw.bin[[1]]))
	# Do linear regression for BInTjct.core(alpha.BMI.tjct.ob
	binregress<-DopseuBINregress(bin.data,regressby="pseudo.index")
	# To get a plotable result
	print("regression finished")
	BINlinear.result<-Tomerge(binregress,geneUMI.dic)
	BINlinear.result.summarized<-GetBirank.g(BINlinear.result,rankby="slope",cutoff=qcut)
	# to generate plots
	print("calculation finished")
	parameter<-paste("Take",binnumber,"bins", "the q value to identify sig genes is 1<=",qcut)
	print("in total identified genes are UP=")
	print(row.names(BINlinear.result.summarized$UP))
	print("in total identified genes are DOWN=")
	print(row.names(BINlinear.result.summarized$DOWN))
	return(list(geneUMI.dic=geneUMI.dic,raw.bin=raw.bin,bin.data=bin.data,binregress=binregress,BINlinear.result=BINlinear.result,BINlinear.result.summarized=BINlinear.result.summarized,parameter=parameter))
}



#' Tjct.core.plot
#'
#' This function is to generate major plots for Repact analysis
#' @param object The first trajectory generated by Prepareforpseudoregress.g
#' @param secondobj  The regressed and bined object by Tjct.core.gen
#' @param pheno  This is the column we are interested for RePACT gotta be the same with previous Prepareforpseudoregress.g
#' @param f1.name  The pdf name of violin plot e.g,"XX.10d.violin.pdf".  10 d means 10 PCs used.
#' @param f2.name  The pdf name of histgram e.g,"XX.his.pdf".
#' @param f3.name  The pdf name of the heatmap e.g,"XX.trj.heatmap.pdf".
#' @param f3.height  The pdf heigh for heatmap, change this when the number of genes shown is changed, default is 12
#' @param f3.tittle  The figure title for heatmap.e.g,"cell type:Changing genes on phenotype trajectory\ntop 0.06"
#' @param table1.name  The name of a csv file  for upregulated genes e.g,"XX.traj.up.genes-q0.05top0.06.csv"
#' @param table2.name  The name of a csv file  for upregulated genes e.g,"XX.traj.up.genes-q0.05top0.06.csv"
#' @param rankcut The percentile to be shown on the heatmap. default is 0.06
#' @param colorset The color pallete"Set1" "Set2" "Set3"
#' @return  A list of upregulated and downregulated trajectory genes.
#' @export
#' @examples
#' Tjct.core.plot(BMI.tjct.ob,BMI.tjct.2nd.ob,pheno="BMI",f1.name="BMI.tjct.10d.violin.pdf",f2.name="BMI.tjct.his.pdf",f3.name="BMI.tjct.trj.heatmap.pdf",f3.height=14,f3.tittle="cell type:Changing genes on phenotype trajectory\ntop6%",table1.name="BMI.tjct.traj.up.genes-q0.05Full.csv",table2.name="BMI.tjct.traj.dowb.genes-q0.05Full.csv",rankcut=0.05)

Tjct.core.plot<-function(object=NULL,secondobj=NULL,pheno=NULL,f1.name="XX.10d.violin.pdf",f2.name="XX.his.pdf",f3.name="XX.trj.heatmap.pdf",f3.height=12,f3.tittle="cell type:Changing genes on phenotype trajectory\ntop6%",table1.name="XX.traj.up.genes-q0.05top0.06.csv",table2.name="XX.traj.dowb.genes-q0.05top0.06.csv",rankcut=0.06,colorset="Set1",do.return=F)
{
	Do_heatmap<-function(bindata,df1,df2,rankname,cutoff=0.3,title,hardadd=NULL,last=T,insertinto=0,doreturn=F){
	require(ggplot2)
	#require(reshape2)
	require(gridExtra)
	DEgenelist<-c(row.names(df1)[as.numeric(as.character(df1[,rankname]))<cutoff],row.names(df2)[as.numeric(as.character(df2[,rankname]))<cutoff])
	if(last==T)
	{
	#DEgenelist<-c(DEgenelist,hardadd)
	DEgenelist.L<-DEgenelist[1:(length(DEgenelist)-insertinto)]
	DEgenelist.R<-setdiff(DEgenelist,DEgenelist.L)
	DEgenelist<-c(DEgenelist.L,hardadd,DEgenelist.R)
	}else
	{
	#DEgenelist<-c(hardadd,DEgenelist)
	DEgenelist.R<-DEgenelist[(insertinto+1):length(DEgenelist)]
	DEgenelist.L<-setdiff(DEgenelist,DEgenelist.R)
	DEgenelist<-c(DEgenelist.L,hardadd,DEgenelist.R)
	}
	bindata<-bindata[,c(DEgenelist,"tag")]
	bindata<-data.frame(apply(bindata[,-ncol(bindata)],2,normalize_01),bin=bindata[,ncol(bindata)])
	bindata.m<-reshape2::melt(bindata,id.vars="bin")
	p<-ggplot(bindata.m)+aes(bin,variable,fill=value)+geom_tile()+scale_fill_gradient2(low="white",high="red",mid="orange",midpoint=0.6)+labs(y="Variable genes")+ggtitle(title)
	if (doreturn)
	return(p)
	}
	normalize_01<-function(vector){
	normed<-(vector-min(vector))/(max(vector)-min(vector))
	return(normed)
	}
	if(is.null(pheno))
	{
		print(head(object$PCanfpheno))
		print("please enter a column name for histograme fill")
	}else if(is.null(object))
	{
		print("please enter parameter like below")
		print("object=tjct.ob,binnumber=20,qcut=0.05,f1.name='XX.10d.violin.pdf',f2.name='XX.his.pdf',f3.name='XX.trj.heatmap.pdf',f3.height=12,f3.tittle='cell type:Changing genes on phenotype trajectory\ntop6%',table1.name='XX.traj.up.genes-q0.05top0.06.csv',table2.name='XX.traj.dowb.genes-q0.05top0.06.csv")
	}else if (is.null(secondobj))
	{
		print ("please enter secoindary object, make sure it is consistant with the primary trajectory object")
	}else
	{
		bin.data<-secondobj$bin.data
		BINlinear.result.summarized<-secondobj$BINlinear.result.summarized
		raw.bin<-secondobj$raw.bin
		print("start to plot f1")
		pdf(f1.name)
		if(!is.factor(object$PCanfpheno[,pheno]))
		{
		p<-ggplot(object$PCanfpheno)+aes_string("Sample","pseudo.index",group="Sample",fill=pheno)+geom_violin()+coord_flip()+scale_fill_gradient(low="white",high="red")+geom_hline(yintercept=raw.bin[[2]],linetype=5,size=0.25)
		}else
		{
			p<-ggplot(object$PCanfpheno)+aes_string("Sample","pseudo.index",group="Sample",fill=pheno)+geom_violin()+coord_flip()+geom_hline(yintercept=raw.bin[[2]],linetype=5,size=0.25)+scale_fill_manual(values=c("blue","red"))
		}
		print(p)
		dev.off()
		print("start to plot f2")
		pdf(f2.name)
		p1<-ggplot(object$PCanfpheno)+aes(pseudo.index)+geom_histogram(fill="orange")+geom_vline(xintercept=raw.bin[[2]],linetype=5,size=0.25)
		p2<-ggplot(object$PCanfpheno)+aes(pseudo.index,residues,color=Sample)+geom_point(size=0.3)+geom_vline(xintercept=raw.bin[[2]],linetype=5,size=0.25)
		p3<-ggplot(object$PCanfpheno)+aes(pseudo.index,fill=Sample)+geom_histogram(position="fill")+geom_vline(xintercept=raw.bin[[2]],linetype=5,size=0.25)+scale_fill_brewer(palette=colorset)
		p4<-ggplot(object$PCanfpheno)+aes(pseudo.index,fill=Sample)+geom_histogram(position="stack")+geom_vline(xintercept=raw.bin[[2]],linetype=5,size=0.25)+scale_fill_brewer(palette=colorset)
		grid.arrange(p1,p2,ncol=1)
		grid.arrange(p3,p4,ncol=1)
		dev.off()
		print("start to plot f3")
		pdf(f3.name,height=f3.height)
		p5<-Do_heatmap(bin.data,df1=BINlinear.result.summarized$UP,df2=BINlinear.result.summarized$DOWN,rankname="rank",cutoff=rankcut,title=f3.tittle,doreturn=do.return)  # cutoff is the rank cutoff
		dev.off()
		print("start to print out tables")
		write.csv(BINlinear.result.summarized$UP,table1.name)
		write.csv(BINlinear.result.summarized$DOWN,table2.name)
		if(do.return){
			return(list(p,p1,p2,p3,p4,p5))
		}

	}
}
