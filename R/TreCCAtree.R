#' mydplyr.percentage
#'
#' This function is an internal function
#' @param df
#' @param by
#' @return
#' @export
#' @examples
#' Primaryclutering(dgelist,txcut=500,cluster.resos=first.reso,RAWinput=israw)

mydplyr.percentage<-function(df,by="stage")
{
  newdf<-c()
  for (i in as.character(unique(df[,by])))
  {
    dftm<-df[which(as.character(df[,by])==i),]
    current.sum<-sum(dftm[,2])
    newdf<-rbind(newdf,cbind(dftm,percent=apply(dftm,1,function(x){as.numeric(as.character(x[2]))/current.sum})))
  }
  return(newdf)
}


#' Primaryclutering
#'
#' This function is an advanced clustering function that when for example input dgeA, dgeB. This will do cluster individually for A, B and also for merged AB
#' @param rawdatapackage a list of several dge files. At least 2 dge. in the case of TreCCA tree. It is alwasy two dge, one from early time, one form latert
#' @param txcut The transcripts cutoff for clustering analysis
#' @param RAWinput If true, it is a raw dge
#' @param cluster.resos A vector of number, each is the cluster resolution for each dge.
#' @param allsamereso If true, then take the first number of cluster.resos for clustering of all dges
#' @param Together.reso The resolution of merged dge(all samples together)
#' @return  This will return a list of seurat result. It is used in function such as Tree.build.prepare
#' @export
#' @examples
#' Primaryclutering(dgelist,txcut=500,cluster.resos=first.reso,RAWinput=israw)

Primaryclutering<-function(rawdatapackage,txcut=500,RAWinput=T,cluster.resos=c(0.6,0.6),allsamereso=F,together.reso=0.6)
{
	Seurat.list<-list()
	for (i in 1:length(rawdatapackage))
	{
		if(allsamereso)
		{
			Seurat<-docluster(dgepreprocess(rawdatapackage[[i]],txcut,norowname=RAWinput),GetinformativeGene(dgepreprocess(rawdatapackage[[i]],txcut,norowname=RAWinput),500),names(rawdatapackage)[i],reso=cluster.resos[1])
		}else
		{
			Seurat<-docluster(dgepreprocess(rawdatapackage[[i]],txcut,norowname=RAWinput),GetinformativeGene(dgepreprocess(rawdatapackage[[i]],txcut,norowname=RAWinput),500),names(rawdatapackage)[i],reso=cluster.resos[i])
		}
		Seurat.list<-c(Seurat.list,list(Seurat))
	}
	names(Seurat.list)<-names(rawdatapackage)
	Sample.dict<-c()
	Seurat.raw.list<-list()
	for (i in 1:length(rawdatapackage))
	{
		donor_cell<-data.frame(Sample=Seurat.list[[i]]@data.info[,6],row.names=row.names(Seurat.list[[i]]@data.info))
		Sample.dict<-rbind(Sample.dict,donor_cell)
		Seurat.raw.list<-c(Seurat.raw.list,list(as.matrix(Seurat.list[[i]]@raw.data)))
	}
	names(Seurat.raw.list)<-names(Seurat.list)
	#step4  Make original clustering
	SeuratALL.origin<-docluster.multi(500,txcutoff=500,Seurat.raw.list,names(Seurat.list),reso=together.reso)
	return(list(Seurat.list=Seurat.list,Sample.dict=Sample.dict,Seurat.raw.list=Seurat.raw.list,SeuratALL.origin=SeuratALL.origin))
}

#' myddply.center
#'
#' This is an internal function used by dist.matrix.prep and dist.matrix.prep.v3
#' @param df
#' @param Col
#' @param is.originalcluster
#' @param combinecol
#' @return   This will return an object for dist.matrix.prep and dist.matrix.prep.v3
#' @export
#' @examples
#' Sample1.centers<-myddply.center(ob1.PCAS.clst,cluster.col1)

myddply.center<-function(df,Col,is.originalcluster=F,combinecol="Sample")  # only when you want to make a combination of two columns into a new clusrer name#  Note the myfun should be a function that return a vector
{
	df<-df[complete.cases(df),]
	Centers<-c()
	cluster<-unique(df[,Col])
	cellNumber<-c()
	size.matrix<-c()
	df<-df[complete.cases(df),]
	for (i in unique(df[,Col]))
	{
		if(is.originalcluster)
		{
			cluster<-c(cluster,paste(unique(df[,combinecol]),i,sep="_"))
		}
		df.sub<-df[df[,Col]==i,]
		curCenter<-colMeans(df.sub[,1:40])
		Centers<-rbind(Centers,c(colMeans(df.sub[,1:40])))
		cellNumber<-c(cellNumber,nrow(df.sub))
		cur.size.matrix<-as.matrix(summary(apply(df.sub[,1:40],1,function(x){dist(rbind(x,curCenter))})))
		colnames(cur.size.matrix)<-paste("",i,sep="")
		size.matrix<-cbind(size.matrix,cur.size.matrix)
	}
	Centers<-cbind(Centers,cellNumber=cellNumber)
	row.names(Centers)<-cluster
	if (nrow(Centers)==1)
	{
		centers.dist.mtx<-"No distance matrix, because only one clust is found"
	}else
	{
		centers.dist.mtx<-as.matrix(dist(Centers[,1:40]))
	}
	return(list(Centers=Centers,centers.dist.mtx=centers.dist.mtx,size.matrix=size.matrix))
}

#' myddply.center.v2
#'
#' This is another internal function used by dist.matrix.prep and dist.matrix.prep.v3
#' @param df.1
#' @param df.2
#' @param Col1
#' @param Col2
#' @return   This will return an object for dist.matrix.prep and dist.matrix.prep.v3
#' @export
#' @examples
#' Sample1.centers<-myddply.center(ob1.PCAS.clst,cluster.col1)

myddply.center.v2<-function(df.1,df.2,Col1,Col2)  #  Note the myfun should be a function that return a vector
{
	require(rdist)
	Centers.1<-c()
	cluster.1<-c()
	cellNumber.1<-c()
	size.matrix.1<-c()
	for (i in unique(df.1[,Col1]))
	{
		cluster.1<-c(cluster.1,paste("",i,sep=""))
		df.sub.1<-df.1[as.character(df.1[,Col1])==i,]
		curCenter.1<-colMeans(df.sub.1[,1:40])
		df.sub.1<-df.sub.1[complete.cases(df.sub.1),]
		Centers.1<-rbind(Centers.1,c(colMeans(df.sub.1[,1:40])))
		cellNumber.1<-c(cellNumber.1,nrow(df.sub.1))
		cur.size.matrix.1<-as.matrix(summary(apply(df.sub.1[,1:40],1,function(x){dist(rbind(x,curCenter.1))})))
		colnames(cur.size.matrix.1)<-paste("",i,sep="")
		size.matrix.1<-cbind(size.matrix.1,cur.size.matrix.1)
	}
	Centers.1<-cbind(Centers.1,cellNumber.1=cellNumber.1)
	row.names(Centers.1)<-cluster.1
	#``` for df.2
	Centers.2<-c()
	cluster.2<-c()
	cellNumber.2<-c()
	size.matrix.2<-c()
	for (i in unique(df.2[,Col2]))
	{
		cluster.2<-c(cluster.2,paste("",i,sep=""))
		df.sub.2<-df.2[df.2[,Col2]==i,]
		df.sub.2<-df.sub.2[complete.cases(df.sub.2),]
		curCenter.2<-colMeans(df.sub.2[,1:40])
		Centers.2<-rbind(Centers.2,c(colMeans(df.sub.2[,1:40])))
		cellNumber.2<-c(cellNumber.2,nrow(df.sub.2))
		cur.size.matrix.2<-as.matrix(summary(apply(df.sub.2[,1:40],1,function(x){dist(rbind(x,curCenter.2))})))
		colnames(cur.size.matrix.2)<-paste("",i,sep="")
		size.matrix.2<-cbind(size.matrix.2,cur.size.matrix.2)
	}
	Centers.2<-cbind(Centers.2,cellNumber.2=cellNumber.2)
	row.names(Centers.2)<-cluster.2[cluster.2!="NA"]
	center.1.PCs<-Centers.1[,1:40]
	center.2.PCs<-Centers.2[,1:40]
	if(!is.matrix(center.1.PCs))
	{
		center.1.PCs<-t(center.1.PCs)
	}
	if(!is.matrix(center.2.PCs))
	{
		center.2.PCs<-t(center.2.PCs)
	}
	if(is.null(row.names(center.1.PCs)))
	{
		row.names(center.1.PCs)<-cluster.1
	}
	if(is.null(row.names(center.1.PCs)))
	{
		row.names(center.2.PCs)<-cluster.2
	}
	centers.dist.mtx<-cdist(center.1.PCs,center.2.PCs)
	row.names(centers.dist.mtx)<-cluster.1
	colnames(centers.dist.mtx)<-cluster.2[cluster.2!="NA"]
	df1_df2.relation<-paste(row.names(centers.dist.mtx),colnames(centers.dist.mtx)[apply(centers.dist.mtx,1,function(x){which(x==min(x))})],sep="--")
	df2_df1.relation<-paste(colnames(centers.dist.mtx),row.names(centers.dist.mtx)[apply(centers.dist.mtx,2,function(x){which(x==min(x))})],sep="--")
	##Creat the segData which is critical for plotting the connection plot seg.Data.main  is a seg data seeking from later stage toward earlier stage which is a major considerarion
	seg.Data.main<-c()
	if(dim(centers.dist.mtx)[1]==1)
	{
		for (older.cluster.names in colnames(centers.dist.mtx))
		{
			average.dist<-mean(centers.dist.mtx)
			if(min(centers.dist.mtx[,older.cluster.names])<=average.dist)
			{
				seg<-c(older.cluster.names,row.names(centers.dist.mtx))
				seg.Data.main<-rbind(seg.Data.main,seg)
			}
		}
	}else
	{
		seg.Data.main<-c()
		for (older.cluster.names in colnames(centers.dist.mtx))
		{
			average.dist<-mean(centers.dist.mtx)
			if(min(centers.dist.mtx[,older.cluster.names])<=average.dist)
			{
				seg<-c(older.cluster.names,names(which(centers.dist.mtx[,older.cluster.names]==min(centers.dist.mtx[,older.cluster.names]))))
				seg.Data.main<-rbind(seg.Data.main,seg)
			}
		}
	}
	row.names(seg.Data.main)<-NULL
	seg.Data.main<-as.data.frame(seg.Data.main)
	names(seg.Data.main)<-c("SeekFrom.cluster.names","SeekToward.cluster.names")
	seg.Data.main<-cbind(seg.Data.main,SeekFrom.stage.names=unlist(lapply(strsplit(as.character(seg.Data.main$SeekFrom.cluster.names),"_"),function(x){x[1]})))
	seg.Data.main<-cbind(seg.Data.main,SeekToward.stage.names=unlist(lapply(strsplit(as.character(seg.Data.main$SeekToward.cluster.names),"_"),function(x){x[1]})))
	##Creat the segData which is critical for plotting the connection plot seg.Data.main  is a seg data seeking from earlier stage toward later stage which is a alternative considerarion
	seg.Data.alt<-c()
	if(dim(centers.dist.mtx)[2]==1)
	{
		for (earlier.cluster.names in row.names(centers.dist.mtx))
		{
			average.dist<-mean(centers.dist.mtx)
			if(min(centers.dist.mtx[earlier.cluster.names,])<=average.dist)
			{
				seg<-c(earlier.cluster.names,colnames(centers.dist.mtx))
				seg.Data.alt<-rbind(seg.Data.alt,seg)
			}
		}
	}else
	{
		seg.Data.alt<-c()
		for (earlier.cluster.names in row.names(centers.dist.mtx))
		{
			average.dist<-mean(centers.dist.mtx)
			if(min(centers.dist.mtx[earlier.cluster.names,])<=average.dist)
			{
				seg<-c(earlier.cluster.names,names(which(centers.dist.mtx[earlier.cluster.names,]==min(centers.dist.mtx[earlier.cluster.names,]))))
				seg.Data.alt<-rbind(seg.Data.alt,seg)
			}
		}
	}
	row.names(seg.Data.alt)<-NULL
	seg.Data.alt<-as.data.frame(seg.Data.alt)
	names(seg.Data.alt)<-c("SeekFrom.cluster.names","SeekToward.cluster.names")
	seg.Data.alt<-cbind(seg.Data.alt,SeekFrom.stage.names=unlist(lapply(strsplit(as.character(seg.Data.alt$SeekFrom.cluster.names),"_"),function(x){x[1]})))
	seg.Data.alt<-cbind(seg.Data.alt,SeekToward.stage.names=unlist(lapply(strsplit(as.character(seg.Data.alt$SeekToward.cluster.names),"_"),function(x){x[1]})))
	return(list(Centers.1=Centers.1,Centers.2=Centers.2,centers.dist.mtx=centers.dist.mtx,size.matrix.1=size.matrix.1,size.matrix.2=size.matrix.2,df1_df2.relation=df1_df2.relation,df2_df1.relation=df2_df1.relation,seg.Data.main=seg.Data.main,seg.Data.alt=seg.Data.alt))
}

#' Tree.build.prepare
#'
#' This function is the first function to start calculate the subtree from dgeA, dgeB. This  basically do the clustering and fullplot at low resolution
#' @param dge1 the dge file (usually from earlier timepoint), either raw or processed is acceptable
#' @param dge2 the dge file (usually from later timepoint), either raw or processed is acceptable
#' @param name1 The name for the dge1 (usually from earlier timepoint)
#' @param name2 The name for the dge2 (usually from later timepoint)
#' @param first.reso  The cluster resolution for dge1 and dge2 respectively. A vector of two number , for example c(0.03,0.03)
#' @param treemake.dir  A path. The workspace where the whole tree making project is on
#' @return  A list of stuff including seurat clutsering result; fullplots; the resolution used and  the path
#' @export
#' @examples
#'  S4_S5ROCK.tree.prep<-Tree.build.prepare(dge1=s4.B.dge,dge2=S5.ROCKII.dge,name1="S4.B",name2="S5.rock",first.reso=c(0.03,0.03),treemake.dir="/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/ROCKII/Treemake/")

Tree.build.prepare<-function(dge1,dge2,name1,name2,first.reso=c(0.03,0.03),treemake.dir="./")  # S4.B  apply resolusion 0.015
{
	setwd(treemake.dir)
	if(!any(grepl(paste(name1,name2,sep="_"),list.files()))){
		dir.create(paste(name1,name2,sep="_"))
	}
	setwd(paste(c("./",name1,"_",name2),collapse=""))
	path<-getwd()
	first.resos.used<-paste("first.reso is",paste(first.reso,collapse="_"))
	print(first.resos.used)
	#	dir.create(paste(name1,name2,sep="_"))
	#	setwd(paste(c("./",name1,"_",name2),collapse=""))
	#sink(paste(c(name1,"_",name2,".log"),collapse=""))
	dgelist<-list(dge1,dge2)
	names(dgelist)<-c(name1,name2)
	print(paste(c("Start to do clustering between",name1,name2),collapse=" "))
	## Primary clustering
	#
	israw<-colnames(dge1)[1]=="GENE"
	primary.total<-Primaryclutering(dgelist,txcut=500,cluster.resos=first.reso,RAWinput=israw)
	primary.total$Seurat.list[[1]]@data.info<-cbind(primary.total$Seurat.list[[1]]@data.info,Sample.2nd=paste(primary.total$Seurat.list[[1]]@data.info$Sample,primary.total$Seurat.list[[1]]@data.info[,paste("res.",first.reso[1],sep="")],sep="_"))
	primary.total$Seurat.list[[2]]@data.info<-cbind(primary.total$Seurat.list[[2]]@data.info,Sample.2nd=paste(primary.total$Seurat.list[[2]]@data.info$Sample,primary.total$Seurat.list[[2]]@data.info[,paste("res.",first.reso[2],sep="")],sep="_"))
	print("Total clustering with low resolution is completed")
	print(paste(c("with low resolution ",name1," can be clustered into ",length(unique(primary.total$Seurat.list[[1]]@data.info[,paste("res.",first.reso[1],sep="")]))," ; ",name2," can be clustered into ",length(unique(primary.total$Seurat.list[[2]]@data.info[,paste("res.",first.reso[2],sep="")]))),collapse=""))
	##
	#saveRDS(primary.total,file=paste(c(name1,name2,"totalprimary.RDS"),collapse="_"))
	## plot full plot for each two stage
	dir.create(paste(c("LowresCluster",name1,name2),collapse="."))
	setwd(paste(c("LowresCluster",name1,name2),collapse="."))
	print("Start to plot clustering into pdf...")
	total.fullplot.1<-Fullplot_v2(primary.total$Seurat.list[[1]],paste(c(name1,".total.",first.reso[1],".pdf"),collapse=""),signiture=NULL,resolusion=paste("res.",first.reso[1],sep=""),doreturn=T)
	print(paste(c("The yonger: ",name1," Fullplot has completed, the file name is ",paste(c(name1,".total.",first.reso[1],".pdf"),collapse="")),collapse=""))
	total.fullplot.2<-Fullplot_v2(primary.total$Seurat.list[[2]],paste(c(name2,".total.",first.reso[2],".pdf"),collapse=""),signiture=NULL,resolusion=paste("res.",first.reso[2],sep=""),doreturn=T)
	print(paste(c("The older: ",name2," Fullplot has completed, the file name is ",paste(c(name2,".total.",first.reso[2],".pdf"),collapse="")),collapse=""))
	return(list(primary.total=primary.total,total.fullplot.1=total.fullplot.1,total.fullplot.2=total.fullplot.2,path=path,first.resos.used=first.reso))
}


#' Prep.res.adjust
#'
#' This function is to adjust the clustering resolution generated from Tree.build.prepare
#' @param prep the output from Tree.build.prepare
#' @param first.reso.adjust    This is the ADJUSTED CLUSTER RESOLUTION.The format is the same as First.reso in Tree.build.prepare
#' @return  This will return the same object as that generated from Tree.build.prepare
#' @export
#' @examples
#'  S4_S5ROCK.tree.prep<-Prep.res.adjust(S4_S5ROCK.tree.prep,first.reso.adjust=c(0.015,0.015))

Prep.res.adjust<-function(prep,first.reso.adjust=c(0.03,0.03))  # S4.B  apply resolusion 0.015
{
	setwd(prep$path)
  path=prep$path
	name1<-names(prep$primary.total$Seurat.list)[1]
	name2<-names(prep$primary.total$Seurat.list)[2]
	setwd(paste(c("LowresCluster",name1,name2),collapse="."))
	prep$primary.total$Seurat.list[[1]]<- FindClusters(prep$primary.total$Seurat.list[[1]], pc.use = 1:10, resolution = first.reso.adjust[1], print.output = 0, save.SNN = T)
	prep$primary.total$Seurat.list[[1]]@data.info<-prep$primary.total$Seurat.list[[1]]@data.info[,c("nGene","nUMI","orig.ident","percent.mito",paste("res",first.reso.adjust[1],sep="."),"Sample")] %>% cbind(.,Sample.2nd=paste(.$Sample,.[,paste("res",first.reso.adjust[1],sep=".")],sep="_"))

	prep$primary.total$Seurat.list[[2]]<- FindClusters(prep$primary.total$Seurat.list[[2]], pc.use = 1:10, resolution = first.reso.adjust[2], print.output = 0, save.SNN = T)
	prep$primary.total$Seurat.list[[2]]@data.info<-prep$primary.total$Seurat.list[[2]]@data.info[,c("nGene","nUMI","orig.ident","percent.mito",paste("res",first.reso.adjust[2],sep="."),"Sample")] %>% cbind(.,Sample.2nd=paste(.$Sample,.[,paste("res",first.reso.adjust[2],sep=".")],sep="_"))

	total.fullplot.1<-Fullplot_v2(prep$primary.total$Seurat.list[[1]],paste(c(name1,".total.",first.reso.adjust[1],".pdf"),collapse=""),signiture=NULL,resolusion=paste("res.",first.reso.adjust[1],sep=""),doreturn=T)
	print(paste(c("The yonger: ",name1," Fullplot has completed, the file name is ",paste(c(name1,".total.",first.reso.adjust[1],".pdf"),collapse="")),collapse=""))

	total.fullplot.2<-Fullplot_v2(prep$primary.total$Seurat.list[[2]],paste(c(name2,".total.",first.reso.adjust[2],".pdf"),collapse=""),signiture=NULL,resolusion=paste("res.",first.reso.adjust[2],sep=""),doreturn=T)
	print(paste(c("The older: ",name2," Fullplot has completed, the file name is ",paste(c(name2,".total.",first.reso.adjust[2],".pdf"),collapse="")),collapse=""))
	return(list(primary.total=prep$primary.total,total.fullplot.1=total.fullplot.1,total.fullplot.2=total.fullplot.2,path=path,first.resos.used=first.reso.adjust))
}

#' dist.matrix.prep
#'
#' This is an important internal function to calculate distance matrix, it is used in Tree.build.1 and Tree.build.2nd.clustering
#' @param binary.primary This is the primary clustering result that can be retrieved from Tree.build.prepare object. (Because Tree.build.prepare mainly do the primary clustering.)The input here is tree.prep.ob$primary.total
#' @param datainfo.col1   default is c("res.0.6","nUMI","Sample")
#' @param cluster.col1 default is "res.0.6"
#' @param datainfo.col2 default is c("res.0.6","nUMI","Sample")
#' @param cluster.col2 default is "res.0.6"
#' @param res1  default is "res.0.06
#' @param res2 default is "res.0.06"
#' @param datainfo.col3 c("res.0.6","nUMI","Sample")
#' @return   This will return an object used for Tree.build.1
#' @export
#' @examples
#'  	total.dist.matrxes<-dist.matrix.prep(primary.total,datainfo.col1=c(paste("res.",first.reso[1],sep=""),"nUMI","Sample"),cluster.col1=paste("res.",first.reso[1],sep=""),datainfo.col2=c(paste("res.",first.reso[2],sep=""),"nUMI","Sample"),cluster.col2=paste("res.",first.reso[2],sep=""),res1=paste("res.",first.reso[1],sep=""),res2=paste("res.",first.reso[2],sep=""))

dist.matrix.prep<-function(binary.primary,datainfo.col1=c("res.0.6","nUMI","Sample"),cluster.col1="res.0.6",datainfo.col2=c("res.0.6","nUMI","Sample"),cluster.col2="res.0.6",res1="res.0.06",res2="res.0.06",datainfo.col3=c("res.0.6","nUMI","Sample"))
{
	print(paste("sample1 is:", names(binary.primary$Seurat.list)[1]))
	print(paste("sample2 is:", names(binary.primary$Seurat.list)[2]))
	samplename1<-names(binary.primary$Seurat.list)[1]
	samplename2<-names(binary.primary$Seurat.list)[2]
	##part1 to prepare the pca and information for each of the two sample
	####prepare for sample1
	ob1.PCAS.clst<-Tomerge_v2(binary.primary$Seurat.list[[1]]@pca.rot,binary.primary$Seurat.list[[1]]@data.info[,datainfo.col1])
	ob1.info<-binary.primary$Seurat.list[[1]]@data.info[,datainfo.col1]
	row.names(ob1.info)<-paste(samplename1,substr(row.names(ob1.info),1,12),sep="_")
	ob1.info[,cluster.col1]<-paste(samplename1,ob1.info[,res1],sep="_")
	####prepare for sample2
	ob2.PCAS.clst<-Tomerge_v2(binary.primary$Seurat.list[[2]]@pca.rot,binary.primary$Seurat.list[[2]]@data.info[,datainfo.col2])
	ob2.info<-binary.primary$Seurat.list[[2]]@data.info[,datainfo.col2]
	ob2.info<-ob2.info[!duplicated(substr(row.names(ob2.info),1,12)),]
	row.names(ob2.info)<-paste(samplename2,substr(row.names(ob2.info),1,12),sep="_")
	ob2.info[,cluster.col2]<-paste(samplename2,ob2.info[,res2],sep="_")
	####prepare for sample1--sample2
	obALL.PCAS.clst<-Tomerge_v2(binary.primary$SeuratALL.origin@pca.rot,binary.primary$SeuratALL.origin@data.info[,datainfo.col3])
	obALL.PCAS.clst.SP1<-subset(obALL.PCAS.clst,Sample==samplename1)
	obALL.PCAS.clst.SP2<-subset(obALL.PCAS.clst,Sample==samplename2)
	obALL.PCAS.clst.SP1<-obALL.PCAS.clst.SP1[which(!duplicated(paste(samplename1,substr(row.names(obALL.PCAS.clst.SP1),1,12),sep="_"))),]
	row.names(obALL.PCAS.clst.SP1)<-paste(samplename1,substr(row.names(obALL.PCAS.clst.SP1),1,12),sep="_")
	obALL.PCAS.clst.SP2<-obALL.PCAS.clst.SP2[which(!duplicated(paste(samplename2,substr(row.names(obALL.PCAS.clst.SP2),1,12),sep="_"))),]
	row.names(obALL.PCAS.clst.SP2)<-paste(samplename2,substr(row.names(obALL.PCAS.clst.SP2),1,12),sep="_")
	obALL.PCAS.clst.SP1<-Tomerge_v2(obALL.PCAS.clst.SP1[,1:40],ob1.info[,1:2])
	obALL.PCAS.clst.SP2<-Tomerge_v2(obALL.PCAS.clst.SP2[,1:40],ob2.info[,1:2])
	##part2 to calculate the center and distance matrix
	Sample1.centers<-myddply.center(ob1.PCAS.clst,cluster.col1)
	Sample2.centers<-myddply.center(ob2.PCAS.clst,cluster.col2)
	Sample1cross2.centers<-myddply.center.v2(obALL.PCAS.clst.SP1,obALL.PCAS.clst.SP2,cluster.col1,cluster.col2)
	return(list(Sample1.centers=Sample1.centers,Sample2.centers=Sample2.centers,Sample1cross2.centers=Sample1cross2.centers))
}


#' dist.matrix.prep.v3
#'
#' This is an important internal function to calculate distance matrix, it is used particularly in Tree.build.2nd.treemaking
#' @param binary.primary This is the primary clustering result that can be retrieved from Tree.build.prepare object. (Because Tree.build.prepare mainly do the primary clustering.)The input here is tree.prep.ob$primary.total
#' @param datainfo.col1   default is c("res.0.6","nUMI","Sample")
#' @param cluster.col1 default is "res.0.6"
#' @param datainfo.col2 default is c("res.0.6","nUMI","Sample")
#' @param datainfo.col.dog   c("res.0.6","nUMI","Sample")
#' @param cluster.col2 default is "res.0.6"
#' @param cluster.col.dog  "res.0.6"
#' @param withsingledog Default is F
#' @param dog.seurat.ob =NULL
#' @param neighborwithdog.seurat.ob =NULL
#' @return   This will return an object used for Tree.build.2nd.treemaking
#' @export
#' @examples
#'  	total.dist.matrxes<-dist.matrix.prep(primary.total,datainfo.col1=c(paste("res.",first.reso[1],sep=""),"nUMI","Sample"),cluster.col1=paste("res.",first.reso[1],sep=""),datainfo.col2=c(paste("res.",first.reso[2],sep=""),"nUMI","Sample"),cluster.col2=paste("res.",first.reso[2],sep=""),res1=paste("res.",first.reso[1],sep=""),res2=paste("res.",first.reso[2],sep=""))

dist.matrix.prep.v3<-function(binary.primary,datainfo.col1=c("res.0.6","nUMI","Sample"),datainfo.col2=c("res.0.6","nUMI","Sample"),datainfo.col.dog=c("res.0.6","nUMI","Sample"),cluster.col1="res.0.6",cluster.col2="res.0.6",cluster.col.dog="res.0.6",withsingledog=F,dog.seurat.ob=NULL,neighborwithdog.seurat.ob=NULL)
{
	require(ggdendro)
	#require(reshape2)
	print(paste("sample1 is:", names(binary.primary$Seurat.list)[1]))
	print(paste("sample2 is:", names(binary.primary$Seurat.list)[2]))
	samplename1<-names(binary.primary$Seurat.list)[1]
	samplename2<-names(binary.primary$Seurat.list)[2]
	##part1 to prepare the pca and information for each of the two sample
	####prepare for sample1
	ob1.PCAS.clst<-Tomerge_v2(binary.primary$Seurat.list[[1]]@pca.rot,binary.primary$Seurat.list[[1]]@data.info[,datainfo.col1])
	ob1.info<-binary.primary$Seurat.list[[1]]@data.info[,datainfo.col1]
	ob1.info<-cbind(ob1.info,Sample.cluster=paste(as.character(ob1.info$Sample),ob1.info[,cluster.col1],sep="_"))
	####prepare for sample2
	ob2.PCAS.clst<-Tomerge_v2(binary.primary$Seurat.list[[2]]@pca.rot,binary.primary$Seurat.list[[2]]@data.info[,datainfo.col2])
	ob2.info<-binary.primary$Seurat.list[[2]]@data.info[,datainfo.col2]
	row.names(ob2.info)<-paste(row.names(ob2.info),"b",sep="_")
	ob2.info<-cbind(ob2.info,Sample.cluster=paste(as.character(ob2.info$Sample),ob2.info[,cluster.col2],sep="_"))
	####prepare for singledog if it is there
	if (withsingledog)
	{
		print(paste("There is a single dog---",unique(dog.seurat.ob@data.info$Sample)))
		ob.dog.PCAS.clst<-Tomerge_v2(dog.seurat.ob@pca.rot,dog.seurat.ob@data.info[,datainfo.col.dog])
		ob.dog.info<-dog.seurat.ob@data.info
		row.names(ob.dog.info)<-paste(row.names(ob.dog.info),"c",sep="_")
		ob.dog.info<-cbind(ob.dog.info,Sample.cluster=paste(as.character(ob.dog.info$Sample),ob.dog.info[,cluster.col.dog],sep="_"))
		### put together with dog
		ob.1.2.dog.info<-rbind(ob1.info[,c("nUMI","Sample.cluster")],ob2.info[,c("nUMI","Sample.cluster")],ob.dog.info[,c("nUMI","Sample.cluster")])
		####prepare for sample1--sample2 plus singledog
		obALL.dog.PCAS.clst<-Tomerge_v2(neighborwithdog.seurat.ob@pca.rot,ob.1.2.dog.info)
		SampleALL.centers.dog<-myddply.center(obALL.dog.PCAS.clst,"Sample.cluster")
		hc.dog<-hclust(as.dist(SampleALL.centers.dog$centers.dist.mtx))
		p.dendro.dog<-ggdendrogram(hc.dog)
		label.order.dog<-hc.dog$labels[hc.dog$order]
		#Adjust the order of heatmap
		matrix.dog.m<-reshape2::melt(SampleALL.centers.dog$centers.dist.mtx)
		matrix.dog.m$Var1<-factor(matrix.dog.m$Var1,levels=label.order.dog)
		matrix.dog.m$Var2<-factor(matrix.dog.m$Var2,levels=label.order.dog)
		p.heat.dog<-ggplot(matrix.dog.m)+aes(Var1,Var2,fill=value)+geom_tile()+scale_fill_gradient(low="red",high="white")+geom_text(aes(label=format(value,digits=2)))+theme(axis.text.x=element_text(angle=45,hjust=1))
	}
	#### put together for that without singledog
	ob.1.2.info<-rbind(ob1.info[,c("nUMI","Sample.cluster")],ob2.info[,c("nUMI","Sample.cluster")])
	####prepare for sample1--sample2
	obALL.PCAS.clst<-Tomerge_v2(binary.primary$SeuratALL.origin@pca.rot,ob.1.2.info)
	##part2 to calculate the center and distance matrix
	#obALL.PCAS.clst<-cbind(obALL.PCAS.clst,sample_cluster=paste(as.character(obALL.PCAS.clst$Sample),obALL.PCAS.clst[,"Sample.cluster"],sep="_"))
	SampleALL.centers<-myddply.center(obALL.PCAS.clst,"Sample.cluster")
	# Do the dendro
	hc<-hclust(as.dist(SampleALL.centers$centers.dist.mtx))
	p.dendro<-ggdendrogram(hc)
	label.order<-hc$labels[hc$order]
	#Adjust the order of heatmap
	matrix.m<-reshape2::melt(SampleALL.centers$centers.dist.mtx)
	matrix.m$Var1<-factor(matrix.m$Var1,levels=label.order)
	matrix.m$Var2<-factor(matrix.m$Var2,levels=label.order)
	p.heat<-ggplot(matrix.m)+aes(Var1,Var2,fill=value)+geom_tile()+scale_fill_gradient(low="red",high="white")+geom_text(aes(label=format(value,digits=2)))+theme(axis.text.x=element_text(angle=45,hjust=1))
	#
	if(withsingledog)
	{
		part1.withoutdog<-list(p.dendro=p.dendro,p.heat=p.heat,matrix=SampleALL.centers$centers.dist.mtx)
		part2.withdog=list(p.dendro.dog=p.dendro.dog,p.heat.dog=p.heat.dog,matrix=SampleALL.centers.dog$centers.dist.mtx)
		return(list(part1.withoutdog=part1.withoutdog,part2.withdog=part2.withdog))
	}else
	{
		part1.withoutdog=list(p.dendro=p.dendro,p.heat=p.heat,matrix=SampleALL.centers$centers.dist.mtx)
		return(list(part1.withoutdog=part1.withoutdog))
	}
}




#' Tree.build.1
#'
#' This is the subtree relationship builder for low resolution
#' @param tree.prep.ob This is an object generated from Tree.build.prepare or Prep.res.adjust
#' @return   This will return an object used for Tree.build.1
#' @export
#' @examples
#'  	S4_S5ROCK.tree.1.ob<-Tree.build.1(S4_S5ROCK.tree.prep)

Tree.build.1<-function(tree.prep.ob)
{
  setwd(tree.prep.ob$path)
  name1<-names(tree.prep.ob$primary.total$Seurat.list)[1]
  name2<-names(tree.prep.ob$primary.total$Seurat.list)[2]
  setwd(paste(c("LowresCluster",name1,name2),collapse="."))
	primary.total<-tree.prep.ob$primary.total
	first.reso<-tree.prep.ob$first.resos.used
	name1<-names(tree.prep.ob$primary.total$Seurat.list)[1]
	name2<-names(tree.prep.ob$primary.total$Seurat.list)[2]
	total.dist.matrxes<-dist.matrix.prep(primary.total,datainfo.col1=c(paste("res.",first.reso[1],sep=""),"nUMI","Sample"),cluster.col1=paste("res.",first.reso[1],sep=""),datainfo.col2=c(paste("res.",first.reso[2],sep=""),"nUMI","Sample"),cluster.col2=paste("res.",first.reso[2],sep=""),res1=paste("res.",first.reso[1],sep=""),res2=paste("res.",first.reso[2],sep=""))   # This gave back distance matrix as well as segData that is critical for relationship drawing.  and importantly, the segData determines what will be further explored.  I reccomend manual inspection although integrated should be able to well tell the relationship
	dendro.heap.p<-dist.matrix.prep.v3(primary.total,datainfo.col1=c(paste("res.",first.reso[1],sep=""),"nUMI","Sample"),datainfo.col2=c(paste("res.",first.reso[2],sep=""),"nUMI","Sample"),cluster.col1=paste("res.",first.reso[1],sep=""),cluster.col2=paste("res.",first.reso[2],sep=""))
	##
	#  to creat internal relationship
	unpaired.cluster.names<-setdiff(colnames(dendro.heap.p$part1.withoutdog$matrix)[grepl(name1,colnames(dendro.heap.p$part1.withoutdog$matrix))],total.dist.matrxes$Sample1cross2.centers$seg.Data.main$SeekToward.cluster.names)
	total.dist.matrxes$Sample1cross2.centers$seg.Data.main[,1]<-factor(total.dist.matrxes$Sample1cross2.centers$seg.Data.main[,1],levels=colnames(dendro.heap.p$part1.withoutdog$matrix))
	total.dist.matrxes$Sample1cross2.centers$seg.Data.main[,2]<-factor(total.dist.matrxes$Sample1cross2.centers$seg.Data.main[,2],levels=colnames(dendro.heap.p$part1.withoutdog$matrix))
	total.dist.matrxes$Sample1cross2.centers$seg.Data.main.internal<-c()
	total.dist.matrxes<-Refine.relationship(total.dist.matrxes,dendro.heap.p$part1.withoutdog$matrix)
	total.relation.plots<-tree.plot.v2(dendrodata=dendro.heap.p,matrixes=total.dist.matrxes,textsize=2.5,youngername=name1,maturername=name2,nocellnumber=T)  #  This will give back a relationship network (low resolution)
	#dist.cutoff<-mean(dendro.heap.p$part1.withoutdog$matrix)   #  Take the whole matrix, including 0
	dist.cutoff<-mean(dendro.heap.p$part1.withoutdog$matrix[dendro.heap.p$part1.withoutdog$matrix>0])   #  Calculate a distance cutoff, out of which will be considered no connection
	for (unclust in unpaired.cluster.names)
	{
		if(min(dendro.heap.p$part1.withoutdog$matrix[unclust,][dendro.heap.p$part1.withoutdog$matrix[unclust,]>0])<dist.cutoff)
		{
			internaltarget.current<-names(dendro.heap.p$part1.withoutdog$matrix[unclust,][dendro.heap.p$part1.withoutdog$matrix[unclust,]>0])[dendro.heap.p$part1.withoutdog$matrix[unclust,][dendro.heap.p$part1.withoutdog$matrix[unclust,]>0]==min(dendro.heap.p$part1.withoutdog$matrix[unclust,][dendro.heap.p$part1.withoutdog$matrix[unclust,]>0])]
		}else
		{
			internaltarget.current<-NA
		}
		total.dist.matrxes$Sample1cross2.centers$seg.Data.main.internal<-rbind(total.dist.matrxes$Sample1cross2.centers$seg.Data.main.internal,c(SeekFrom.cluster.names=internaltarget.current,SeekToward.cluster.names=unclust,SeekFrom.stage.names=levels(total.dist.matrxes$Sample1cross2.centers$seg.Data.main[,3]),SeekToward.stage.names=levels(total.dist.matrxes$Sample1cross2.centers$seg.Data.main[,4])))
	}
	##
	pdf(paste(c(name1,".",first.reso[1],".dendrogam.pdf"),collapse=""),height=9,width=9)
	print(dendro.heap.p$part1.withoutdog[[1]])
	print(dendro.heap.p$part1.withoutdog[[2]])
	grid.newpage()
	grid.table(as.data.frame(total.dist.matrxes$Sample1cross2.centers$df2_df1.relation))
	grid.newpage()
	grid.draw(tableGrob(total.dist.matrxes$Sample1cross2.centers$seg.Data.main, theme = ttheme_default(base_size = 9), vp = NULL))
	print(total.relation.plots$p1)
	print(total.relation.plots$p2)
	dev.off()
	return(list(total.dist.matrxes=total.dist.matrxes,total.relation.plots=total.relation.plots,dendro.heap.p=dendro.heap.p))
}



#' Refine.relationship
#'
#' This is and internal function used in Tree.build.2nd.treemaking to refine the relationship if one upstreme nod is connected to more than one nods in the downstreme stage.  Input is Input\    deeper.dist.matrxes.current  and p.dendro
#' @param dist.matrix   deeper.dist.matrxes.current
#' @param all.matrix   p.dendro[[length(p.dendro)]]$matrix
#' @return   This will return an object used in Tree.build.2nd.treemaking.  return the new "XX.dist.matrxes.current",  assign it back
#' @export
#' @examples
#' Refine.relationship(deeper.dist.matrxes.current,p.dendro[[length(p.dendro)]]$matrix)

Refine.relationship<-function(dist.matrix=deeper.dist.matrxes.current,all.matrix=p.dendro[[length(p.dendro)]]$matrix)
{
	reverse.seg<-dist.matrix$Sample1cross2.centers$seg.Data.alt
	main.seg=dist.matrix$Sample1cross2.centers$seg.Data.main
	cutoff<-sum(all.matrix)/(nrow(all.matrix)*(ncol(all.matrix)-1))
	main.seg.new<-c()
	internal.seg.new<-c()
	for (seekFrom in as.character(main.seg$SeekFrom.cluster.names))  #  to loop around the seek from cluster(aka maturer cluster which supposedto find its upstreme in most cases)
	{
		current.toward.clust<-main.seg$SeekToward.cluster.names[main.seg$SeekFrom.cluster.names==seekFrom]   # To show the seek-toward clusters that is before refine. It couild be multiple seekToward
		if(length(which(main.seg$SeekToward.cluster.names==current.toward.clust))==1)    # If the seektoward is only one, which means only one line is connected. Then this will directly pass my refine.
		{
			main.seg.new<-rbind(main.seg.new,main.seg[which(main.seg$SeekToward.cluster.names==current.toward.clust),])
		}else          # else, firstly to check if the current seekfrom is also the closest in terms of reverse relationship, if it is then this relationship will be maintained
		{
			multi.strongest.current<-reverse.seg[reverse.seg$SeekFrom.cluster.names==as.character(current.toward.clust),]
			multi.strongest.current<-multi.strongest.current[,c(2,1,4,3)]
			names(multi.strongest.current)<-names(multi.strongest.current)[c(2,1,4,3)]
			if(seekFrom==as.character(multi.strongest.current$SeekFrom.cluster.names))   # This is the situation where tis current seekfrom is exactly the closest in terms of reverse relationship
			{
				main.seg.new<-rbind(main.seg.new,multi.strongest.current)
			}else if(min(all.matrix[,seekFrom][all.matrix[,seekFrom]>0])<=cutoff)     # If the current seekfrom cluster is not the closest then check out the all matrix, to see if the closest to this current seekfrom cluster is upstreme or internal
			{
				if( grepl(unique(main.seg$SeekToward.stage.names),names(all.matrix[,seekFrom][all.matrix[,seekFrom]>0])[all.matrix[,seekFrom][all.matrix[,seekFrom]>0]==min(all.matrix[,seekFrom][all.matrix[,seekFrom]>0])]))
				{
					main.seg.new<-rbind(main.seg.new,main.seg[which(main.seg$SeekFrom.cluster.names==seekFrom),]	)
				}else
				{
					internal.seg.current<-c(SeekToward.cluster.names=names(all.matrix[,seekFrom][all.matrix[,seekFrom]>0])[all.matrix[,seekFrom][all.matrix[,seekFrom]>0]==min(all.matrix[,seekFrom][all.matrix[,seekFrom]>0])],SeekFrom.cluster.names=seekFrom,SeekFrom.stage.names=strsplit(seekFrom,"_")[[1]][1],SeekToward.stage.names=strsplit(seekFrom,"_")[[1]][1])
					internal.seg.new<-rbind(internal.seg.new,internal.seg.current)
				}
			}
		}
	}
	row.names(internal.seg.new)<-NULL
	internal.seg.new<-as.data.frame(internal.seg.new)
	if(!all(dim(internal.seg.new)==0))
	{
		dist.matrix$Sample1cross2.centers$seg.Data.main.internal<-internal.seg.new
	}
	dist.matrix$Sample1cross2.centers$seg.Data.main<-main.seg.new
	return(dist.matrix)
}






#' tree.plot.v2
#'
#' This is and internal function used in Tree.build.2nd.treemaking To plot the tree
#' @param maturername =downstremename,
#' @param youngername =upstremename,
#' @param withsingledog =F,dendrodata,
#' @param matrixes =deeper.dist.matrxes.current,
#' @param print.later.names =T,
#' @param print.earlier.names =T,
#' @param print.stage.level =T,
#' @param inputorder.mature r=F,
#' @param inputorder.yonger =F,
#' @param order.maturer =NULL,
#' @param order.yonger =NULL,
#' @param textsize =2,
#' @param ballsize =5,
#' @param percent =2,
#' @param nocellnumber =F,
#' @param thecellnumber =cellnumber
#' @return   This will return an object used in Tree.build.2nd.treemaking.  return the new "XX.dist.matrxes.current",  assign it back
#' @export
#' @examples
#' Refine.relationship(deeper.dist.matrxes.current,p.dendro[[length(p.dendro)]]$matrix)

tree.plot.v2<-function(maturername=downstremename,youngername=upstremename,withsingledog=F,dendrodata,matrixes=deeper.dist.matrxes.current,print.later.names=T,print.earlier.names=T,print.stage.level=T,inputorder.maturer=F,inputorder.yonger=F,order.maturer=NULL,order.yonger=NULL,textsize=2,ballsize=5,percent=2,nocellnumber=F,thecellnumber=cellnumber)
	{
		strsplit(youngername,"_") %>%  unlist ->youngername
		strsplit(maturername,"_") %>%  unlist ->maturername
		#### First of all  I will fix the self loop problem here  on seg.Data.main.weak
		if(!is.null(matrixes$Sample1cross2.centers$seg.Data.main.weak))
		{
			loop.df<-c()
			for (i in 1:nrow(matrixes$Sample1cross2.centers$seg.Data.main.weak))
			{
				if(matrixes$Sample1cross2.centers$seg.Data.main.weak[i,2] %in% matrixes$Sample1cross2.centers$seg.Data.main.weak[,1])
				{
					if(which(matrixes$Sample1cross2.centers$seg.Data.main.weak[,1]==matrixes$Sample1cross2.centers$seg.Data.main.weak[i,2]) %>% matrixes$Sample1cross2.centers$seg.Data.main.weak[.,2]  %>% as.character(.)==as.character(matrixes$Sample1cross2.centers$seg.Data.main.weak[i,1]))
					loopcur<-c(as.character(matrixes$Sample1cross2.centers$seg.Data.main.weak[i,1]),as.character(matrixes$Sample1cross2.centers$seg.Data.main.weak[i,2]))
					loop.df<-rbind(loop.df,loopcur)
				}

			}
			if(!is.null(loop.df))
			{
				loop.df<-t(apply(loop.df,1,function(x){x[order(x)]})) %>% .[!duplicated(.)]  %>% as.matrix %>% t
				for(pair in 1: nrow(loop.df))
				{
					dendrodata$part2.withdog$matrix[loop.df[pair,1],!colnames(dendrodata$part2.withdog$matrix) %in% loop.df[pair,]]  %>% min -> A.out.dist
					dendrodata$part2.withdog$matrix[loop.df[pair,2],!colnames(dendrodata$part2.withdog$matrix) %in% loop.df[pair,]]  %>% min -> B.out.dist
					toOut<-c(A.out.dist,B.out.dist) %>% which.min    ## 1 is A ;  2 is B
					Out.seektoward<-dendrodata$part2.withdog$matrix[loop.df[pair,toOut],!colnames(dendrodata$part2.withdog$matrix) %in% loop.df[pair,]] %>% which.min %>% names
					matrixes$Sample1cross2.centers$seg.Data.main.weak[matrixes$Sample1cross2.centers$seg.Data.main.weak[,1]== loop.df[pair,toOut] ,2]<-Out.seektoward
				}
			}
		}
		### Do the same same on seg.Data.main.internal
		if(!is.null(matrixes$Sample1cross2.centers$seg.Data.main.internal))
		{
			loop.df<-c()
			for (i in 1:nrow(matrixes$Sample1cross2.centers$seg.Data.main.internal))
			{
				if(matrixes$Sample1cross2.centers$seg.Data.main.internal[i,2] %in% matrixes$Sample1cross2.centers$seg.Data.main.internal[,1])
				{
					if(which(matrixes$Sample1cross2.centers$seg.Data.main.internal[,1]==matrixes$Sample1cross2.centers$seg.Data.main.internal[i,2]) %>% matrixes$Sample1cross2.centers$seg.Data.main.internal[.,2]  %>% as.character(.)==as.character(matrixes$Sample1cross2.centers$seg.Data.main.internal[i,1]))
					{
						loopcur<-c(as.character(matrixes$Sample1cross2.centers$seg.Data.main.internal[i,1]),as.character(matrixes$Sample1cross2.centers$seg.Data.main.internal[i,2]))
						loop.df<-rbind(loop.df,loopcur)
					}
				}

			}
			if(!is.null(loop.df))
			{
				loop.df<-t(apply(loop.df,1,function(x){x[order(x)]})) %>% .[!duplicated(.)]  %>% as.matrix %>% t
				for(pair in 1: nrow(loop.df))
				{
					dendrodata$part2.withdog$matrix[loop.df[pair,1],!colnames(dendrodata[[length(dendrodata)]]$matrix) %in% loop.df[pair,]]  %>% min -> A.out.dist
					dendrodata$part2.withdog$matrix[loop.df[pair,2],!colnames(dendrodata[[length(dendrodata)]]$matrix) %in% loop.df[pair,]]  %>% min -> B.out.dist
					toOut<-c(A.out.dist,B.out.dist) %>% which.min    ## 1 is A ;  2 is B
					Out.seektoward<-dendrodata$part2.withdog$matrix[loop.df[pair,toOut],!colnames(dendrodata$part2.withdog$matrix) %in% loop.df[pair,]] %>% which.min %>% names
					matrixes$Sample1cross2.centers$seg.Data.main.internal[matrixes$Sample1cross2.centers$seg.Data.main.internal[,2]== loop.df[pair,toOut] ,1]<-Out.seektoward
				}
			}
		}
		### ------------------------------------------------------------------------------------------------------------------------------------------------ self loopo problem solved.
		if (!nocellnumber)
		{
			cellnumber.info=thecellnumber
		}
		#merge.end<-"end.cluster.names"
		#merge.start<-"start.cluster.names"
		#merge.end.stage<-"end.stage.names"
		#merge.start.stage<-"start.stage.names"
		#if(length(matrixes)==3)
		#{
		matrixes<-matrixes$Sample1cross2.centers
		names(matrixes)[3]<-"dist.matrix"
		names(matrixes)[8]<-"segData.2_1"  #  This was the seg.Data.main
		names(matrixes)[9]<-"segData.1_2"
		merge.end<-"SeekFrom.cluster.names"
		merge.start<-"SeekToward.cluster.names"
		merge.end.stage<-"SeekFrom.stage.names"
		merge.start.stage<-"SeekToward.stage.names"
		#}
		#Calculate the tree dataframe
		if(withsingledog)
		{
			tree.df<-data.frame(cluster.names=colnames(dendrodata$part2.withdog$matrix),order=1:length(colnames(dendrodata$part2.withdog$matrix)))   #  here the order is just to make a column. it will be corrected later
		}else
		{
			tree.df<-data.frame(cluster.names=colnames(dendrodata$part1.withoutdog$matrix),order=1:length(colnames(dendrodata$part1.withoutdog$matrix)))
		}
		tree.df<-cbind(tree.df,stage=unlist(lapply(strsplit(as.character(tree.df$cluster.names),"_"),function(x){x[1]})))
		tree.df.dn<-subset(tree.df,stage==maturername)
		tree.df.up<-subset(tree.df,stage!=maturername)
		tree.df.dn$order<-1:nrow(tree.df.dn)
		tree.df.up$order<-1:nrow(tree.df.up)
		if(inputorder.maturer)
		{
			tree.df.up$cluster.names=order.maturer
		}
		if(inputorder.yonger)
		{
			tree.df.dn$cluster.names=order.yonger
		}
		matrixes$tree.df<-rbind(tree.df.up,tree.df.dn)
		#To prepare the segData(the point to point link informations)
		matrixes$segData.2_1.ordered<-merge(matrixes$segData.2_1,matrixes$tree.df[,1:2],by.x=merge.end,by.y="cluster.names")  ##  order.x indicate the order for "seekFrom"
		matrixes$segData.2_1.ordered<-merge(matrixes$segData.2_1.ordered,matrixes$tree.df[,1:2],by.x=merge.start,by.y="cluster.names")  ##  order.x indicate the order for "seekToward"
		if(any(grepl("internal",names(matrixes))))
		{
			matrixes$segData.2_1.internal.ordered<-matrixes$seg.Data.main.internal[complete.cases(matrixes$seg.Data.main.internal),]
			if(is.atomic(matrixes$segData.2_1.internal.ordered))
			{
				matrixes$segData.2_1.internal.ordered<-rbind(matrixes$segData.2_1.internal.ordered,matrixes$segData.2_1.internal.ordered)
				matrixes$segData.2_1.internal.ordered<-merge(matrixes$segData.2_1.internal.ordered,matrixes$tree.df[,1:2],by.x=merge.end,by.y="cluster.names")
				matrixes$segData.2_1.internal.ordered<-merge(matrixes$segData.2_1.internal.ordered,matrixes$tree.df[,1:2],by.x=merge.start,by.y="cluster.names")
				matrixes$segData.2_1.internal.ordered<-unique(matrixes$segData.2_1.internal.ordered)
			}else
			{
				matrixes$segData.2_1.internal.ordered<-merge(matrixes$segData.2_1.internal.ordered,matrixes$tree.df[,1:2],by.x=merge.end,by.y="cluster.names")
				matrixes$segData.2_1.internal.ordered<-merge(matrixes$segData.2_1.internal.ordered,matrixes$tree.df[,1:2],by.x=merge.start,by.y="cluster.names")
			}

		}
		if(any(grepl("weak",names(matrixes))))
		{
			matrixes$segData.2_1.weak.ordered<-matrixes$seg.Data.main.weak[complete.cases(matrixes$seg.Data.main.weak),]
			if(is.atomic(matrixes$segData.2_1.weak.ordered))
			{
				matrixes$segData.2_1.weak.ordered<-rbind(matrixes$segData.2_1.weak.ordered,matrixes$segData.2_1.weak.ordered)
				matrixes$segData.2_1.weak.ordered<-merge(matrixes$segData.2_1.weak.ordered,matrixes$tree.df[,1:2],by.x=merge.end,by.y="cluster.names")
				matrixes$segData.2_1.weak.ordered<-merge(matrixes$segData.2_1.weak.ordered,matrixes$tree.df[,1:2],by.x=merge.start,by.y="cluster.names")
				matrixes$segData.2_1.weak.ordered<-unique(matrixes$segData.2_1.weak.ordered)
			}else
			{
				matrixes$segData.2_1.weak.ordered<-merge(matrixes$segData.2_1.weak.ordered,matrixes$tree.df[,1:2],by.x=merge.end,by.y="cluster.names")
				matrixes$segData.2_1.weak.ordered<-merge(matrixes$segData.2_1.weak.ordered,matrixes$tree.df[,1:2],by.x=merge.start,by.y="cluster.names")
			}
		}
		matrixes$segData.1_2.ordered<-merge(matrixes$segData.1_2,matrixes$tree.df[,1:2],by.x=merge.end,by.y="cluster.names")
		matrixes$segData.1_2.ordered<-merge(matrixes$segData.1_2.ordered,matrixes$tree.df[,1:2],by.x=merge.start,by.y="cluster.names")
		#To store the stage level into stage.levels. the mature group order; the yonger group order
		stage.level<-c(maturername,youngername)
		matrixes$tree.df$stage<-factor(matrixes$tree.df$stage,levels=stage.level)
		yonger.order<-as.character(tree.df.up$cluster.names)
		maturer.order<-as.character(tree.df.dn$cluster.names)
		matrixes$tree.df<-cbind(matrixes$tree.df,First.cluster.tag=unlist(lapply(strsplit(as.character(matrixes$tree.df$cluster.names),"_"),function(x){x[2]})))  #  Got a tag to tell the first cluster, which will be coded by color.
		# Print out the order information
		if(print.later.names)
		{
			print("Here are the ordered cluster names that are relatively maturer")
			print(maturer.order)
			print("")
		}
		if(print.earlier.names)
		{
			print("Here are the ordered cluster names that are relatively yonger\n")
			print(yonger.order)
			print("")
		}
		if(print.stage.level)
		{
			print("Here is the level of stages\n:The first is on bottom, the laste is on top")
			print(stage.level)
			print("")
		}
		if (!nocellnumber)
		{
			matrixes$tree.df<-merge(matrixes$tree.df,cellnumber.info[,c(1,2,4)],by.x="cluster.names",by.y="Var1")
		}
		p1<-ggplot(matrixes$tree.df)+aes(order,stage,color=First.cluster.tag,size=percent)+geom_point()+geom_text(aes(label=cluster.names),vjust=1.5,hjust="left",color="blue",size=textsize,fontface="bold",nudge_y=0.1)+geom_segment(data=matrixes$segData.2_1.ordered,aes_string(x="order.x",xend="order.y",y=merge.end.stage,yend=merge.start.stage),arrow=arrow(),lwd=1,col="red")+xlim(1,max(matrixes$tree.df$order)+1)
		if(any(grepl("internal",names(matrixes))))
		{
			p1<-p1+geom_curve(data=matrixes$segData.2_1.internal.ordered,aes_string(x="order.x",xend="order.y",y=merge.end.stage,yend=merge.start.stage),arrow=arrow(length=unit(0.03, "npc"),ends="last", type = "closed"),lwd=0.5,linetype=1,col="black")
		}
		if(any(grepl("weak",names(matrixes))))
		{
			p1<-p1+geom_curve(data=matrixes$segData.2_1.weak.ordered,aes_string(x="order.x",xend="order.y",y=merge.end.stage,yend=merge.start.stage),arrow=arrow(length=unit(0.03, "npc"),ends="last", type = "closed"),lwd=0.5,linetype=2,col="black")
		}
		p2<-ggplot(matrixes$tree.df)+aes(order,stage,color=First.cluster.tag)+geom_point(size=8)+geom_text(aes(label=cluster.names),vjust=1.5,hjust="left",color="blue",size=textsize,fontface="bold",nudge_y=0.1)+geom_segment(data=matrixes$segData.1_2.ordered,aes_string(x="order.x",xend="order.y",y=merge.end.stage,yend=merge.start.stage),arrow=arrow(),lwd=1,col="black")+xlim(1,max(matrixes$tree.df$order)+1)
		return(list(p1=p1,p2=p2,matrix.advance=matrixes))
	}


	###############################################  Tree.build.2nd.clustering     #################################################################
	### Introduction\   This function did seurat primary clustering for the second layer
	### Input\   tree.prep,tree.1.ob, which is the output from Tree.build.1;  the secondary resolution for clustering
	### Output\   primaries.deeper.lst which is a list primary objects, including multiple pairs of connections
	###  Example\  s4_s5.tree.2nd.primary.list<-Tree.build.2nd.clustering(s4_s5.tree.prep.test,s4_s5.tree.1.test,second.reso=c(0.3,0.3))
	### Function_define



#' Tree.build.2nd.clustering
#'
#' This function did seurat primary clustering for the second layer,  which is higher resolution . Input\   tree.prep,tree.1.ob, which is the output from Tree.build.1;  thesecondary resolution for clustering. ### Output\   primaries.deeper.lst which is a list primary objects, including multiple pairs of connections.
#' @param tree.prep This is an object generated from Tree.build.prepare or Prep.res.adjust.  This is basically the clustering result.
#' @param tree.1.ob This is an object generated from Tree.build.1.  This is the low resulotion subtree
#' @param second.reso This is the resolution for the high resolution.  For example, c(0.3,0.3)
#' @return   This will return the secondary clustering result. It will be used for Tree.build.2nd.clustering.patch and Tree.build.2nd.treemaking
#' @export
#' @examples
#' S4_S5ROCK.tree.2nd.primary_0.3_0.3.list<-Tree.build.2nd.clustering(S4_S5ROCK.tree.prep,S4_S5ROCK.tree.1.ob,second.reso=c(0.3,0.3))

Tree.build.2nd.clustering<-function(tree.prep,tree.1.ob,second.reso=c(0.3,0.3))
{

		primaries.deeper.lst<-list()
		all.list.names<-c()
		for (i in 1:nrow(tree.1.ob$total.dist.matrxes$Sample1cross2.centers$seg.Data.main))
		{
			cur.list.name<-paste(tree.1.ob$total.dist.matrxes$Sample1cross2.centers$seg.Data.main$SeekToward.cluster.names[i],tree.1.ob$total.dist.matrxes$Sample1cross2.centers$seg.Data.main$SeekFrom.cluster.names[i],sep="_")
			all.list.names<-c(all.list.names,cur.list.name)
			# To assign the dge data
			dge.pair.lst<-list(as.matrix(tree.prep$primary.total$Seurat.list[[1]]@raw.data)[,row.names(tree.prep$primary.total$Seurat.list[[1]]@data.info)[which(tree.prep$primary.total$Seurat.list[[1]]@data.info$Sample.2nd==as.character(tree.1.ob$total.dist.matrxes$Sample1cross2.centers$seg.Data.main[i,2]))]],as.matrix(tree.prep$primary.total$Seurat.list[[2]]@raw.data)[,row.names(tree.prep$primary.total$Seurat.list[[2]]@data.info)[which(tree.prep$primary.total$Seurat.list[[2]]@data.info$Sample.2nd==as.character(tree.1.ob$total.dist.matrxes$Sample1cross2.centers$seg.Data.main[i,1]))]])
			names(dge.pair.lst)<-c(as.character(tree.1.ob$total.dist.matrxes$Sample1cross2.centers$seg.Data.main[i,2]),as.character(tree.1.ob$total.dist.matrxes$Sample1cross2.centers$seg.Data.main[i,1]))
			primary.current<-Primaryclutering(dge.pair.lst,txcut=500,cluster.resos=second.reso,RAWinput=F)
			primaries.deeper.lst<-c(primaries.deeper.lst,list(primary.current))
		}
		names(primaries.deeper.lst)<-all.list.names
		###Print  information about how I did the secondary clustering
		print (paste(c("There are ",length(primaries.deeper.lst)," pairs of connection branches"),collapse=""))
		return(primaries.deeper.lst=primaries.deeper.lst)
}


	###############################################  Tree.build.2nd.clustering.patch     #################################################################
	### Introduction\   This function is to address single dog problem. It will Identify and do seurat clustering for the single dog & dog+neighbours
	### Input\    the clustered result from Tree.build.2nd.clustering as well as tree.prep  and tree.1.ob. singledog.reso should be specified
	### Output\     If there was no single dog, the out put is exactly the same as imput,  otherwise, two seurat object would be added at the end, forming a list as whole
	###  Example\   S2_S2.24.tree.2nd.primary_0.06.list.patched<-Tree.build.2nd.clustering.patch(S2_S2.24.tree.2nd.primary_0.06.list,S2_S2.24.tree.prep,S2_S2.24.tree.1.ob,singledog.reso=0.06)
	### Function_define


#' Tree.build.2nd.clustering.patch
#'
#' Introduction\   This function is to address single dog problem. It will Identify and do seurat clustering for the single dog & dog+neighbours. ### Input\    the clustered result from Tree.build.2nd.clustering as well as tree.prep  and tree.1.ob. singledog.reso should be specified. ### Output\     If there was no single dog, the out put is exactly the same as imput,  otherwise, two seurat object would be added at the end, forming a list as whole
#' @param primaries.deeper.lst   This is generated by Tree.build.2nd.clustering and Second.primary.adjust
#' @param tree.prep  This is the low resolution clustering result from Tree.build.prepare
#' @param tree.1.ob This is the low reslution subtree result from Tree.build.1
#' @param singledog.reso  This is the resolution singledog cluster for 0.3
#' @return   This will return the high resolution result with a pach that consider the single dog.  This is going to be used for Tree.build.2nd.treemaking
#' @export
#' @examples
#' S4_S5ROCK.tree.2nd.primary_0.3_0.3.list<-Tree.build.2nd.clustering(S4_S5ROCK.tree.prep,S4_S5ROCK.tree.1.ob,second.reso=c(0.3,0.3))


Tree.build.2nd.clustering.patch<-function(primaries.deeper.lst,tree.prep,tree.1.ob,singledog.reso=0.3)
{
		Singledog.clusters<-setdiff(tree.1.ob$total.relation.plots$matrix.advance$tree.df$cluster.names[as.character(tree.1.ob$total.relation.plots$matrix.advance$tree.df$stage)==as.character(tree.1.ob$total.relation.plots$matrix.advance$segData.2_1.ordered$SeekFrom.stage.names)],as.character(tree.1.ob$total.relation.plots$matrix.advance$segData.2_1.ordered$SeekFrom.cluster.names))  # To get the cluster name who has no any relationship with upstream clusters.  I will also do clustering for it AND seperartely store it in single.seurat.list
		#	single.seurat.list<-c()   ###  I annotated this part out because most commonly there would be only one single dog cluster.
		#	for(i in length(Singledog.clusters))
		#	{
		#
		#	}
		print(paste("Number of single dog is",length(Singledog.clusters)))
		if (length(Singledog.clusters)>0)
		{
			singledog.dge<-as.matrix(tree.prep$primary.total$Seurat.list[[2]]@raw.data[,row.names(tree.prep$primary.total$Seurat.list[[2]]@data.info)[tree.prep$primary.total$Seurat.list[[2]]@data.info$Sample.2nd==Singledog.clusters]])
			singledog.seurat.ob<-docluster(dgepreprocess(singledog.dge,500,norowname=F),GetinformativeGene(dgepreprocess(singledog.dge,500,norowname=F),500),Singledog.clusters,reso=singledog.reso)

			neighbourWithsingledog.seurat.ob<-docluster.multi(500,sets=list(as.matrix(primaries.deeper.lst[[1]]$Seurat.list[[1]]@raw.data),as.matrix(primaries.deeper.lst[[1]]$Seurat.list[[2]]@raw.data),singledog.dge),nms=c(names(primaries.deeper.lst[[1]]$Seurat.list)[1:2],Singledog.clusters))
			singledog.twoobs<-list(singledog.seurat.ob=singledog.seurat.ob,neighbourWithsingledog.seurat.ob=neighbourWithsingledog.seurat.ob)
			primaries.deeper.lst<-c(primaries.deeper.lst,singledog=singledog.twoobs)
		}
		return(primaries.deeper.lst)
}




#' Second.primary.adjust
#'
#' This is to adjust the clustering resolution for secondary clusters.
#' @param second.pr This is the secondary clustering object from Tree.build.2nd.clustering
#' @param prep this is the low resolution clutsering object. from Tree.build.prepare
#' @param second.reso.adjust  This is the adjusted resolution    for example: c(0.3,0.6)
#' @return   This will return.   The adjusted secondary clustering, assign back to the same object generated from Tree.build.2nd.clustering
#' @export
#' @examples
#' S4_S5ROCK.tree.2nd.primary_0.3_0.3.list<-Second.primary.adjust(S4_S5ROCK.tree.2nd.primary_0.3_0.3.list,S4_S5ROCK.tree.prep,second.reso.adjust=c(0.3,0.6))

Second.primary.adjust<-function(second.pr,prep,second.reso.adjust=c(0.3,0.6))
{
	setwd(prep$path)
	name1<-names(prep$primary.total$Seurat.list)[1]
	name2<-names(prep$primary.total$Seurat.list)[2]
	setwd(paste(c("Hires",name1,name2),collapse="."))
	for (i in 1: length(second.pr))
	{
		second.pr[[i]]$Seurat.list[[1]]<-FindClusters(second.pr[[i]]$Seurat.list[[1]], pc.use = 1:10, resolution = second.reso.adjust[1], print.output = 0, save.SNN = T)
		second.pr[[i]]$Seurat.list[[1]]@data.info<-second.pr[[i]]$Seurat.list[[1]]@data.info[,c("nGene","nUMI","orig.ident","percent.mito",paste("res",second.reso.adjust[1],sep="."),"Sample")] %>% cbind(.,Sample.2nd=paste(.$Sample,.[,paste("res",second.reso.adjust[1],sep=".")],sep="_"))

		second.pr[[i]]$Seurat.list[[2]]<-FindClusters(second.pr[[i]]$Seurat.list[[2]], pc.use = 1:10, resolution = second.reso.adjust[2], print.output = 0, save.SNN = T)
		second.pr[[i]]$Seurat.list[[2]]@data.info<-second.pr[[i]]$Seurat.list[[2]]@data.info[,c("nGene","nUMI","orig.ident","percent.mito",paste("res",second.reso.adjust[2],sep="."),"Sample")] %>% cbind(.,Sample.2nd=paste(.$Sample,.[,paste("res",second.reso.adjust[2],sep=".")],sep="_"))
	}
	return(second.pr)
}




#' DOCLEAN
#'
#' This is to delete the past files with not-so-good resolution.
#' @param prep  This is the object from Tree.build.prepare
#' @param hiorlo   if == "hi", then go delete the files in hi resolution folder. if == "lo" then go to low resolution folder to delete
#' @param res.tokeep   This is th3e final resolution to keep. for example c(0.3,0.6)
#' @return   No return
#' @export
#' @examples
#' S4_S5ROCK.tree.2nd.primary_0.3_0.3.list<-Second.primary.adjust(S4_S5ROCK.tree.2nd.primary_0.3_0.3.list,S4_S5ROCK.tree.prep,second.reso.adjust=c(0.3,0.6))

DOCLEAN<-function(prep,hiorlo="hi",res.tokeep=c(0.3,0.6)){
	setwd(prep$path)
	name1<-names(prep$primary.total$Seurat.list)[1]
	name2<-names(prep$primary.total$Seurat.list)[2]
	if(hiorlo=="hi"){
		setwd(paste(c("Hires",name1,name2),collapse="."))
		c(grep(paste(res.tokeep,collapse=" _ "),list.files()),grep(paste(name1,"...",res.tokeep[1],sep=""),list.files()),grep(paste(name2,"...",res.tokeep[2],sep=""),list.files()))->indextokeep
    print(paste("To removing following",list.files()[setdiff(1:length(list.files()),indextokeep)]))
		file.remove(list.files()[setdiff(1:length(list.files()),indextokeep)])

	}else if(hiorlo=="lo"){
		setwd(paste(c("LowresCluster",name1,name2),collapse="."))
		c(grep(paste(c(name1,"total",res.tokeep[1]),collapse="."),list.files()), grep(paste(c(name2,"total",res.tokeep[2]),collapse="."),list.files()), grep(paste(c(name1,res.tokeep[1]),collapse="."),list.files()))->indextokeep
    print(paste("To removing following",list.files()[setdiff(1:length(list.files()),indextokeep)]))
		file.remove(list.files()[setdiff(1:length(list.files()),indextokeep)])

	}
}



#' Tree.build.2nd.treemaking
#'
#' Introduction\   Similar with Tree.build.1,  this function after getting the clustering primary object for the second layer. this one is literally building the tree and generate plots
#' @param primaries.deeper.lst  This is the adjusted secondary  clustering results from Tree.build.2nd.clustering.patch
#' @param second.reso   This should be  consistent with what eventually used in secondary clustering,  it should follow Second.primary.adjust
#' @param singledog.reso   This should follow that in Tree.build.2nd.clustering.patch
#' @param downstremename   ="s1.24"
#' @param upstremename   ="s1"
#' @param relationtext   =2.5
#' @param dorefine   =T
#' @param dir  The path of the whole tree generating.  keep consistent with Tree.build.prepare
#' @return   No return
#' @export
#' @examples
#' S4_S5ROCK.tree.2nd.primary_0.3_0.3.list<-Second.primary.adjust(S4_S5ROCK.tree.2nd.primary_0.3_0.3.list,S4_S5ROCK.tree.prep,second.reso.adjust=c(0.3,0.6))

Tree.build.2nd.treemaking<-function(primaries.deeper.lst,second.reso=c(0.3,0.3),singledog.reso=0.06,downstremename="s1.24",upstremename="s1",relationtext=2.5,dorefine=T,dir){
		setwd(paste(dir,paste(upstremename,downstremename,sep="_"),sep=""))
		if(!any(grepl(paste("Hires",upstremename,downstremename,sep="."),list.files()))){
			dir.create(paste("Hires",upstremename,downstremename,sep="."))
		}
		setwd(paste("Hires",upstremename,downstremename,sep="."))
		deeper.dist.matrxes.lst<-list()
		deeper.relation.plots.lst<-list()
		deeper.fullplots.1.lst<-list()
		deeper.fullplots.2.lst<-list()
		deeper.dendro.lst<-list()
		deeper.heat.lst<-list()
		deeper.allmatrix.lst<-list()
		thereisdog<-any(grepl("singledog",names(primaries.deeper.lst)))
		for (i in which(!grepl("singledog",names(primaries.deeper.lst))))
		{
			print(paste(c(i,":", names(primaries.deeper.lst[[i]]$Seurat.list)),collapse=" "))
			primaries.deeper.lst[[i]]$Seurat.list[[1]]@data.info<-primaries.deeper.lst[[i]]$Seurat.list[[1]]@data.info[,1:6]  #  This is to reset
			primaries.deeper.lst[[i]]$Seurat.list[[2]]@data.info<-primaries.deeper.lst[[i]]$Seurat.list[[2]]@data.info[,1:6]
			primaries.deeper.lst[[i]]$Seurat.list[[1]]@data.info<-cbind(primaries.deeper.lst[[i]]$Seurat.list[[1]]@data.info,Sample.2nd=paste(primaries.deeper.lst[[i]]$Seurat.list[[1]]@data.info$Sample,primaries.deeper.lst[[i]]$Seurat.list[[1]]@data.info[,paste("res.",second.reso[1],sep="")],sep="_"))   # give a
			primaries.deeper.lst[[i]]$Seurat.list[[2]]@data.info<-cbind(primaries.deeper.lst[[i]]$Seurat.list[[2]]@data.info,Sample.2nd=paste(primaries.deeper.lst[[i]]$Seurat.list[[2]]@data.info$Sample,primaries.deeper.lst[[i]]$Seurat.list[[2]]@data.info[,paste("res.",second.reso[2],sep="")],sep="_"))
			if(thereisdog)
			{
				primaries.deeper.lst$singledog.singledog.seurat.ob@data.info<-primaries.deeper.lst$singledog.singledog.seurat.ob@data.info[,1:6]
				primaries.deeper.lst$singledog.singledog.seurat.ob@data.info<-cbind(primaries.deeper.lst$singledog.singledog.seurat.ob@data.info,Sample.2nd=paste(primaries.deeper.lst$singledog.singledog.seurat.ob@data.info$Sample,primaries.deeper.lst$singledog.singledog.seurat.ob@data.info[,paste("res.",singledog.reso,sep="")],sep="_"))
				print("Printing single dog fullplot into pdf")
				if(length(unique(primaries.deeper.lst$singledog.singledog.seurat.ob@data.info[,5]))>1)
				{
					fullplot.current.dog<-Fullplot_v2(primaries.deeper.lst$singledog.singledog.seurat.ob,paste(unique(primaries.deeper.lst$singledog.singledog.seurat.ob@data.info$Sample),singledog.reso,".pdf",collapse=""),signiture=NULL,resolusion=paste("res.",singledog.reso,sep=""),heatmapannosize=0.5,doreturn=T)
				}else
				{
					fullplot.current.dog<-Fullplot_v2(primaries.deeper.lst$singledog.singledog.seurat.ob,paste(unique(primaries.deeper.lst$singledog.singledog.seurat.ob@data.info$Sample),singledog.reso,".pdf",collapse=""),signiture=NULL,resolusion=paste("res.",singledog.reso,sep=""),heatmapannosize=0.5,doreturn=T)
				}
			}
			# Print out the cluster
			print("Start to plot clustering into pdf...")
			fullplot.current.1<-Fullplot_v2(primaries.deeper.lst[[i]]$Seurat.list[[1]],paste(names(primaries.deeper.lst[[i]]$Seurat.list)[1],second.reso[1],".pdf",collapse=""),signiture=NULL,resolusion=paste("res.",second.reso[1],sep=""),heatmapannosize=0.5,doreturn=T)
			fullplot.current.2<-Fullplot_v2(primaries.deeper.lst[[i]]$Seurat.list[[2]],paste(names(primaries.deeper.lst[[i]]$Seurat.list)[2],second.reso[2],".pdf",collapse=""),signiture=NULL,resolusion=paste("res.",second.reso[2],sep=""),heatmapannosize=0.5,doreturn=T)
			deeper.fullplots.1.lst<-c(deeper.fullplots.1.lst,list(fullplot.current.1))
			deeper.fullplots.2.lst<-c(deeper.fullplots.2.lst,list(fullplot.current.2))
			deeper.dist.matrxes.current<-dist.matrix.prep(primaries.deeper.lst[[i]],datainfo.col1=c(paste("res.",second.reso[1],sep=""),"nUMI","Sample"),cluster.col1=paste("res.",second.reso[1],sep=""),datainfo.col2=c(paste("res.",second.reso[2],sep=""),"nUMI","Sample"),cluster.col2=paste("res.",second.reso[2],sep=""),res1=paste("res.",second.reso[1],sep=""),res2=paste("res.",second.reso[2],sep=""))
			# This gave back distance matrix as well as segData that is critical for relationship drawing.  and importantly, the segData determines what will be further explored.  I reccomend manual inspection although integrated should be able to well tell the relationship
			if(thereisdog)
			{
				p.dendro<-dist.matrix.prep.v3(binary.primary=primaries.deeper.lst[[i]],datainfo.col1=c(paste("res.",second.reso[1],sep=""),"nUMI","Sample"),datainfo.col2=c(paste("res.",second.reso[2],sep=""),"nUMI","Sample"),cluster.col1=paste("res.",second.reso[1],sep=""),cluster.col2=paste("res.",second.reso[2],sep=""),cluster.col.dog=paste("res.",singledog.reso,sep=""),datainfo.col.dog=c(paste("res.",singledog.reso,sep="")),dog.seurat.ob=primaries.deeper.lst$singledog.singledog.seurat.ob,neighborwithdog.seurat.ob=primaries.deeper.lst$singledog.neighbourWithsingledog.seurat.ob,withsingledog=thereisdog)
			}else
			{
				p.dendro<-dist.matrix.prep.v3(binary.primary=primaries.deeper.lst[[i]],datainfo.col1=c(paste("res.",second.reso[1],sep=""),"nUMI","Sample"),datainfo.col2=c(paste("res.",second.reso[2],sep=""),"nUMI","Sample"),cluster.col1=paste("res.",second.reso[1],sep=""),cluster.col2=paste("res.",second.reso[2],sep=""))
			}
			#  to creat internal relationship, or (if not < cutoff) weak relationship which could be inter or intra for unpaired clusters
			unpaired.cluster.names<-as.character(setdiff(colnames(p.dendro$part1.withoutdog$matrix)[grepl(downstremename,colnames(p.dendro$part1.withoutdog$matrix))],deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main$SeekFrom.cluster.names))
			deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main[,1]<-factor(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main[,1],levels=colnames(p.dendro$part1.withoutdog$matrix))
			deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main[,2]<-factor(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main[,2],levels=colnames(p.dendro$part1.withoutdog$matrix))
			deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.internal<-c()
			deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak<-c()
			if(dorefine)
			{
				deeper.dist.matrxes.current<-Refine.relationship(deeper.dist.matrxes.current,p.dendro[[length(p.dendro)]]$matrix)
			}	#dist.cutoff<-mean(p.dendro$part1.withoutdog$matrix)   #  Take the whole matrix, including 0
			dist.cutoff<-mean(p.dendro$part1.withoutdog$matrix[p.dendro$part1.withoutdog$matrix>0])   #  Calculate a distance cutoff, out of which will be considered no connection
			for(j in 1:4)
			{
				deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.internal[,j]<-as.character(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.internal[,j])
			}
			for (unclust in unpaired.cluster.names)
			{
				current.unclust.distances<-p.dendro$part1.withoutdog$matrix[unclust,][p.dendro$part1.withoutdog$matrix[unclust,]>0]   #  a vector of distances from this uncluster to either other cluster in this stage or the one from corresponding upstreme stage clusters(accept itself)
				current.unclust.distances.intra<-current.unclust.distances[grepl(downstremename,names(current.unclust.distances))]  #  Only the same stage distances are left
				if(min(current.unclust.distances.intra)<dist.cutoff)
				{
					internaltarget.current<-names(current.unclust.distances.intra)[which(current.unclust.distances.intra==min(current.unclust.distances.intra))]
					deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.internal<-rbind(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.internal,c(SeekToward.cluster.names=internaltarget.current,SeekFrom.cluster.names=unclust,SeekFrom.stage.names=unlist(strsplit(unclust,"_"))[1],SeekToward.stage.names=unlist(strsplit(internaltarget.current,"_"))[1]))
				}else
				{
					weaktarget.current<-names(current.unclust.distances)[which(current.unclust.distances==min(current.unclust.distances))]
					deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak<-rbind(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak,c(SeekFrom.cluster.names=unclust,SeekToward.cluster.names=weaktarget.current,SeekFrom.stage.names=unlist(strsplit(unclust,"_"))[1],SeekToward.stage.names=unlist(strsplit(weaktarget.current,"_"))[1]))
				}
			}
			if(thereisdog)
			{
				dist.cutoff<-mean(p.dendro$part2.withdog$matrix[p.dendro$part2.withdog$matrix>0])
				for(doggy in as.character(unique(primaries.deeper.lst$singledog.singledog.seurat.ob@data.info$Sample.2nd)))
				{
					current.doggy.distances<-p.dendro$part2.withdog$matrix[doggy,][p.dendro$part2.withdog$matrix[doggy,]>0]   #  a vector of distances from this uncluster to either other cluster in this stage or the one from corresponding upstreme stage clusters(accept itself)
					current.doggy.distances.intra<-current.doggy.distances[grepl(downstremename,names(current.doggy.distances))]  #  Only the same stage distances are left
					if(grepl(upstremename,names(which(current.doggy.distances==min(current.doggy.distances)))))   # If the minimum relationship is from the previous stage
					{
						if(min(current.doggy.distances)<dist.cutoff)      #  if the minimum relationship is from the the same stage stage
						{
							target.current<-names(current.doggy.distances)[which(current.doggy.distances==min(current.doggy.distances))]
							deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main$SeekFrom.cluster.names<-factor(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main$SeekFrom.cluster.names,levels=unique(c(levels(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main$SeekFrom.cluster.names),doggy)))
							deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main<-rbind(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main,c(SeekFrom.cluster.names=doggy,SeekToward.cluster.names=target.current,SeekFrom.stage.names=unlist(strsplit(doggy,"_"))[1],SeekToward.stage.names=unlist(strsplit(target.current,"_"))[1]))
						}else
						{
							if(nrow(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak)>1)
							{
								deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak$SeekFrom.cluster.names<-factor(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak$SeekFrom.cluster.names,levels=unique(c(levels(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak$SeekFrom.cluster.names),doggy)))
							}
							weaktarget.current<-names(current.doggy.distances)[which(current.doggy.distances==min(current.doggy.distances))]
							deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak<-rbind(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak,c(SeekFrom.cluster.names=doggy,SeekToward.cluster.names=weaktarget.current,SeekFrom.stage.names=unlist(strsplit(doggy,"_"))[1],SeekToward.stage.names=unlist(strsplit(weaktarget.current,"_"))[1]))
						}
					}else
					{
						if(min(current.doggy.distances.intra)<dist.cutoff)      #  if the minimum relationship is from the the same stage stage
						{
							internaltarget.current<-names(current.doggy.distances.intra)[which(current.doggy.distances.intra==min(current.doggy.distances.intra))]
							deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.internal<-rbind(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.internal,c(SeekFrom.cluster.names=doggy,SeekToward.cluster.names=internaltarget.current,SeekFrom.stage.names=unlist(strsplit(doggy,"_"))[1],SeekToward.stage.names=unlist(strsplit(internaltarget.current,"_"))[1]))
						}else
						{
							weaktarget.current<-names(current.doggy.distances)[which(current.doggy.distances==min(current.doggy.distances))]
							deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak<-rbind(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main.weak,c(SeekFrom.cluster.names=doggy,SeekToward.cluster.names=weaktarget.current,SeekFrom.stage.names=unlist(strsplit(doggy,"_"))[1],SeekToward.stage.names=unlist(strsplit(weaktarget.current,"_"))[1]))
						}
					}

				}
			}
			## Here,  to summarize the cell numbers
			cellnumber<-c()
			for (x in which(!grepl("singledog.neighbourWithsingledog",names(primaries.deeper.lst))))
			{
				if(grepl("singledog",names(primaries.deeper.lst)[x]))
				{
					cellnumber<-rbind(cellnumber,as.data.frame(table(primaries.deeper.lst[[x]]@data.info$Sample.2nd)))
				}else
				{
					cellnumber<-rbind(cellnumber,as.data.frame(table(primaries.deeper.lst[[x]]$Seurat.list[[1]]@data.info$Sample.2nd)))
					cellnumber<-rbind(cellnumber,as.data.frame(table(primaries.deeper.lst[[x]]$Seurat.list[[2]]@data.info$Sample.2nd)))
				}
			}
			cellnumber<-cbind(cellnumber,stage=unlist(lapply(strsplit(as.character(cellnumber$Var1),"_"),function(x){x[[1]]})))
			cellnumber<-mydplyr.percentage(cellnumber,by="stage")
			deeper.relation.plots<-tree.plot.v2(withsingledog=thereisdog,dendrodata=p.dendro,deeper.dist.matrxes.current,textsize=relationtext,maturername=downstremename,youngername=upstremename,thecellnumber=cellnumber)   #  This will give back a relationship network (low resolution)
			deeper.dist.matrxes.lst<-c(deeper.dist.matrxes.lst,list(deeper.dist.matrxes.current))
			deeper.relation.plots.lst<-c(deeper.relation.plots.lst,list(deeper.relation.plots))
			names(deeper.dist.matrxes.lst)[i]<-paste(c(names(primaries.deeper.lst[[i]]$Seurat.list)),collapse="__")
			names(deeper.relation.plots.lst)[i]<-paste(c(names(primaries.deeper.lst[[i]]$Seurat.list)),collapse="__")
			pdf(paste(paste(c(names(primaries.deeper.lst[[i]]$Seurat.list)),collapse="__"),second.reso[1],"_",second.reso[2],"dendrogam.pdf",collapse=""),height=9,width=9)

			if(thereisdog)
			{
				print(p.dendro$part2.withdog$p.heat.dog)
				print(p.dendro$part2.withdog$p.dendro.dog)
			}else
			{
				print(p.dendro$part1.withoutdog$p.heat)
				print(p.dendro$part1.withoutdog$p.dendro)
			}
			grid.newpage()
			grid.table(as.data.frame(deeper.dist.matrxes.current$Sample1cross2.centers$df1_df2.relation))
			grid.newpage()
			grid.draw(tableGrob(deeper.dist.matrxes.current$Sample1cross2.centers$seg.Data.main, theme = ttheme_default(base_size = 9), vp = NULL))
			print(deeper.relation.plots$p1)
			print(deeper.relation.plots$p2)
			dev.off()
			if(thereisdog)
			{
				deeper.dendro.lst<-c(deeper.dendro.lst,list(p.dendro$part2.withdog$p.dendro.dog))
				deeper.heat.lst<-c(deeper.heat.lst,list(p.dendro$part2.withdog$p.heat.dog))
				deeper.allmatrix.lst<-c(deeper.allmatrix.lst,list(p.dendro$part2.withdog$p.heat.dog$matrix))
			}else
			{
				deeper.dendro.lst<-c(deeper.dendro.lst,list(p.dendro$part1.withoutdog$p.dendro))
				deeper.heat.lst<-c(deeper.heat.lst,list(p.dendro$part1.withoutdog$p.heat))
				deeper.allmatrix.lst<-c(deeper.allmatrix.lst,list(p.dendro$part1.withoutdog$p.heat$matrix))
			}
			###
			#pdf(paste(c(name1,name2,"tree.pdf"),collapse="_"),height=14,width=14)
			#grid.arrange(total.fullplot.1[[1]],total.fullplot.2[[1]],dendro.heap.p[[1]],dendro.heap.p[[2]],total.relation.plots$p1,total.relation.plots$p2)
			#for (i in 1:length(primaries.deeper.lst))
			#{
			#grid.arrange(deeper.fullplots.1.lst[[i]][[1]],deeper.fullplots.2.lst[[i]][[1]],deeper.dendro.lst[[i]],deeper.heat.lst[[i]],deeper.relation.plots.lst[[i]][[1]],deeper.relation.plots.lst[[i]][[2]])
			#}
			#dev.off()
			##  For ones with S5 involved.  DO NOT DO IT
			#tm<-deeper.dist.matrxes.lst[[1]]
			#deeper.dist.matrxes.lst[[1]]<-deeper.dist.matrxes.lst[[2]]
			#deeper.dist.matrxes.lst[[2]]<-tm
			#tm<-deeper.relation.plots.lst[[1]]
			#deeper.relation.plots.lst[[1]]<-deeper.relation.plots.lst[[2]]
			#deeper.relation.plots.lst[[2]]<-tm
		}
		if(thereisdog)
		{
			return(list(deeper.dist.matrxes.lst=deeper.dist.matrxes.lst,deeper.relation.plots.lst=deeper.relation.plots.lst,deeper.fullplots.1.lst=deeper.fullplots.1.lst,deeper.fullplots.2.lst=deeper.fullplots.2.lst,deeper.dendro.lst=deeper.dendro.lst,deeper.heat.lst=deeper.heat.lst,fullplot.current.dog=fullplot.current.dog))
		}else
		{
			return(list(deeper.dist.matrxes.lst=deeper.dist.matrxes.lst,deeper.relation.plots.lst=deeper.relation.plots.lst,deeper.fullplots.1.lst=deeper.fullplots.1.lst,
				deeper.fullplots.2.lst=deeper.fullplots.2.lst,deeper.dendro.lst=deeper.dendro.lst,deeper.heat.lst=deeper.heat.lst))
		}
	}
