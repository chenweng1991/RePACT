#' datapair.mk
#'
#' This function is to make a datapair that includes UMI/normalized UMI ,atrix as well as information table, out of one or multiple seurat objects
#' @param Seurat.ob.list a list of seurat object that I will extract data from
#' @param cols a vector of columns that will be considered for each object in order
#' @param pick.list a list of vectors that set the rules to take data
#' @param normalizecellsize, if true, then normalize the matrix by total cell umi. Default is F
#' '@param randomizecelloirder If T, randomize the cell order in matrix, default is T
#' @return  return data pair that can be used for DE, bubble plot and others.
#' @export
#' @examples
#' ROCKvsnorock.endo.paired<-datapair.mk(list(S7rock=S7rock_1.ob,S7=S7.ob),cols=c("Sample","Sample.2nd"),pick.list=list(c("s7.RockII_1"),c("s7.B_1")),normalizecellsize=F,randomizecelloirder=T)

datapair.mk<-function(Seurat.ob.list,cols,pick.list,normalizecellsize=F,randomizecelloirder=T)
{
		data.collect<-c()
		info.collect<-c()
		info.first<-Seurat.ob.list[[1]]@data.info[Seurat.ob.list[[1]]@data.info[,cols[1]] %in% pick.list[[1]],]
		picked.cells.names<-row.names(Seurat.ob.list[[1]]@data.info)[Seurat.ob.list[[1]]@data.info[,cols[1]] %in% pick.list[[1]]]
		if(randomizecelloirder)
		{
			picked.cells.names<-sample(picked.cells.names)
		}
		if(normalizecellsize)
		{
			data.first<-as.matrix(Seurat.ob.list[[1]]@raw.data[,picked.cells.names])
			size<-colSums(data.first)
			data.first<-t(apply(data.first,1,function(x){x/size*10^4}))
		}else
		{
			data.first<-as.matrix(Seurat.ob.list[[1]]@raw.data[,picked.cells.names])
		}
		info.first<-info.first[picked.cells.names,]
		row.names(info.first)<-paste(picked.cells.names,1,sep="_")
		colnames(data.first)<-paste(picked.cells.names,1,sep="_")
		info.first<-info.first[,c("nGene","nUMI",cols[1])]
		info.first[,1]<-as.numeric(as.character(info.first[,1]))
		info.first[,2]<-as.numeric(as.character(info.first[,2]))
		colnames(info.first)[3]<-"name"
		data.collect<-data.first
		info.collect<-info.first
		if(length(Seurat.ob.list)>1)
		{
			for (i in 2: length(Seurat.ob.list))
			{
				info.current<-Seurat.ob.list[[i]]@data.info[Seurat.ob.list[[i]]@data.info[,cols[i]] %in% pick.list[[i]],]
				picked.cells.names<-row.names(Seurat.ob.list[[i]]@data.info)[Seurat.ob.list[[i]]@data.info[,cols[i]] %in% pick.list[[i]]]
				if(randomizecelloirder)
				{
					picked.cells.names<-sample(picked.cells.names)
				}
				if(normalizecellsize)
				{
					data.current<-as.matrix(Seurat.ob.list[[i]]@raw.data[,picked.cells.names])
					size<-colSums(data.current)
					data.current<-t(apply(data.current,1,function(x){x/size*10^4}))
				}else
				{
					data.current<-as.matrix(Seurat.ob.list[[i]]@raw.data[,picked.cells.names])
				}
				info.current<-info.current[picked.cells.names,]
				row.names(info.current)<-paste(picked.cells.names,i,sep="_")
				colnames(data.current)<-paste(picked.cells.names,i,sep="_")
				info.current<-info.current[,c("nGene","nUMI",cols[i])]
				info.current[,1]<-as.numeric(as.character(info.current[,1]))
				info.current[,2]<-as.numeric(as.character(info.current[,2]))
				colnames(info.current)[3]<-"name"
				data.collect<-Tomerge_v2(data.collect,data.current)
				data.collect[is.na(data.collect)]<-0
				info.collect<-rbind(info.collect,info.current)
			}
		}
		info.collect$name<-factor(info.collect$name,levels=unlist(pick.list))
		return(list(data=data.collect,info=info.collect))
}

#' DE.gettripple
#'
#' This function is to prepare the data format that is used to differentially expression calling. It include the raw matrix; data.info and size effect
#' @param datapair    tyhe datapair generated from datapair.mk
#' @param cpcol The column name for comparison.
#' @param withscran  if true, use deconvolution to calculate size effect.
#' @return  This will return .tri.dummy file that is the input for DE analysis
#' @export
#' @examples
#' ROCKvsnorock.endo.tri.dummy<-DE.gettripple(ROCKvsnorock.endo.paired,cpcol="name")
DE.gettripple<-function(datapair,cpcol,withscran=F)
{
#	library(scran)
	x<-datapair$data
	x[is.na(x)]<-0
	x.info<-datapair$info
	if(!all(colnames(x) ==row.names(x.info)))
	{
		stop("The datapair is not match well")
	}else
	{
		print("datapair match well......continue")
	}
#
	wherezero<-function(vector)
	{

		if (all(vector==0))
		{
			return(FALSE)
		}
		else
		{
			return(TRUE)
		}
	}
#
	makedummycell<-function(vector)
	{
		dummycell<-0
		if(all(vector==0))
		{
			dummycell<-1
		}
		return(dummycell)
	}
#
	x<-x[,apply(x,2,wherezero)]  # get rid of cells with no expression
	if(withscran)
	{
		sf<-computeSumFactors(as.matrix(x), positive=F)  # calculate  size factor
		x<-x[,which(sf>0)]
		x.info<-x.info[which(sf>0),]
		sf<-sf[which(sf>0)]
	}else
	{
		sf<-colSums(x)
	}
	for (celltype in unique(x.info[,cpcol]) )
	{
		subdata<-x[,row.names(x.info)[which(x.info[,cpcol]==celltype)]]
		subdata[is.na(subdata)]<-0
		x<-cbind(x,apply(subdata,1,makedummycell))
		colnames(x)[ncol(x)]<-celltype
		newinfoline<-x.info[1,]
		newinfoline[,cpcol]<-celltype
		x.info<-rbind(x.info,newinfoline)
		row.names(x.info)[nrow(x.info)]<-celltype
		sf<-c(sf,median(sf))
		names(sf)[length(sf)]<-celltype
	}
	return(list(data=x,info=x.info,sf=sf))
}


#' DoDE
#'
#' This is the main function for calculating differentially expressed genes
#' @param tri.dummy this is generated from  DE.gettripple
#' @param cpcol  the column in tri.dummy$info, the contents of which are used for iteratively compare with one another
#' @param onlyoneSample If true, regress out batch effect. Notice, there should be a "Sample" column in in tri.dummy$info that indicate sample or donor or batch
#' @param cpus a number of cpus being used for calculation, default is 16
#' @return  return a list that includes all DE result iteratively
#' @export
#' @examples
#' ROCKvsnorock.endo.de<-DoDE(ROCKvsnorock.endo.tri.dummy,"name",onlyoneSample=T,cpus=16)

DoDE<-function(tri.dummy,cpcol,onlyoneSample=F,cpus=16)
{

Donbregression<-function(data,info,cpcol,reftype,sf,gene,onlyoneSample){
		require(MASS)
		info[,cpcol]<-factor(info[,cpcol],levels=c(reftype,setdiff(levels(info[,cpcol]),reftype)))
		if(onlyoneSample)
		{
			md<-try(glm.nb(data[gene,] ~ info[,cpcol]+log(sf)),silent=T)
		}else
		{
			md<-try(glm.nb(data[gene,] ~ info[,cpcol]+info[,"Sample"]+log(sf)),silent=T)
		}
					coeff.table<-try(summary(md)$coefficients,silent=T)
					N<-2
					output<-c()
					for(i in setdiff(levels(info[,cpcol]),reftype))
					{
							p.value<-try(coeff.table[N,4],silent=T)
							if(grepl("Error",p.value))
							{
									p.value<-NA
							}
				ref.mean<-mean(data[,row.names(info)[which(info[,cpcol]==reftype)]])
				i.mean<-mean(data[,row.names(info)[which(info[,cpcol]==i)]])
							assign(gene,data.frame(
							ref.sum=as.integer(sum(data[gene,row.names(info)[which(info[,cpcol]==reftype)]])),
							alt.sum=as.integer(sum(data[gene,row.names(info)[which(info[,cpcol]==i)]])),
							ref.nonzero=as.integer(length(which(data[gene,row.names(info)[which(info[,cpcol]==reftype)]]>0))),
							alt.nonzero=as.integer(length(which(data[gene,row.names(info)[which(info[,cpcol]==i)]]>0))),
							nonzero.ratio=(length(which(data[gene,row.names(info)[which(info[,cpcol]==reftype)]]>0))/as.integer(as.matrix(table(info[,cpcol]))[reftype,]))/(length(which(data[gene,row.names(info)[which(info   [,cpcol]==i)]]>0))/as.integer(as.matrix(table(info[,cpcol]))[i,])),
							mean.ratio=(mean(data[gene,row.names(info)[which(info[,cpcol]==reftype)]])/ref.mean)/(mean(data[gene,row.names(info)[which(info[,cpcol]==i)]])/i.mean),
							p.value=p.value
							))
							cline<-get(gene)
							row.names(cline)<-gene
							clist<-list(cline)
							names(clist)<-i
							output<-c(output,clist)
							N<-N+1
					}
					return(output)
	}

calculateDE.sg<-function(data,info,cpcol,sf,gene,onlyoneSample){
	singlegene.result<-list()
	length(singlegene.result)<-length(levels(info[,cpcol]))
	names(singlegene.result)<-levels(info[,cpcol])
	for (reftype in levels(info[,cpcol])){
		singlegene.result[[reftype]]<-Donbregression(data,info,cpcol,reftype,sf,gene,onlyoneSample)
		}
		return(singlegene.result)
}

rbindlist<-function(pardoresult,info,cpcol){
		print("Start binding......")
		DEresult<-list()
		length(DEresult)<-length(levels(info[,cpcol]))
		names(DEresult)<-levels(info[,cpcol])
		for (reftype in levels(info[,cpcol]))
		{
			DEresult[[reftype]]<-list(DEresult[[reftype]])
			length(DEresult[[reftype]])<-length(levels(info[,cpcol]))-1
			names(DEresult[[reftype]])<-setdiff(levels(info[,cpcol]),reftype)
		}
		n<-0
		for(reftype in levels(info[,cpcol]))
		{
			print(n)
			n<-n+1
	    		for (alttype in setdiff(levels(info[,cpcol]),reftype))
	    		{
	        		for (i in 1:length(pardoresult))
	        		{
	            			DEresult[[reftype]][[alttype]]<-rbind(DEresult[[reftype]][[alttype]],pardoresult[[i]][[reftype]][[alttype]])
	        		}
	    		}
		cellnumber<-list(as.matrix(table(info[,cpcol]))[reftype,])
		names(cellnumber)="cellnumber"
		DEresult[[reftype]]<-c(DEresult[[reftype]],cellnumber)
		}
		return(DEresult)
	}

# main function starting from here
data<-tri.dummy[[1]]
info<-tri.dummy[[2]]
sf<-tri.dummy[[3]]
data<-as.matrix(data)
if(!all(colnames(data) ==row.names(info)))
{
	stop("datapair doesnt match well....data:info doesn't match")
}else
{
	print("data::info match well......continue")
}
if(!all(colnames(data) ==row.names(sf)))
{
	print("datapair doesnt match well....data:sf doesn't match")
}else
{
	print("data::sf match well......continue")
}
print(paste("The compared column is",cpcol))
info[,cpcol]<-factor(info[,cpcol],levels=as.character(unique(info[,cpcol])))
print(table(info[,cpcol]))
timestart<-Sys.time()
require(doParallel)
cl <- makeCluster(cpus)
doParallel::registerDoParallel(cl)
time.start<-Sys.time()
for (reftype in levels(info[,cpcol]))
{
	assign(reftype,list())
}
pardoresult<-list()
pardoresult<-foreach (gene=row.names(data)) %dopar%
{
	print(which(row.names(data)==gene))
	print(calculateDE.sg(data,info,cpcol,sf,gene,onlyoneSample))
}
print("pardoresult Get!!")
time.ends<-Sys.time()
DEresult<-rbindlist(pardoresult,info,cpcol)
timeend<-Sys.time()
print(timeend-timestart)
return(DEresult)
}
