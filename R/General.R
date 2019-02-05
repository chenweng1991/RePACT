#' format2
#'
#' This function is to make a ggplot figure into a square and clean
#' @param p  a ggplot plot
#' @param data the data ggplot took for plotting
#' @param x the column name of the x
#' @param y the column name of the y
#' @param nolegend if true then take out the legend
#' @return this will return the clean and square plot
#' @export
#' @examples
#' Gettopgenes(XX.ob,number)

format2<-function(p,data,x,y,nolegend=T){
require(ggplot2)
data<-data[complete.cases(data),]
p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+coord_fixed( (max(data[,x])-min(data[,x]))/(max(data[,y])-min(data[,y])))
if(nolegend)
{
p<-p+theme(legend.position="none")
}
return(p)
}



#' GetPCAcelldata_v2
#'
#' This function is to quickly extract data info together with  PCA
#' @param object  The seurat object that has been analyzed.
#' @return  return a data frame with both sample/cell information and PCA information
#' @export
#' @examples
#' Gettopgenes(XX.ob,number)

GetPCAcelldata_v2<-function(object){
Sample.info<-object@data.info[complete.cases(object@data.info),]
cellPCAstate<-Tomerge(object@pca.rot,Sample.info)
cellPCAstate<-cellPCAstate[complete.cases(cellPCAstate),]
return(cellPCAstate)
}

#' Tomerge
#'
#' This function is to quickly merge two dataframe by rownames, but only leave the rows both dataframe have
#' @param A  dataframe A
#' @param B  dataframe B
#' @return  return a data frame with merged information
#' @export
#' @examples
#' Gettopgenes(XX.ob,number)

Tomerge<-function(A,B){
mergeAB<-merge(A,B,by="row.names",all=TRUE)
row.names(mergeAB)<-mergeAB[,1]
mergeAB<-mergeAB[,-1]
return(mergeAB)
}



#' Tomerge_v2
#'
#' This function is to quickly merge two dataframe by rownames, but can choose to leave A or B all information
#' @param A  dataframe A
#' @param B  dataframe B
#' @return  return a data frame with merged information
#' @export
#' @examples
#' Gettopgenes(XX.ob,number)

Tomerge_v2<-function(A,B,leavex=T,leavey=F){
	mergeAB<-merge(A,B,by="row.names",all.x=leavex,all.y=leavey)
	row.names(mergeAB)<-mergeAB[,1]
	mergeAB<-mergeAB[,-1]
	return(mergeAB)
}
