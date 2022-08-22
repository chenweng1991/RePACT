#' Do_heatmap
#'
#' This function is to generate the heatmap for variable genes across HT->T2D 20 trajactory bins
#' @param bindata, "RepACT.2nd.ob$bin.data" in the repact.scRNA output list
#' @param df1, "RepACT.2nd.ob$BINlinear.result.summarized$UP" in the repact.scRNA output list
#' @param df2, "RepACT.2nd.ob$BINlinear.result.summarized$DOWN" in the repact.scRNA output list
#' @param rankname
#' @param top_gene_num, the number of differentially genes across peusoindex for heatmap visualization
#' @return the output files will be written to disk
#' @import ggplot2 Matrix RColorBrewer gridExtra pscl
#' @export
#' @examples
#' Do_heatmap(T2D.scRNA.RePACT$RepACT.2nd.ob$bin.data,df1=T2D.scRNA.RePACT$RepACT.2nd.ob$BINlinear.result.summarized$UP,df2=T2D.scRNA.RePACT$RepACT.2nd.ob$BINlinear.result.summarized$DOWN,rankname="rank",top_gene_num=20)

Do_heatmap <- function(bindata,df1,df2,rankname,top_gene_num=10,hardadd=NULL,last=T,insertinto=0,doreturn=T){
    normalize_01<-function(vector){
            normed<-(vector-min(vector))/(max(vector)-min(vector))
            return(normed)
    }
    require(ggplot2)
    require(gridExtra)
    top_gene_num_DOWN=top_gene_num
    top_gene_num_UP=top_gene_num
    if(nrow(df1)<top_gene_num){
            top_gene_num_UP=nrow(df1)
    }
    if(nrow(df2)<top_gene_num){
            top_gene_num_DOWN=nrow(df2)
    }
    DEgenelist<-c(row.names(df1)[1:top_gene_num_UP],row.names(df2)[1:top_gene_num_DOWN])
    if(last==T){
            DEgenelist.L<-DEgenelist[1:(length(DEgenelist)-insertinto)]
            DEgenelist.R<-setdiff(DEgenelist,DEgenelist.L)
            DEgenelist<-c(DEgenelist.L,hardadd,DEgenelist.R)
    }else{
            DEgenelist.R<-DEgenelist[(insertinto+1):length(DEgenelist)]
            DEgenelist.L<-setdiff(DEgenelist,DEgenelist.R)
            DEgenelist<-c(DEgenelist.L,hardadd,DEgenelist.R)
    }
    DEgenelist<-na.omit(DEgenelist)
    bindata<-bindata[,c(DEgenelist,"tag")]
    bindata<-data.frame(apply(bindata[,-ncol(bindata),drop=F],2,normalize_01),bin=bindata[,ncol(bindata)])
    bindata.m<-reshape2::melt(bindata,id.vars="bin")
    p<-ggplot(bindata.m)+aes(bin,variable,fill=value)+geom_tile()+scale_fill_gradient2(low="white",high="red",mid="orange",midpoint=0.6)+labs(y="Variable genes")
    return(p)
}

