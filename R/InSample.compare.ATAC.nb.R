#' InSample.compare.ATAC.nb
#'
#' This function is to identify intra-donor genes by negative binomial
#' @param matrix, cellxpeak matrix
#' @param meta
#' @param pseudo
#' @param sample.col
#' @return The function return logistic model
#' @import Seurat plyr dplyr Seurat ggrepel ggplot2 Matrix 
#' @export
#' @examples
#' donorWise.T2D.RePACT <- DonorWise.scRNA.RePACT(OBJ=scRNA.OBJ, pheno="diseaseStat", outputname="T2D.RePACT.donorWise")

InSample.compare.ATAC.nb <- function(matrix=ALL.Beta.ATAC.Mtx,meta=Beta.all.metadata,pseudo=beta.atac.LSI.withinfo[,"pseudo.index.balanced",drop=F],sample.col="pubdonor",peak="chr11_92949418_92950719",GetDatatoplot=T){
        require(MASS)
        require(Matrix)
        require(dplyr)
        require(lmtest)
        ps<-c()
        datatoplot<-c()
        for(sample in levels(meta[,sample.col])){
                 print(sample)
                 meta.tmp<-Tomerge_v2(subset(meta,pubdonor==sample),pseudo) %>% .[complete.cases(.),]
                 quantiles<-quantile(meta.tmp[,"pseudo.index.balanced"]) %>% as.numeric
                 meta.tmp$Q<-sapply(meta.tmp[,"pseudo.index.balanced"],AddQ,quantiles=quantiles)
                 meta.tmp<-cbind(meta.tmp,peak=matrix[peak,row.names(meta.tmp)])
                 Ave.intensity<-ddply(meta.tmp,.(Q),summarize,sum(Fragments),sum(peak)) %>% data.frame(row.names=.$Q,peak=1e3*(.[,3]/.[,2]))
                 meta.tmp$Ave.intensity<-Ave.intensity[as.character(meta.tmp$Q),"peak"]
                 meta.tmp.Q1Q4<-subset(meta.tmp,Q %in% c("Q1","Q4"))
                 md<-try(glm.nb(peak ~ Q + log10(Fragments),data=meta.tmp.Q1Q4),silent=T)
                 md0<-try(glm.nb(peak ~ log10(Fragments),data=meta.tmp.Q1Q4),silent=T)
                 if(all(grepl("Error",md))|all(grepl("Error",md0))){
                         p<-NA
                 }else{
                         p<-lrtest(md,md0)[[5]][2]
                 }
                 ps<-c(ps,p)
                 datatoplot<-rbind(datatoplot,meta.tmp)
        }
        FCs<-datatoplot[,c("pubdonor","Q","Ave.intensity")] %>% unique %>% subset(.,Q %in% c("Q1","Q4")) %>% .[order(.$pubdonor,.$Q),] %>% ddply(.,.(pubdonor),summarize,Ave.intensity[2]
/Ave.intensity[1]) %>% .[,2]
        names(ps)<-levels(meta[,sample.col])
        names(FCs)<-levels(meta[,sample.col])
        if(GetDatatoplot){
                return(list(ps=ps,FCs=FCs,datatoplot=datatoplot))
        }else{
                return(list(ps=ps,FCs=FCs))
        }
}

