

InSample.compare.nb<-function(OBJ, pseudo=Islet12.BetaPCA.pseudo[,"pseudo.index.balanced",drop=F],gene="INS",GetDatatoplot=F){
        require(MASS)
        require(dplyr)
        require(lmtest)
        require(EZsinglecell)
	matrix <- OBJ@assays$RNA@counts
	meta <- OBJ@meta.data
	sample.col <- "Sample"

        ps<-c()
        datatoplot<-c()
        for(sample in levels(meta[,sample.col])){
                 meta.tmp<-Tomerge_v2(subset(meta, Sample==sample),pseudo) %>% .[complete.cases(.),]
                 quantiles<-quantile(meta.tmp[,"pseudo.index.balanced"]) %>% as.numeric
                 meta.tmp$Q<-sapply(meta.tmp[,"pseudo.index.balanced"],AddQ,quantiles=quantiles)
                 meta.tmp<-cbind(meta.tmp,gene=matrix[gene,row.names(meta.tmp)])
                 Ave.expr<-ddply(meta.tmp,.(Q),summarize,sum(nCount_RNA),sum(gene)) %>% data.frame(row.names=.$Q,gene=1e3*(.[,3]/.[,2]))
                 meta.tmp$Ave.expr<-Ave.expr[as.character(meta.tmp$Q),"gene"]
                 meta.tmp.Q1Q4<-subset(meta.tmp,Q %in% c("Q1","Q4"))
                 md<-try(glm.nb(gene ~ Q + log10(nCount_RNA),data=meta.tmp.Q1Q4),silent=T)
                 md0<-try(glm.nb(gene ~ log10(nCount_RNA),data=meta.tmp.Q1Q4),silent=T)
                 if(all(grepl("Error",md))){
                         p<-NA
                 }else{
                         p<-lrtest(md,md0)[[5]][2]
                 }
                 ps<-c(ps,p)
                 datatoplot<-rbind(datatoplot,meta.tmp)
        }
        FCs<-datatoplot[,c("Sample","Q","Ave.expr")] %>% unique %>% subset(.,Q %in% c("Q1","Q4")) %>% .[order(.$Sample,.$Q),] %>% ddply(.,.(Sample),summarize,Ave.expr[2]/Ave.expr[1]) %>% .[,2]
        names(ps)<-levels(meta[,sample.col])
        names(FCs)<-levels(meta[,sample.col])
        if(GetDatatoplot){
                return(list(ps=ps,FCs=FCs,datatoplot=datatoplot))
        }else{
                return(list(ps=ps,FCs=FCs))
        }
}

