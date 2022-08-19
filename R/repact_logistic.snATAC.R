#' repact_logistic.snATAC
#'
#' This function is to run regression (linear/logistic) based on the number of PCs, and characteristics of samples.
#' @param OBJ, A snATAC Seurat (V3/V4), with LSIs and clustering info, the metadata should contain Sample column(it refers to donor in this study), characteristics to compare(it referes to if the donor/cells are from healthy or T2D, or it can be continuous, e.g. BMI)
#' @param LSIrange, Specified PCs for regression
#' @param pheno, the column name of the "characteristics to compare" in the phenodic.use dataframe
#' @param outputname, the prefix of output files
#' @return The function will write output files .pdf, .rds to disk
#' @import plot3D reshape2 ggplot2 Seurat Signac DESeq2 gridExtra dplyr plyr qvalue rlist
#' @export
#' @examples
#' repact_logistic.snATAC(snATAC.OBJ, LSIrange=1:20, pheno="diseaseStat", outputname="T2D_Beta.snATAC.RePACT")

repact_logistic.snATAC <- function(OBJ, LSIrange, pheno, outputname){
    require(plot3D)
    require(reshape2)
    require(ggplot2)
    require(Seurat)
    require(Signac)
    require(DESeq2)
    require(gridExtra)
    require(dplyr)
    require(plyr)
    require(qvalue)
    require(rlist)
    LSIrange <- paste("LSI_",LSIrange,sep='')
    LSI.withinfo <- merge(OBJ@reductions$lsi@cell.embeddings, OBJ@meta.data, by=0)
    rownames(LSI.withinfo) <- LSI.withinfo$Row.names
    LSI.withinfo <- LSI.withinfo[,-1]
    LSI.withinfo <- LSI.withinfo[!is.na(LSI.withinfo$Sample),]
    pheno.1 <- unique(as.character(LSI.withinfo[,pheno]))[1]
    pheno.2 <- unique(as.character(LSI.withinfo[,pheno]))[2]
    # p1 <- ggplot(melt(as.matrix(LSI.withinfo[,c(LSIrange, pheno)])))+aes(Var2,Var1,fill=value)+geom_boxplot(outlier.shape=NA)+theme_classic()+scale_fill_manual(values=c("dodgerblue","firebrick1"))+theme(axis.title=element_blank(),axis.text.x=element_text(angle=45,vjust=0.75))
    Select_LSI <- matrix(0,length(LSIrange),2)
    colnames(Select_LSI) <- c("FC","pval")
    rownames(Select_LSI) <- LSIrange
    for(i in 1:length(LSIrange)){
        Select_LSI[i,1] <- mean(LSI.withinfo[which(LSI.withinfo[, pheno]==pheno.1),i]) - mean(LSI.withinfo[which(LSI.withinfo[,pheno]==pheno.2),i])
        Select_LSI[i,2] <- wilcox.test(LSI.withinfo[which(LSI.withinfo[, pheno]==pheno.1),i], LSI.withinfo[which(LSI.withinfo[,pheno]==pheno.2),i])$p.value
    }
    Select_LSI.sig <- Select_LSI[which(Select_LSI[,"pval"]<0.0000000001),]
    Select_LSI.sig <- Select_LSI.sig[order(abs(Select_LSI.sig[,"FC"]),decreasing=T),]
    if(nrow(Select_LSI.sig)>=3){
        LSI_Top <- rownames(Select_LSI.sig)[1:3]
    }else{
        x <- Select_LSI[order(abs(Select_LSI[,"FC"]),decreasing=T),]
        LSI_Top <- rownames(x)[1:3]
    }
    p2 <- ggplot(LSI.withinfo)+aes(get(LSI_Top[1]),get(LSI_Top[2]),color=get(pheno))+geom_point(size=1.5)+theme_classic()+scale_color_manual(values=c("grey","red")) + xlab(LSI_Top[1]) + ylab(LSI_Top[2])
    p3 <- ggplot(LSI.withinfo)+aes(get(LSI_Top[1]),get(LSI_Top[3]),color=get(pheno))+geom_point(size=1.5)+theme_classic()+scale_color_manual(values=c("grey","red")) + xlab(LSI_Top[1]) + ylab(LSI_Top[3])
    p4 <- ggplot(LSI.withinfo)+aes(get(LSI_Top[2]),get(LSI_Top[3]),color=get(pheno))+geom_point(size=1.5)+theme_classic()+scale_color_manual(values=c("grey","red")) + xlab(LSI_Top[2]) + ylab(LSI_Top[3])
    pdf(paste(outputname,".CARePACT.pdf",sep=''),12,5)
#    print(p1)
    grid.arrange(p2,p3,p4,ncol=3)
    scatter3D(LSI.withinfo[,LSI_Top[1]], LSI.withinfo[,LSI_Top[2]],LSI.withinfo[,LSI_Top[3]],ticktype = "detailed", pch = 20, theta = 90, phi = 30, colvar = ifelse(LSI.withinfo[,pheno]==pheno.2,1,0), bty = "b2", cex = 0.3, col = alpha.col(col = c("steelblue", "red"), 0.6))
    LSI.withinfo.subs <- list()
    for(i in 1:100){
        LSI.withinfo.sub<-c()
        for (d in unique(LSI.withinfo$Sample)){
        tmp<-subset(LSI.withinfo,Sample==d)
        if(nrow(tmp)<=500){
            LSI.withinfo.sub <- rbind(LSI.withinfo.sub,tmp)
        }else{
            LSI.withinfo.sub <- rbind(LSI.withinfo.sub,tmp[sample(1:nrow(tmp),500),])
        }
        }
        LSI.withinfo.subs <- c(LSI.withinfo.subs,list(LSI.withinfo.sub))
    }
    pseudo.indexes<-list()
    for(i in 1:100){
        md <- GetRePACTmodel.cca(ccaWithinfo=LSI.withinfo.subs[[i]],prefix="LSI",pheno=pheno,CCrange=c(1:10))
        trainingdata <- LSI.withinfo.subs[[i]]
        Restdata <- LSI.withinfo[setdiff(row.names(LSI.withinfo),row.names(LSI.withinfo.subs[[i]])),]
        Alldata <- LSI.withinfo
        LSI.withinfo.subs[[i]]$pseudo.index<-apply(trainingdata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        Restdata$pseudo.index<-apply(Restdata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        pseudo.indexes<-c(pseudo.indexes,list(apply(Alldata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})))
    }
    md.all<-GetRePACTmodel.cca(ccaWithinfo=LSI.withinfo,prefix="LSI",pheno=pheno,CCrange=c(1:10))
    LSI.withinfo$pseudo.index = md.all$linear.predictors
    LSI.withinfo$pseudo.index.balanced<-do.call(cbind,pseudo.indexes) %>% rowMeans()
    ggplot(LSI.withinfo)+aes(Sample,pseudo.index.balanced,fill=get(pheno))+geom_violin()+geom_boxplot(width=0.2,outlier.shape = NA,notch=F,coef = 0,fill="grey25",color="grey75")+coord_flip()+theme_classic()+scale_fill_manual(values=c("steelblue", "red"))+theme_bw()+theme(legend.position="none") 
    ggplot(LSI.withinfo)+aes(pseudo.index.balanced,fill=get(pheno))+geom_density()+scale_fill_manual(values=c("steelblue","red"))+theme_classic()+theme(legend.position="none")
    LSI.withinfo$rank<-rank(LSI.withinfo$pseudo.index.balanced)
    tmp.ATAC.LSI.20bin.ob <- MakeEvenBinBydepth(cellvsPeak.m=t(OBJ@assays$ATAC@counts), data.info=LSI.withinfo, binnumber=20)
    tmpT2D.diffPeaks.20bin.LSI<-CallT2Dpeak_qvalue(tmp.ATAC.LSI.20bin.ob$cellvsPeak.m.aggr,tmp.ATAC.LSI.20bin.ob$depths,tmp.ATAC.LSI.20bin.ob$index,qcut=0.01,slopecut1=0.5,slopecut2=-0.5,doscale=T)
    Evenbin.donorContribute<-tmp.ATAC.LSI.20bin.ob$data.info.withbin %>% .[,c("Sample","evenfragbin")] %>% table %>% as.matrix %>% apply(.,1,function(x){x/sum(x)}) %>% melt
    Evenbin.donorContribute<-cbind(Evenbin.donorContribute,disease=ifelse(grepl(pheno.2,Evenbin.donorContribute$Sample),pheno.2,pheno.1))
    p <- ggplot(Evenbin.donorContribute)+aes(Sample,value,fill=disease)+geom_bar(stat="identity",color="black")+facet_grid(~evenfragbin)+theme(axis.text=element_blank(),axis.ticks=element_blank())+scale_fill_manual(values=c("steelblue","red"))+theme_bw()
    print(p)
    apply(tmpT2D.diffPeaks.20bin.LSI$UPDN.toplot[,c(rownames(tmpT2D.diffPeaks.20bin.LSI$UP), rownames(tmpT2D.diffPeaks.20bin.LSI$DN))],2,function(x){scale(x)}) %>% melt() %>% ggplot()+aes(Var1,Var2,fill=value)+geom_tile()+scale_fill_gradient2(low="steelblue",mid="white",high="red")+theme_classic()+theme(axis.text=element_blank())+ggtitle("LSI-20bins")
    dev.off()
    saveRDS(tmpT2D.diffPeaks.20bin.LSI, file=paste(outputname,".CARePACT.T2D_diffPeaks.20bin.LSI.rds",sep=''))
    ALL.tmp.ATAC.bulk.lis <- list()
    sampleLis <- unique(as.character(OBJ@meta.data$Sample))
    for(sampleInd in sampleLis){
            ALL.tmp.ATAC.bulk.lis[[sampleInd]] <- rowSums(OBJ@assays$ATAC@counts[,row.names(subset(OBJ@meta.data, Sample==sampleInd))])
    }
    ALL.tmp.ATAC.bulk <- list.cbind(ALL.tmp.ATAC.bulk.lis)
    Bulk.meta <- unique(OBJ@meta.data[,c("Sample", pheno)])
    # colnames(Bulk.meta) <- c("row.names", "condition")
    rownames(Bulk.meta) <- NULL
    dds <- DESeqDataSetFromMatrix(countData = ALL.tmp.ATAC.bulk, colData = Bulk.meta, design = ~ diseaseStat)
    dds <- DESeq(dds)
    res<-results(dds)
    summary(res)
    tmpT2D.diffPeaks.repact.bulk <- merge(tmpT2D.diffPeaks.20bin.LSI$pseudoregress.all[,c("slope","qvalue")],as.data.frame(results(dds))[,c("log2FoldChange","pvalue","padj")],by=0)
    rownames(tmpT2D.diffPeaks.repact.bulk) <- tmpT2D.diffPeaks.repact.bulk$Row.names
    tmpT2D.diffPeaks.repact.bulk <- tmpT2D.diffPeaks.repact.bulk[,-1]
    caRePACTpeaks<-c()
    caRePACTpeaks[which(row.names(tmpT2D.diffPeaks.repact.bulk) %in% row.names(tmpT2D.diffPeaks.20bin.LSI$UP))]<-"UP"
    caRePACTpeaks[which(row.names(tmpT2D.diffPeaks.repact.bulk) %in% row.names(tmpT2D.diffPeaks.20bin.LSI$DN))]<-"DN"
    caRePACTpeaks[which(!row.names(tmpT2D.diffPeaks.repact.bulk) %in% c(row.names(tmpT2D.diffPeaks.20bin.LSI$UP),row.names(tmpT2D.diffPeaks.20bin.LSI$DN)))]<-""
    tmpT2D.diffPeaks.repact.bulk$caRePACTpeaks<-caRePACTpeaks
    deseq0.01<-c()
    deseq0.01[which(tmpT2D.diffPeaks.repact.bulk$pvalue<0.01 &tmpT2D.diffPeaks.repact.bulk$log2FoldChange>0)]<-"bulk.up"
    deseq0.01[which(tmpT2D.diffPeaks.repact.bulk$pvalue<0.01 &tmpT2D.diffPeaks.repact.bulk$log2FoldChange<0)]<-"bulk.dn"
    deseq0.01[which(tmpT2D.diffPeaks.repact.bulk$pvalue>0.01)]<-""
    deseq0.05<-c()
    deseq0.05[which(tmpT2D.diffPeaks.repact.bulk$pvalue<0.05 &tmpT2D.diffPeaks.repact.bulk$log2FoldChange>0)]<-"bulk.up"
    deseq0.05[which(tmpT2D.diffPeaks.repact.bulk$pvalue<0.05 &tmpT2D.diffPeaks.repact.bulk$log2FoldChange<0)]<-"bulk.dn"
    deseq0.05[which(tmpT2D.diffPeaks.repact.bulk$pvalue>0.05)]<-""
    tmpT2D.diffPeaks.repact.bulk$deseq0.01<-deseq0.01
    tmpT2D.diffPeaks.repact.bulk$deseq0.05<-deseq0.05
    saveRDS(tmpT2D.diffPeaks.repact.bulk,paste(outputname,".T2D.diffPeaks.repact.bulk.rds",sep=''))
}

