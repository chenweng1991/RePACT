#' DonorWise.scRNA.RePACT
#'
#' This function is to run logistic regression based on the number of LSIs, and characteristics of samples.
#' @param OBJ, a scRNA-seq Seurat object
#' @param PCrange, Specified PCs for regression
#' @param pheno, the column name of the "characteristics to compare" in the phenodic.use dataframe
#' @return The function return logistic model
#' @import Seurat plyr dplyr Seurat ggrepel ggplot2
#' @export
#' @examples
#' donorWise.T2D.RePACT <- DonorWise.scRNA.RePACT(OBJ=scRNA.OBJ, pheno="diseaseStat", outputname="T2D.RePACT.donorWise")

DonorWise.scRNA.RePACT <- function(OBJ, pheno){
    require(plyr)
    require(dplyr)
    require(Seurat)
    require(ggrepel)
    PCAInfo <- OBJ@reductions$pca@cell.embeddings %>% apply(.,2,function(x){x/sqrt(sum(x^2))}) %>% merge(., OBJ@meta.data, by=0)
    rownames(PCAInfo) <- PCAInfo$Row.names
    PCAInfo <- PCAInfo[,-1]
    beta.rna.pca.withinfo.subs<-list()
    for(i in 1:100){
        beta.rna.pca.withinfo.sub<-c()
        for (d in levels(PCAInfo$Sample)){
        tmp<-subset(PCAInfo,Sample==d)
        if(nrow(tmp)<=200){
            beta.rna.pca.withinfo.sub<-rbind(beta.rna.pca.withinfo.sub,tmp)
        }else{
            beta.rna.pca.withinfo.sub<-rbind(beta.rna.pca.withinfo.sub,tmp[sample(1:nrow(tmp),200),])
        }
        }
        beta.rna.pca.withinfo.subs<-c(beta.rna.pca.withinfo.subs,list(beta.rna.pca.withinfo.sub))
    }
    pseudo.indexes<-list()
    for(i in 1:100){
        md<-GetRePACTmodel.cca(ccaWithinfo=beta.rna.pca.withinfo.subs[[i]],prefix="PC",pheno=pheno,CCrange=c(1:10))
        trainingdata<-beta.rna.pca.withinfo.subs[[i]]
        Restdata<-PCAInfo[setdiff(row.names(PCAInfo),row.names(beta.rna.pca.withinfo.subs[[i]])),]
        Alldata<-PCAInfo
        beta.rna.pca.withinfo.subs[[i]]$pseudo.index<-apply(trainingdata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        Restdata$pseudo.index<-apply(Restdata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        pseudo.indexes<-c(pseudo.indexes,list(apply(Alldata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})))
    }
    PCAInfo$pseudo.index.balanced <- do.call(cbind,pseudo.indexes) %>% rowMeans()
    PCAInfo$rank<-rank(PCAInfo$pseudo.index.balanced)
    PCAInfo.20bin.ob <- MakeEvenBinBydepth(cellvsPeak.m=t(as.matrix(OBJ@assays$RNA@counts)),data.info=PCAInfo[,51:ncol(PCAInfo)],binnumber=20)
    PCAInfo.20bin.ob.PCA <- CallT2Dpeak_qvalue(PCAInfo.20bin.ob$cellvsPeak.m.aggr, PCAInfo.20bin.ob$depths, PCAInfo.20bin.ob$index,qcut=0.2,slopecut1=0.3,slopecut2=-0.3,doscale=T)   
    bin20.g.up <- row.names(subset(PCAInfo.20bin.ob.PCA$UP,qvalue<0.005 & MaxRPKM>quantile(PCAInfo.20bin.ob.PCA$pseudoregress.all$MaxRPKM,0.25)))
    bin20.g.dn <- row.names(subset(PCAInfo.20bin.ob.PCA$DN,qvalue<0.005 & MaxRPKM>quantile(PCAInfo.20bin.ob.PCA$pseudoregress.all$MaxRPKM,0.25)))
    RePACT.diff.genes.bydonor<-list()
    for (curDonor in levels(PCAInfo$Sample)){
        print(curDonor)
        PCAInfo.donor <- subset(PCAInfo, Sample==curDonor)
        PCAInfo.donor$rank <- rank(PCAInfo.donor$pseudo.index.balanced)
        PCAInfo.20bin.donor.ob <- MakeEvenBinBydepth(cellvsPeak.m=t(as.matrix(OBJ@assays$RNA@counts[,row.names(PCAInfo.donor)])),data.info=PCAInfo.donor[,51:ncol(PCAInfo.donor)],binnumber=20)
        PCAInfo.20bin.donor.PCA <- CallT2Dpeak_pvalueOneTail(PCAInfo.20bin.donor.ob$cellvsPeak.m.aggr, PCAInfo.20bin.donor.ob$depths, PCAInfo.20bin.donor.ob$index,doscale=T, GlobalSlopes=PCAInfo.20bin.ob.PCA$pseudoregress.all[,1,drop=F])
        RePACT.diff.genes.bydonor <- c(RePACT.diff.genes.bydonor,list(PCAInfo.20bin.donor.PCA))
    }
    names(RePACT.diff.genes.bydonor)<-levels(PCAInfo$Sample)
    ps<-lapply(RePACT.diff.genes.bydonor,function(x){x$pvalue.onetail}) %>% do.call(cbind,.)
    row.names(ps)<-row.names(RePACT.diff.genes.bydonor[[1]])
    ps[is.na(ps)]<-1
    slopes<-lapply(RePACT.diff.genes.bydonor,function(x){x$pseudoregress.all$slope}) %>% do.call(cbind,.)
    row.names(slopes)<-row.names(RePACT.diff.genes.bydonor[[1]]$pseudoregress.all)
    FishersMethod.p<-apply(ps,1,function(x){-2*sum(log(x))}) %>% pchisq(.,2*11,lower.tail=FALSE)  ## k=11, df=2k
    FishersMethod.p<-FishersMethod.p[c(bin20.g.up,bin20.g.dn)]
    FishersMethod.q<-qvalue(FishersMethod.p)$qvalues %>% .[order(.)]
    FishersMethod.q.df<-data.frame(gene=names(FishersMethod.q),qvalueInSample=FishersMethod.q)
    Globaltag<-c()
    Globaltag[which(FishersMethod.q.df$gene %in% bin20.g.up)]<-"UP"
    Globaltag[which(FishersMethod.q.df$gene %in% bin20.g.dn)]<-"DN"
    Globaltag[which(!FishersMethod.q.df$gene %in% c(bin20.g.up,bin20.g.dn))]<-""
    FishersMethod.q.df$Globaltag<-Globaltag
    GeneLabel<-c(subset(FishersMethod.q.df,Globaltag=="UP") %>% .[order(.$qvalueInSample),] %>% head(.,n=25) %>% row.names,subset(FishersMethod.q.df,Globaltag=="DN") %>% .[order(.$qvalueInSample),] %>% head(.,n=25) %>% row.names)
    FishersMethod.q.df$GeneLabel<-ifelse(FishersMethod.q.df$gene %in% GeneLabel,row.names(FishersMethod.q.df),"")
    FishersMethod.q.df<-Tomerge_v2(FishersMethod.q.df,PCAInfo.20bin.ob.PCA$pseudoregress.all[,"qvalue",drop=F])
    FishersMethod.q.df$Globaltag<-as.factor(FishersMethod.q.df$Globaltag)
    FishersMethod.q.dn.df<-subset(FishersMethod.q.df,Globaltag=="DN")
    FishersMethod.q.up.df<-subset(FishersMethod.q.df,Globaltag=="UP")
    FishersMethod.q.dn.df$rank<-rank(FishersMethod.q.dn.df$qvalueInSample)
    FishersMethod.q.up.df$rank<-rank(FishersMethod.q.up.df$qvalueInSample)
    FishersMethod.q.dn.df$InSampleTag<-ifelse(FishersMethod.q.dn.df$qvalueInSample<0.01,"Intra-donor","Inter-donor")
    FishersMethod.q.up.df$InSampleTag<-ifelse(FishersMethod.q.up.df$qvalueInSample<0.01,"Intra-donor","Inter-donor")
    # pdf(paste(outputname, "pdf",sep='.'),useDingbats=F)
    # ggplot(FishersMethod.q.dn.df)+aes(rank,-log10(qvalueInSample),label=GeneLabel,shape=InSampleTag)+geom_point(color="red")+geom_text_repel(color="red")+geom_hline(yintercept=2,linetype=2)+theme_classic()+scale_shape_manual(values=c(20,2))
    # ggplot(FishersMethod.q.up.df)+aes(rank,-log10(qvalueInSample),label=GeneLabel,shape=InSampleTag)+geom_point(color="darkgreen")+geom_text_repel(color="darkgreen")+geom_hline(yintercept=2,linetype=2)+theme_classic()+scale_shape_manual(values=c(20,2))
    # ggplot(FishersMethod.q.df)+aes(-log10(qvalue),-log10(qvalueInSample),color=Globaltag,label=GeneLabel)+geom_point(size=0.5)+geom_text_repel(aes(color=Globaltag))+scale_color_manual(values=c("grey","red","green"))+geom_hline(yintercept=2,linetype=2)+geom_vline(xintercept=-log10(0.005),linetype=2)+theme_classic()+xlab("Global qvalue")
    # ggplot(FishersMethod.q.df)+aes(-log10(qvalue),-log10(qvalueInSample),color=Globaltag,label=GeneLabel)+geom_text_repel(aes(color=Globaltag))+scale_color_manual(values=c("grey","red","green"))+geom_hline(yintercept=2,linetype=2)+geom_vline(xintercept=-log10(0.005),linetype=2)+theme_classic()+xlab("Global qvalue")
    # dev.off()
    DN.hetero<-names(FishersMethod.q[bin20.g.dn])[FishersMethod.q[bin20.g.dn]<=0.01]
    DN.homo<-names(FishersMethod.q[bin20.g.dn])[FishersMethod.q[bin20.g.dn]>0.01]
    UP.hetero<-names(FishersMethod.q[bin20.g.up])[FishersMethod.q[bin20.g.up]<=0.01]
    UP.homo<-names(FishersMethod.q[bin20.g.up])[FishersMethod.q[bin20.g.up]>0.01]
    geneLis <- list(DN.hetero=DN.hetero, DN.homo=DN.homo, UP.hetero=UP.hetero, UP.homo=UP.homo)
    # geneLis <- c("FOS","RPL36AL","CDKN2A","ID1","DLK1","SPP1","PDE4B","SIX3","RGS16","ID3","LSAMP","INHBA","NPY","NEAT1","PPP1R1A","S100A10","FXYD2","A1CF","RBP4","FFAR4")
    # InSample.compare.nb(OBJ, pseudo=PCAInfo[,"pseudo.index.balanced",drop=F], gene=gene, GetDatatoplot=T)[[3]] %>% ggplot(.)+geom_violin(aes(Q,pseudo.index.balanced,fill=Ave.expr),alpha=0.9) + facet_grid(cols = vars(Sample),space="free")+scale_fill_gradient(low="white",high="red")+theme_classic()+ggtitle(gene)
    return(list(PCAInfo=PCAInfo,PCAInfo.20bin.ob.PCA=PCAInfo.20bin.ob.PCA, Global.bin20.g.up=bin20.g.up, Global.bin20.g.dn=bin20.g.dn, FishersMethod.q.df=FishersMethod.q.df, FishersMethod.q.dn.df=FishersMethod.q.dn.df, FishersMethod.q.up.df=FishersMethod.q.up.df,geneLis_qval001=geneLis))
}

