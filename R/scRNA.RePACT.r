#' scRNA.RePACT
#'
#' This function is to run regression based on the number of PCs, and characteristics of samples.
#' @param OBJ, a scRNA-seq Seurat object
#' @param Sample, colnames of donor or sample infomation in OBJ@meta.data
#' @param pheno, the column name of the "characteristics to compare" in OBJ@meta.data, e.g. diseaseStatus or BMI or HBA1C
#' @param is_continuous, if pheno is continous variable. default is F. e.g. diseaseStatus is F, BMI is T
#' @param if_donorWise, if perform donor wise RePACT, default is F
#' @param RePACT_qvalCut, qvalue cutoff to determine the significant genes along the disease or other characteristics trajactory, certain cutoff for slope can be further applied
#' @param donorWise_qvalCut, qvalue cutoff to determine the significant genes varying within donors or across donors. 
#' @return The function return a list of objects 
#' @import Seurat plyr dplyr Seurat ggrepel ggplot2 plot3D pscl
#' @export
#' @examples
#' T2D.scRNA.RePACT <- scRNA.RePACT(OBJ=scRNA.OBJ,Sample="Donor", pheno="diseaseStat", is_continuous=F, if_donorWise=F)
scRNA.RePACT <- function(OBJ, Sample, pheno, is_continuous=F, if_donorWise=F, RePACT_qvalCut=0.005, donorWise_qvalCut=0.01){
    GetRePACTLinearmodel.cca<-function(ccaWithinfo=cca.L2.info, prefix="CC",pheno="Disease",CCrange=1:10){
        CCnames <- paste("ccaWithinfo$",prefix,"_",CCrange, sep = "")
        ccaWithinfo[,pheno] <- as.factor(ccaWithinfo[,pheno])
        ccaWithinfo[,pheno] <- as.numeric(levels(ccaWithinfo[,pheno]))[ccaWithinfo[,pheno]]
        form <- formula(paste("ccaWithinfo[,pheno]", paste(CCnames, collapse = "+"), sep = "~"))
        model <- lm(form)
        return(model)
    }
    BetaPCA<-OBJ@reductions$pca@cell.embeddings %>% Tomerge_v2(.,OBJ@meta.data)
    BetaPCA<-BetaPCA[,1:50] %>% apply(.,2,function(x){x/sqrt(sum(x^2))}) %>% cbind(.,BetaPCA[,51:ncol(BetaPCA)])
    BetaPCA[,Sample] <- factor(BetaPCA[,Sample] ,levels=unique(BetaPCA[,Sample] ))
    beta.rna.pca.withinfo.subs<-list()
    for(i in 1:100){
        beta.rna.pca.withinfo.sub<-c()
        for (d in levels(BetaPCA[,Sample])){
            tmp<- BetaPCA[which(BetaPCA[,Sample]==d),]
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
        if(is_continuous==F){
            md<-GetRePACTmodel.cca(ccaWithinfo=beta.rna.pca.withinfo.subs[[i]],prefix="PC",pheno=pheno,CCrange=c(1:10))
        }else{
            md<-GetRePACTLinearmodel.cca(ccaWithinfo=beta.rna.pca.withinfo.subs[[i]],prefix="PC",pheno=pheno,CCrange=c(1:10))
        }
        trainingdata<-beta.rna.pca.withinfo.subs[[i]]
        Restdata<-BetaPCA[setdiff(row.names(BetaPCA),row.names(beta.rna.pca.withinfo.subs[[i]])),]
        Alldata<-BetaPCA
        beta.rna.pca.withinfo.subs[[i]]$pseudo.index<-apply(trainingdata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        Restdata$pseudo.index<-apply(Restdata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        pseudo.indexes<-c(pseudo.indexes,list(apply(Alldata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})))
    }
    BetaPCA$pseudo.index.balanced<-do.call(cbind,pseudo.indexes) %>% rowMeans()
    adjustrange = seq(0, 180, length.out = 13)
    ks.test(BetaPCA[which(BetaPCA[,pheno]==unique(BetaPCA[,pheno])[1]),]$pseudo.index.balanced,BetaPCA[which(BetaPCA[,pheno]==unique(BetaPCA[,pheno])[2]),]$pseudo.index.balanced)
    BetaPCA$rank<-rank(BetaPCA$pseudo.index.balanced)
    beta.RNA.PCA.20bin.ob<-MakeEvenBinBydepth(cellvsPeak.m=t(as.matrix(OBJ@assays$RNA@counts)),data.info=BetaPCA[,51:ncol(BetaPCA)],binnumber=20)
    betaT2D.diffGene.20bin.PCA<-CallT2Dpeak_qvalue(beta.RNA.PCA.20bin.ob$cellvsPeak.m.aggr, beta.RNA.PCA.20bin.ob$depths, beta.RNA.PCA.20bin.ob$index, qcut=0.2,slopecut1=0.3,slopecut2=-0.3,doscale=T)
    bin20.g.up <- row.names(subset(betaT2D.diffGene.20bin.PCA$UP,qvalue<RePACT_qvalCut & MaxRPKM>quantile(betaT2D.diffGene.20bin.PCA$pseudoregress.all$MaxRPKM,0.25)))
    bin20.g.dn <- row.names(subset(betaT2D.diffGene.20bin.PCA$DN,qvalue<RePACT_qvalCut & MaxRPKM>quantile(betaT2D.diffGene.20bin.PCA$pseudoregress.all$MaxRPKM,0.25)))
    RePACT_call <- list(UP=bin20.g.up, DN=bin20.g.dn)
    if(if_donorWise==T){
        PCAInfo <- BetaPCA
        RePACT.diff.genes.bydonor<-list()
        for (curDonor in levels(PCAInfo[,Sample])){
            print(curDonor)
            PCAInfo.donor <- subset(PCAInfo, Sample==curDonor)
            PCAInfo.donor$rank <- rank(PCAInfo.donor$pseudo.index.balanced)
            PCAInfo.20bin.donor.ob <- MakeEvenBinBydepth(cellvsPeak.m=t(as.matrix(OBJ@assays$RNA@counts[,row.names(PCAInfo.donor)])),data.info=PCAInfo.donor[,51:ncol(PCAInfo.donor)],binnumber=20)
            PCAInfo.20bin.donor.PCA <- CallT2Dpeak_pvalueOneTail(PCAInfo.20bin.donor.ob$cellvsPeak.m.aggr, PCAInfo.20bin.donor.ob$depths, PCAInfo.20bin.donor.ob$index,doscale=T, GlobalSlopes=PCAInfo.20bin.ob.PCA$pseudoregress.all[,1,drop=F])
            RePACT.diff.genes.bydonor <- c(RePACT.diff.genes.bydonor,list(PCAInfo.20bin.donor.PCA))
        }
        names(RePACT.diff.genes.bydonor)<-levels(PCAInfo[,Sample])
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
        FishersMethod.q.dn.df$InSampleTag<-ifelse(FishersMethod.q.dn.df$qvalueInSample<donorWise_qvalCut,"Intra-donor","Inter-donor")
        FishersMethod.q.up.df$InSampleTag<-ifelse(FishersMethod.q.up.df$qvalueInSample<donorWise_qvalCut,"Intra-donor","Inter-donor")
        DN.hetero<-names(FishersMethod.q[bin20.g.dn])[FishersMethod.q[bin20.g.dn]<=donorWise_qvalCut]
        DN.homo<-names(FishersMethod.q[bin20.g.dn])[FishersMethod.q[bin20.g.dn]>donorWise_qvalCut]
        UP.hetero<-names(FishersMethod.q[bin20.g.up])[FishersMethod.q[bin20.g.up]<=donorWise_qvalCut]
        UP.homo<-names(FishersMethod.q[bin20.g.up])[FishersMethod.q[bin20.g.up]>donorWise_qvalCut]
        RePACT_donorWise_intermediate <- list(FishersMethod.q.df=FishersMethod.q.df, FishersMethod.q.up.df=FishersMethod.q.up.df, FishersMethod.q.dn.df=FishersMethod.q.dn.df)
        RePACT_donorWise_call <- list(DN.interdonor=DN.hetero, DN.intradonor=DN.homo, UP.interdonor=UP.hetero, UP.intradonor=UP.homo)
        return(list(BetaPCA=BetaPCA,beta.RNA.PCA.20bin.ob=beta.RNA.PCA.20bin.ob, betaT2D.diffGene.20bin.PCA=betaT2D.diffGene.20bin.PCA, RePACT_call=RePACT_call, RePACT_donorWise_intermediate=RePACT_donorWise_intermediate, RePACT_donorWise_call=RePACT_donorWise_call))
    }else{
        return(list(BetaPCA=BetaPCA,beta.RNA.PCA.20bin.ob=beta.RNA.PCA.20bin.ob, betaT2D.diffGene.20bin.PCA=betaT2D.diffGene.20bin.PCA, RePACT_call=RePACT_call))
    }
}
