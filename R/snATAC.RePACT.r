#' snATAC.RePACT
#'
#' This function is to run logistic regression based on the number of LSIs, and characteristics of samples.
#' @param OBJ, a scRNA-seq Seurat object
#' @param Sample, colnames of donor or sample infomation in OBJ@meta.data
#' @param pheno, the column name of the "characteristics to compare" in OBJ@meta.data, e.g. diseaseStatus or BMI or HBA1C
#' @param is_continuous, if pheno is continous variable. default is F. e.g. diseaseStatus is F, BMI is T
#' @param if_donorWise, if perform donor wise RePACT, default is F
#' @param RePACT_qvalCut, qvalue cutoff to determine the significant genes along the disease or other characteristics trajactory, certain cutoff for slope can be further applied
#' @param donorWise_qvalCut, qvalue cutoff to determine the significant genes varying within donors or across donors. 
#' @return, a list of objects
#' @import Seurat plyr dplyr Seurat ggrepel ggplot2 qvalue reshape2
#' @export
#' @examples
#' snATAC.RePACT=(OBJ, Sample, pheno, is_continuous=F, if_donorWise=F, RePACT_qvalCut=0.005, donorWise_qvalCut=0.01){


snATAC.RePACT <- function(OBJ, Sample, pheno, is_continuous=F, if_donorWise=F, RePACT_qvalCut=0.01, donorWise_qvalCut=0.01){
    GetRePACTLinearmodel.cca<-function(ccaWithinfo=cca.L2.info, prefix="CC",pheno="Disease",CCrange=1:10){
        CCnames <- paste("ccaWithinfo$",prefix,"_",CCrange, sep = "")
        ccaWithinfo[,pheno] <- as.factor(ccaWithinfo[,pheno])
        ccaWithinfo[,pheno] <- as.numeric(levels(ccaWithinfo[,pheno]))[ccaWithinfo[,pheno]]
        form <- formula(paste("ccaWithinfo[,pheno]", paste(CCnames, collapse = "+"), sep = "~"))
        model <- lm(form)
        return(model)
    }

    LSIInfo <- OBJ@reductions$lsi@cell.embeddings %>% apply(.,2,function(x){x/sqrt(sum(x^2))}) %>% merge(., OBJ@meta.data, by=0)
    rownames(LSIInfo) <- LSIInfo$Row.names
    LSIInfo <- LSIInfo[,-1]
    LSIInfo.subs<-list()
    for(i in 1:100){
    LSIInfo.sub<-c()
    for (d in levels(LSIInfo[,Sample])){
        tmp<-subset(LSIInfo, Sample==d)
        if(nrow(tmp)<=500){
        LSIInfo.sub<-rbind(LSIInfo.sub,tmp)
        }else{
        LSIInfo.sub<-rbind(LSIInfo.sub,tmp[sample(1:nrow(tmp),500),])
        }
    }
    LSIInfo.subs<-c(LSIInfo.subs,list(LSIInfo.sub))
    }
    pseudo.indexes<-list()
    for(i in 1:100){
        if(is_continuous==F){
            md <- GetRePACTmodel.cca(ccaWithinfo=LSIInfo.subs[[i]],prefix="LSI",pheno=pheno,CCrange=c(1:10))
        }else{
            md <- GetRePACTLinearmodel.cca(ccaWithinfo=LSIInfo.subs[[i]],prefix="LSI",pheno=pheno,CCrange=c(1:10))
        }
        trainingdata <- LSIInfo.subs[[i]]
        Restdata <- LSIInfo[setdiff(row.names(LSIInfo),row.names(LSIInfo.subs[[i]])),]
        Alldata <- LSIInfo
        LSIInfo.subs[[i]]$pseudo.index<-apply(trainingdata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        Restdata$pseudo.index<-apply(Restdata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        pseudo.indexes<-c(pseudo.indexes,list(apply(Alldata[,1:10],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})))
    }
    md.all<-GetRePACTmodel.cca(ccaWithinfo=LSIInfo,prefix="LSI",pheno=pheno,CCrange=c(1:10))
    LSIInfo$pseudo.index = md.all$linear.predictors
    LSIInfo$pseudo.index.balanced<-do.call(cbind,pseudo.indexes) %>% rowMeans()
    LSIInfo$rank<-rank(LSIInfo$pseudo.index.balanced)
    ## Plot violin pseudoindex
    # p1<-ggplot(LSIInfo)+aes(Sample,pseudo.index.balanced,fill=get(pheno))+geom_violin()+geom_boxplot(width=0.2,outlier.shape = NA,notch=F,coef = 0,fill="grey25",color="grey75")+coord_flip()+theme_classic()+scale_fill_manual(values=c("steelblue", "red"))+theme_bw()+theme(legend.position="none")
    # p2<-ggplot(LSIInfo)+aes(pseudo.index.balanced,fill=get(pheno))+geom_density(alpha=0.75)+scale_fill_manual(values=c("steelblue","red"))+theme_classic()+theme(legend.position="none")
    LSIInfo.20bin.ob <- MakeEvenBinBydepth(cellvsPeak.m=t(as.matrix(OBJ@assays$ATAC@counts)),data.info=LSIInfo,binnumber=20)
    LSIInfo.20bin.ob.LSI <- CallT2Dpeak_qvalue(LSIInfo.20bin.ob$cellvsPeak.m.aggr,LSIInfo.20bin.ob$depths,LSIInfo.20bin.ob$index,qcut=RePACT_qvalCut,slopecut1=0.5,slopecut2=-0.5,doscale=T)
    Evenbin.donorContribute<- LSIInfo.20bin.ob$data.info.withbin %>% .[,c("Sample","evenfragbin")] %>% table %>% as.matrix %>% apply(.,1,function(x){x/sum(x)}) %>% melt
    Evenbin.donorContribute <- merge(Evenbin.donorContribute, unique(LSIInfo.20bin.ob$data.info.withbin[,c("Sample",pheno)]), by="Sample", all.x=TRUE)
    # ggplot(Evenbin.donorContribute)+aes(Sample,value,fill=get(pheno))+geom_bar(stat="identity",color="black")+facet_grid(~evenfragbin)+theme(axis.text=element_blank(),axis.ticks=element_blank())+scale_fill_manual(values=c("steelblue","red"))+theme_bw()
    #   p1<-apply(LSIInfo.20bin.ob.LSI$UPDN.toplot,2,function(x){scale(x)}) %>% melt() %>% ggplot()+aes(Var1,Var2,fill=value)+geom_tile()+scale_fill_gradient2(low="steelblue",mid="white",high="red")+theme_classic()+theme(axis.text=element_blank())+ggtitle("LSI-20bins")+labs(x="trajectory")
    if(if_donorWise==T){
        RePACT.diff.peaks.bydonor<-list()
        RePACT.diff.peaks.bydonor.names <- c()
        for (curDonor in levels(LSIInfo[,"Sample"])){
            print(curDonor)
            BetaLSI.donor <- subset(LSIInfo, Sample==curDonor)
            if(nrow(BetaLSI.donor)>=50){
                RePACT.diff.peaks.bydonor.names <- c(RePACT.diff.peaks.bydonor.names, curDonor)
                BetaLSI.donor$rank <- rank(BetaLSI.donor$pseudo.index.balanced)
                beta.ATAC.LSI.20bin.donor.ob <- MakeEvenBinBydepth(cellvsPeak.m=t(as.matrix(OBJ@assays$ATAC@counts[,row.names(BetaLSI.donor)])),data.info=BetaLSI.donor[,51:ncol(BetaLSI.donor)],binnumber=20)
                betaT2D.diffPeak.20bin.donor <- CallT2Dpeak_pvalueOneTail(beta.ATAC.LSI.20bin.donor.ob$cellvsPeak.m.aggr,beta.ATAC.LSI.20bin.donor.ob$depths,beta.ATAC.LSI.20bin.donor.ob$index,doscale=T,GlobalSlopes=LSIInfo.20bin.ob.LSI$pseudoregress.all[,1,drop=F])
                RePACT.diff.peaks.bydonor <- c(RePACT.diff.peaks.bydonor,list(betaT2D.diffPeak.20bin.donor))
            }
        }
        names(RePACT.diff.peaks.bydonor) <- RePACT.diff.peaks.bydonor.names  
        ps <- lapply(RePACT.diff.peaks.bydonor,function(x){x$pvalue.onetail}) %>% do.call(cbind,.)
        row.names(ps) <- row.names(RePACT.diff.peaks.bydonor[[1]])
        ps[is.na(ps)]<-1
        slopes <- lapply(RePACT.diff.peaks.bydonor,function(x){x$slope}) %>% do.call(cbind,.)
        row.names(slopes) <- row.names(RePACT.diff.peaks.bydonor[[1]]$pseudoregress.all)
        FishersMethod.p <- apply(ps,1,function(x){-2*sum(log(x))}) %>% pchisq(.,2*11,lower.tail=FALSE)  ## k=11, df=2k
        FishersMethod.q<-qvalue(FishersMethod.p)$qvalues %>% .[order(.)]
        FishersMethod.q<-FishersMethod.q[c(row.names(LSIInfo.20bin.ob.LSI$UP),row.names(LSIInfo.20bin.ob.LSI$DN))] 
        FishersMethod.q.df<-data.frame(peak=names(FishersMethod.q),qvalueInSample=FishersMethod.q)
        FishersMethod.q.df<-merge(FishersMethod.q.df,rbind(LSIInfo.20bin.ob.LSI$UP, LSIInfo.20bin.ob.LSI$DN), by=0)
        rownames(FishersMethod.q.df) <- FishersMethod.q.df$Row.names
        FishersMethod.q.df <- FishersMethod.q.df[,-1]
        FishersMethod.q.dn.df<-subset(FishersMethod.q.df,slope<0) %>% cbind(.,rank=rank(.$qvalueInSample)) %>% .[order(.$rank),]
        FishersMethod.q.dn.df$tag1<-ifelse(FishersMethod.q.dn.df$qvalueInSample<0.05,"Intra-donor","NS")
        FishersMethod.q.dn.df$tag2<-ifelse(FishersMethod.q.dn.df$rank<=20,row.names(FishersMethod.q.dn.df),"")
        FishersMethod.q.up.df<-subset(FishersMethod.q.df,slope>0) %>% cbind(.,rank=rank(.$qvalueInSample)) %>% .[order(.$rank),]
        FishersMethod.q.up.df$tag1<-ifelse(FishersMethod.q.up.df$qvalueInSample<0.05,"Intra-donor","NS")
        FishersMethod.q.up.df$tag2<-ifelse(FishersMethod.q.up.df$rank<=20,row.names(FishersMethod.q.up.df),"")
        DN.hetero.peaks<-subset(FishersMethod.q.dn.df,qvalueInSample<donorWise_qvalCut) %>% row.names
        DN.homo.peaks<-subset(FishersMethod.q.dn.df,qvalueInSample>=donorWise_qvalCut) %>% row.names
        UP.hetero.peaks<-subset(FishersMethod.q.up.df,qvalueInSample<donorWise_qvalCut) %>% row.names
        UP.homo.peaks<-subset(FishersMethod.q.up.df,qvalueInSample>=donorWise_qvalCut) %>% row.names
        RePACT_donorWise_intermediate <- list(FishersMethod.q.df=FishersMethod.q.df, FishersMethod.q.up.df=FishersMethod.q.up.df, FishersMethod.q.dn.df=FishersMethod.q.dn.df)
        RePACT_donorWise_call <- list(DN.interDonorPeak=DN.hetero.peaks, DN.intraDonorPeak=DN.homo.peaks, UP.interDonorPeak=UP.hetero.peaks, UP.intraDonorPeak=UP.homo.peaks)
        return(list(LSIInfo=LSIInfo, LSIInfo.20bin.ob=LSIInfo.20bin.ob,LSIInfo.20bin.ob.LSI=LSIInfo.20bin.ob.LSI, RePACT_donorWise_intermediate=RePACT_donorWise_intermediate, RePACT_donorWise_call=RePACT_donorWise_call))
    }else{
        return(list(LSIInfo=LSIInfo, LSIInfo.20bin.ob=LSIInfo.20bin.ob, LSIInfo.20bin.ob.LSI=LSIInfo.20bin.ob.LSI))

    }
}
