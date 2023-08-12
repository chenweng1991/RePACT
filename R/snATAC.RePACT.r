
#' MakeEvenBinBydepth_SpeedUP_ATAC
#'
#' This function is to make even bins for cells, to keep the cell number and total cell fragments even.
#' @param OBJ.tmp, Seurat OBJ
#' @param data.info
#' @param binnumber
#' @return The function return a list: # data.info.withbin:(data.info + evenfragbin), cellvsPeak.m.aggr: gene ~ traj1-traj20 , index: mean pseudoindex of every traj_bin, depth: total frags of cells within each traj_bin
#' @import Seurat pscl rlist plyr
#' @export
#' @examples
#' beta.RNA.PCA.20bin.ob <- MakeEvenBinBydepth_SpeedUP(OBJ=OBJ, data.info=BetaPCA[,51:ncol(BetaPCA)], binnumber=20)

MakeEvenBinBydepth_SpeedUP_ATAC <- function(OBJ.tmp, data.info=BetaPeak.data.info, binnumber=20){
    cellvsPeak.m <- OBJ.tmp@assays$peak@counts[,row.names(data.info[order(data.info$rank),])]
    cell_frags <- colSums(cellvsPeak.m)
    cell_frags.binLis <- splitter(cell_frags, binnumber)
    names(cell_frags.binLis) <- 1:binnumber
    cellvsPeak.m.aggr.Lis <- list()
    for(tmp in 1:length(cell_frags.binLis)){
          cell_frags.binLis[[tmp]] <- data.frame(cell=names(cell_frags.binLis[[tmp]]), frags=cell_frags.binLis[[tmp]], evenfragbin=tmp)
          cellvsPeak.m.aggr.Lis[[paste("traj",tmp,sep='')]] <- rowSums(OBJ.tmp@assays$peak@counts[,cell_frags.binLis[[tmp]][,"cell"]])
    }
    data.info.withbin <- merge(data.info, list.rbind(cell_frags.binLis), by.x=0, by.y='cell',all.x=TRUE)
    rownames(data.info.withbin) <- data.info.withbin$Row.names
    data.info.withbin <- data.info.withbin[order(data.info.withbin$evenfragbin),]
    data.info.withbin$evenfragbin <- factor(data.info.withbin$evenfragbin, levels=unique(data.info.withbin$evenfragbin))
    index<-ddply(as.data.frame(data.info.withbin),.(evenfragbin),summarise,mean(pseudo.index.balanced))[,2]
    cellvsPeak.m.aggr <- t(list.rbind(cellvsPeak.m.aggr.Lis))
    depths <- colSums(cellvsPeak.m.aggr)
    return(list(data.info.withbin=data.info.withbin,cellvsPeak.m.aggr=cellvsPeak.m.aggr,index=index,depths=depths))
}

#' CallT2Dpeak_qvalue_SpeedUP_ATAC
#'
#' This function is to call dynamic peaks along the phenotype trajectory
#' @param cellvsPeak.m.aggr, a vector of peaks (":" "-")
#' @param depths, RePACT reesult OBJ
#' @param index, RePACT reesult OBJ
#' @param qcut, qvalue cutoff, default qcut=0.1
#' @param slopecut1, slope cutoff for UP-regulated genes/peaks, default slopecut1=0.5
#' @param slopecut2, slope cutoff for DN-regulated genes/peaks,default slopecut2=-0.5
#' @param doscale, if the regression is performed on scaled index and expr data, default doscale=T
#' @return The function return a list of plots, pseudoregress.all contains summary statistics for all genes/peaks, 
#' @import ggplot2
#' @export
#' @examples

CallT2Dpeak_qvalue_SpeedUP_ATAC <- function(cellvsPeak.m.aggr=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$cellvsPeak.m.aggr, depths=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$depths, index=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$index,qcut=0.1,slopecut1=0.5,slopecut2=-0.5,doscale=T){
    require(parallel)
    require(qvalue)
    cellvsPeak.m.aggr.norm <- t(t(cellvsPeak.m.aggr)/depths) # normalize by total cell depth per traj_bin, result is gene~traj_bin
    cellvsPeak.m.aggr.norm.scale <- apply(cellvsPeak.m.aggr.norm,1,function(x){(x-mean(x))/sd(x)}) %>% t # standardize, result is gene~traj_bin
    # cellvsPeak.m.aggr.scale<- apply(cellvsPeak.m.aggr,1,function(x){(x-mean(x))/sd(x)}) %>% t
    index.scale <- scale(index) # standadize avg_pseudo-index
    Get_stats <- function(gene, cellvsPeak.m.aggr, cellvsPeak.m.aggr.norm.scale, index.scale){
        if(any(is.na(as.numeric(cellvsPeak.m.aggr.norm.scale[gene,])))){
            return(c(NA,NA,NA,NA,NA))
        }else{
            gmodel1<-glm(as.numeric(cellvsPeak.m.aggr.norm.scale[gene,] )~index.scale,family=gaussian) #build generalized linear models, traj_scale_expression~index_scale
            cur.slope1 <- summary(gmodel1)$coefficients[,1][2]   # to get slope of pseudoBMIindex
            cur.pvalues1 <- summary(gmodel1)$coefficients[,4][2]  # to get p value for pseudoBMIindex
            cur.corr <- cor(1e4*as.numeric(cellvsPeak.m.aggr.norm.scale[gene,]),index.scale) # correlation between traj_expr and index
            cur.MaxMedian.FC <- max(cellvsPeak.m.aggr.norm[gene,])/median(cellvsPeak.m.aggr.norm[gene,])
            cur.Max <- 1e6*max(cellvsPeak.m.aggr.norm[gene,])
            names(cur.slope1)<-"slope"
            names(cur.pvalues1)<-"pvalue"
            names(cur.corr)<-"cor"
            names(cur.MaxMedian.FC)<-"MaxMedianFC"
            names(cur.Max)<-"MaxRPKM"
            return(c(cur.slope1,cur.pvalues1,cur.corr,cur.MaxMedian.FC,cur.Max))
        }
    }
    if(doscale==T){
        startTime <- Sys.time()
        pseudoregress.all.lis <- mclapply(rownames(cellvsPeak.m.aggr), Get_stats, cellvsPeak.m.aggr, cellvsPeak.m.aggr.norm.scale, index.scale, mc.cores = 40)
        endTime <- Sys.time()
        print(endTime - startTime)
    }else{
        startTime <- Sys.time()
        pseudoregress.all.lis <- mclapply(rownames(cellvsPeak.m.aggr), Get_stats, cellvsPeak.m.aggr, cellvsPeak.m.aggr.norm, index, mc.cores = 40)
        endTime <- Sys.time()
        print(endTime - startTime)
    }
    pseudoregress.all <- as.data.frame(list.rbind(pseudoregress.all.lis))
    row.names(pseudoregress.all) <- row.names(cellvsPeak.m.aggr)
    pseudoregress.all$qvalue <- qvalue(pseudoregress.all$pvalue)$qvalues # qvalue
    pseudoregress.all <- pseudoregress.all[complete.cases(pseudoregress.all),]
    UP <- subset(pseudoregress.all, qvalue<qcut & slope>slopecut1) %>% .[order(.$slope,decreasing=T),]
    DN <- subset(pseudoregress.all, qvalue<qcut & slope<slopecut2) %>% .[order(.$slope,decreasing=F),]
    UPDN.toplot<-t(cellvsPeak.m.aggr)/depths
    return(list(pseudoregress.all=pseudoregress.all,UP=UP,DN=DN,UPDN.toplot=UPDN.toplot))
}
#' CallT2Dpeak_pvalueOneTail_SpeedUP
#' @param cellvsPeak.m.aggr
#' @param depths
#' @param index
#' @param doscale
#' @param GlobalSlopes
#' @return 
#' @export
#' @examples PCAInfo.20bin.donor.PCA <- CallT2Dpeak_pvalueOneTail_SpeedUP(PCAInfo.20bin.donor.ob$cellvsPeak.m.aggr, PCAInfo.20bin.donor.ob$depths, PCAInfo.20bin.donor.ob$index,doscale=T, GlobalSlopes=PCAInfo.20bin.ob.PCA$pseudoregress.all[,1,drop=F])

CallT2Dpeak_pvalueOneTail_SpeedUP <- function(cellvsPeak.m.aggr, depths, index, doscale=T, GlobalSlopes){
    cellvsPeak.m.aggr.norm <- t(t(cellvsPeak.m.aggr)/depths)
    cellvsPeak.m.aggr.norm.scale <- apply(cellvsPeak.m.aggr.norm,1,function(x){(x-mean(x))/sd(x)}) %>% t
    # cellvsPeak.m.aggr.scale<- apply(cellvsPeak.m.aggr,1,function(x){(x-mean(x))/sd(x)}) %>% t
    index.scale <- scale(index)
    Get_stats <- function(gene, cellvsPeak.m.aggr, cellvsPeak.m.aggr.norm.scale, index.scale){
        if(any(is.na(as.numeric(cellvsPeak.m.aggr.norm.scale[gene,])))){
            return(c(NA,NA,NA,NA,NA))
        }else{
            gmodel1<-glm(as.numeric(cellvsPeak.m.aggr.norm.scale[gene,] )~index.scale,family=gaussian) #build generalized linear models, traj_scale_expression~index_scale
            OneTail.p <- pt(summary(gmodel1)$coefficients[2,3], gmodel1$df.residual, lower.tail=GlobalSlopes[row.names(cellvsPeak.m.aggr.norm.scale)[gene],]<0) # individual slope vs global slope
            cur.slope1 <- summary(gmodel1)$coefficients[,1][2]   # to get slope of pseudoBMIindex
            cur.pvalues1 <- OneTail.p  # to get p value for pseudoBMIindex, which is one tailed against global
            cur.corr <- cor(1e4*as.numeric(cellvsPeak.m.aggr.norm.scale[gene,]),index.scale)
            cur.MaxMedian.FC <- max(cellvsPeak.m.aggr.norm[gene,])/median(cellvsPeak.m.aggr.norm[gene,])
            cur.Max <- 1e6*max(cellvsPeak.m.aggr.norm[gene,])
            names(cur.slope1) <- "slope"
            names(cur.pvalues1) <- "pvalue.onetail"
            names(cur.corr) <- "cor"
            names(cur.MaxMedian.FC) <- "MaxMedianFC"
            names(cur.Max) <- "MaxRPKM"
            return(c(cur.slope1,cur.pvalues1,cur.corr,cur.MaxMedian.FC,cur.Max))
        }
    }
    if(doscale==T){
        startTime <- Sys.time()
        pseudoregress.all.lis <- mclapply(rownames(cellvsPeak.m.aggr), Get_stats, cellvsPeak.m.aggr, cellvsPeak.m.aggr.norm.scale, index.scale, mc.cores = 40)
        endTime <- Sys.time()
        print(endTime - startTime)
    }else{
        startTime <- Sys.time()
        pseudoregress.all.lis <- mclapply(rownames(cellvsPeak.m.aggr), Get_stats, cellvsPeak.m.aggr, cellvsPeak.m.aggr.norm, index, mc.cores = 40)
        endTime <- Sys.time()
        print(endTime - startTime)
    }
    pseudoregress.all <- as.data.frame(list.rbind(pseudoregress.all.lis))
    row.names(pseudoregress.all) <- row.names(cellvsPeak.m.aggr)
    return(pseudoregress.all)
}

#' snATAC.RePACT
#'
#' This function is to run logistic regression based on the number of LSIs, and characteristics of samples.
#' @param OBJ, a snATAC-seq Seurat object, assay name is "peak"
#' @param Sample, colnames of donor or sample infomation in OBJ@meta.data
#' @param pheno, the column name of the "characteristics to compare" in OBJ@meta.data, e.g. diseaseStatus or BMI or HBA1C
#' @param pheno_levels, specify the levels of pheno column. e.g. c("HT", "T2D"), specify for binary pheno
#' @param is_continuous, if pheno is continous variable. default is F. e.g. diseaseStatus is F, BMI is T
#' @param if_donorWise, if perform donor wise RePACT, default is F
#' @param binnumber, number of bins for cell grouping
#' @param LSIrange. default is "", automatically take top 10 LSIs that are significant with phenotype, otherwise specify 10 LSIs, e.g. 1:10 or 2:11
#' @param RePACT_qvalCut, qvalue cutoff to determine the significant genes along the disease or other characteristics trajactory, certain cutoff for slope can be further applied
#' @param donorWise_qvalCut, qvalue cutoff to determine the significant genes varying within donors or across donors. 
#' @return, a list of objects
#' @import Matrix.utils Signac Seurat plyr dplyr ggrepel ggplot2 qvalue reshape2
#' @export
#' @examples
#' snATAC.RePACT=(OBJ, Sample, pheno, is_continuous=F, if_donorWise=F, RePACT_qvalCut=0.005, donorWise_qvalCut=0.01){

snATAC.RePACT <- function(OBJ, Sample, pheno, pheno_levels, is_continuous=F, if_donorWise=F, binnumber=20, LSIrange="", RePACT_qvalCut=0.01, donorWise_qvalCut=0.01){
    require(Matrix.utils)
    require(Signac)
    require(Seurat)
    require(pscl)
    require(reshape2)
    GetRePACTLinearmodel.cca<-function(ccaWithinfo=cca.L2.info, prefix="CC",pheno="Disease",CCrange=1:10){
        CCnames <- paste("ccaWithinfo$",prefix,"_",CCrange, sep = "")
        ccaWithinfo[,pheno] <- as.factor(ccaWithinfo[,pheno])
#        ccaWithinfo[,pheno] <- as.numeric(levels(ccaWithinfo[,pheno]))[ccaWithinfo[,pheno]]
        form <- formula(paste("ccaWithinfo[,pheno]", paste(CCnames, collapse = "+"), sep = "~"))
        model <- lm(form)
        return(model)
    }
    Tomerge_v2 <- function (A, B, leavex = T, leavey = F){
        mergeAB <- merge(A, B, by = "row.names", all.x = leavex, all.y = leavey)
        row.names(mergeAB) <- mergeAB[, 1]
        mergeAB <- mergeAB[, -1]
        return(mergeAB)
    }

    LSIInfo <- OBJ@reductions$lsi@cell.embeddings %>% apply(.,2,function(x){x/sqrt(sum(x^2))}) %>% merge(., OBJ@meta.data, by=0)
    rownames(LSIInfo) <- LSIInfo$Row.names
    LSIInfo <- LSIInfo[,-1]
    if(is_continuous==F){
        LSIInfo[,pheno] <- as.character(LSIInfo[,pheno])
        LSIInfo[,pheno] <- factor(LSIInfo[,pheno] ,levels=pheno_levels)
    }
    LSIInfo.subs<-list()
    for(i in 1:100){
        LSIInfo.sub<-c()
        for (d in unique(LSIInfo[,Sample])){
            tmp <- LSIInfo[which(LSIInfo[,Sample]==d),]
            if(nrow(tmp)<=500){
                LSIInfo.sub<-rbind(LSIInfo.sub,tmp)
            }else{
                LSIInfo.sub<-rbind(LSIInfo.sub,tmp[sample(1:nrow(tmp),500),])
            }
        }
    LSIInfo.subs<-c(LSIInfo.subs,list(LSIInfo.sub))
    }
    pseudo.indexes<-list()
    if(LSIrange==""){
        sigLSIs <- SelectSigPCs(LSIInfo, Sample, pheno, is_continuous=F)
        sigLSIs_number <- as.numeric(sigLSIs[,1])
    }else{
        sigLSIs_number <- LSIrange
    }
    for(i in 1:100){
        if(is_continuous==F){
            md <- GetRePACTmodel.cca(ccaWithinfo=LSIInfo.subs[[i]],prefix="LSI",pheno=pheno,CCrange=c(sigLSIs_number))
        }else{
            md <- GetRePACTLinearmodel.cca(ccaWithinfo=LSIInfo.subs[[i]],prefix="LSI",pheno=pheno,CCrange=c(sigLSIs_number))
        }
        trainingdata <- LSIInfo.subs[[i]]
        Restdata <- LSIInfo[setdiff(row.names(LSIInfo),row.names(LSIInfo.subs[[i]])),]
        Alldata <- LSIInfo
        LSIInfo.subs[[i]]$pseudo.index<-apply(trainingdata[,sigLSIs_number],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        Restdata$pseudo.index<-apply(Restdata[,sigLSIs_number],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        pseudo.indexes<-c(pseudo.indexes,list(apply(Alldata[,sigLSIs_number],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})))
    }
    md.all<-GetRePACTmodel.cca(ccaWithinfo=LSIInfo, prefix="LSI", pheno=pheno, CCrange=c(sigLSIs_number))
    LSIInfo$pseudo.index = md.all$linear.predictors
    LSIInfo$pseudo.index.balanced<-do.call(cbind,pseudo.indexes) %>% rowMeans()
    LSIInfo$rank<-rank(LSIInfo$pseudo.index.balanced)
    LSIInfo.20bin.ob <- MakeEvenBinBydepth_SpeedUP_ATAC(OBJ.tmp=OBJ, data.info=LSIInfo[,51:ncol(LSIInfo)], binnumber=binnumber)
    LSIInfo.20bin.ob.LSI <- CallT2Dpeak_qvalue_SpeedUP_ATAC(LSIInfo.20bin.ob$cellvsPeak.m.aggr, LSIInfo.20bin.ob$depths, LSIInfo.20bin.ob$index, qcut=RePACT_qvalCut, slopecut1=0.5, slopecut2=-0.5, doscale=T)
    bin20.g.up <- row.names(subset(LSIInfo.20bin.ob.LSI$UP,qvalue<RePACT_qvalCut & MaxRPKM>quantile(LSIInfo.20bin.ob.LSI$pseudoregress.all$MaxRPKM,0.25)))
    bin20.g.dn <- row.names(subset(LSIInfo.20bin.ob.LSI$DN,qvalue<RePACT_qvalCut & MaxRPKM>quantile(LSIInfo.20bin.ob.LSI$pseudoregress.all$MaxRPKM,0.25)))
    RePACT_call <- list(UP=bin20.g.up, DN=bin20.g.dn)
    # Evenbin.donorContribute<- LSIInfo.20bin.ob$data.info.withbin %>% .[,c(Sample,"evenfragbin")] %>% table %>% as.matrix %>% apply(.,1,function(x){x/sum(x)}) %>% melt
    # Evenbin.donorContribute <- merge(Evenbin.donorContribute, unique(LSIInfo.20bin.ob$data.info.withbin[,c(Sample,pheno)]), by=Sample, all.x=TRUE)
    if(if_donorWise==T){
        RePACT.diff.peaks.bydonor<-list()
        RePACT.diff.peaks.bydonor.names <- c()
        for (curDonor in unique(LSIInfo[,Sample])){
            print(curDonor)
            BetaLSI.donor <- LSIInfo[which(LSIInfo[,Sample]==curDonor),]
            if(nrow(BetaLSI.donor)>=100){
                print(paste("run RePACT for donor",curDonor))
                RePACT.diff.peaks.bydonor.names <- c(RePACT.diff.peaks.bydonor.names, curDonor)
                BetaLSI.donor$rank <- rank(BetaLSI.donor$pseudo.index.balanced)
                beta.ATAC.LSI.20bin.donor.ob <- MakeEvenBinBydepth_SpeedUP_ATAC(OBJ=OBJ, data.info=BetaLSI.donor[,51:ncol(LSIInfo)], binnumber=binnumber)
                betaT2D.diffPeak.20bin.donor <- CallT2Dpeak_pvalueOneTail_SpeedUP(beta.ATAC.LSI.20bin.donor.ob$cellvsPeak.m.aggr, beta.ATAC.LSI.20bin.donor.ob$depths, beta.ATAC.LSI.20bin.donor.ob$index, doscale=T, GlobalSlopes=LSIInfo.20bin.ob.LSI$pseudoregress.all[,1,drop=F])
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
        return(list(LSIInfo=LSIInfo, LSIInfo.20bin.ob=LSIInfo.20bin.ob, LSIInfo.20bin.ob.LSI=LSIInfo.20bin.ob.LSI, RePACT_call=RePACT_call, RePACT_donorWise_intermediate=RePACT_donorWise_intermediate, RePACT_donorWise_call=RePACT_donorWise_call))
    }else{
        return(list(LSIInfo=LSIInfo, LSIInfo.20bin.ob=LSIInfo.20bin.ob, LSIInfo.20bin.ob.LSI=LSIInfo.20bin.ob.LSI, RePACT_call=RePACT_call))
    }
}


#' plot_scATAC_RePACT_heatmap
#'
#' This function is to take a vetor of peaks and plot heatmaps from RePACT result
#' @param peakSet, a vector of peaks (":" "-")
#' @param RePACT_OBJ, RePACT reesult OBJ
#' @return The function return a plot
#' @import ggplot2
#' @export
#' @examples

plot_scATAC_RePACT_heatmap <- function(peakSet, RePACT_OBJ){
	   normalize_01<-function(vector){
            normed<-(vector-min(vector))/(max(vector)-min(vector))
            return(normed)
    }

    T2D.bindata <- RePACT_OBJ$LSIInfo.20bin.ob.LSI$UPDN.toplot[, intersect(peakSet, colnames(RePACT_OBJ$LSIInfo.20bin.ob.LSI$UPDN.toplot))]
    T2D.bindata <- data.frame(apply(T2D.bindata[, -ncol(T2D.bindata)-1, drop = F],2, normalize_01), bin = 1:nrow(T2D.bindata))
    T2D.bindata.m <- reshape2::melt(T2D.bindata, id.vars = "bin")
    T2D.bindata.m$bin <- factor(T2D.bindata.m$bin,levels=1:nrow(T2D.bindata))
    T2D.bindata.m$variable <- factor(T2D.bindata.m$variable, levels=gsub("-",".",gsub(":",".",peakSet)))
    T2D.bindata.m <- T2D.bindata.m[which(!is.na(T2D.bindata.m$variable)),]
    p1 <- ggplot(T2D.bindata.m) + aes(bin, variable, fill = value) + geom_tile() + scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 0.6) + 
            labs(y = "Variable peaks")+theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ theme(legend.position = "none")
    return(p1)
}

