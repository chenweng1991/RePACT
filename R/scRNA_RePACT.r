#' MakeEvenBinBydepth
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

MakeEvenBinBydepth_SpeedUP <- function(OBJ.tmp, data.info=BetaPeak.data.info,binnumber=20){
    # require(rlist)
    # require(plyr)
    splitter <- function(values, N){
        inds = c(0, sapply(1:N, function(i) which.min(abs(cumsum(as.numeric(values)) - sum(as.numeric(values))/N*i))))
        dif = diff(inds)
        re = rep(1:length(dif), times = dif)
        return(split(values, re))
    }
			   print("MakeEvenBinBydepth_SpeedUP")
    cellvsPeak.m <- OBJ.tmp@assays$RNA@counts[,row.names(data.info[order(data.info$rank),])]
    print(nrow(cellvsPeak.m))
    print(ncol(cellvsPeak.m))
    cell_frags <- colSums(cellvsPeak.m)
    cell_frags.binLis <- splitter(cell_frags, binnumber)
    names(cell_frags.binLis) <- 1:binnumber
    cell_frags.binLis.copy <- cell_frags.binLis		   
    cellvsPeak.m.aggr.Lis <- list()
    for(tmp in 1:length(cell_frags.binLis)){
          cell_frags.binLis.copy[[tmp]] <- data.frame(cell=names(cell_frags.binLis[[tmp]]), frags=cell_frags.binLis[[tmp]], evenfragbin=tmp)
          cellvsPeak.m.aggr.Lis[[paste("traj",tmp,sep='')]] <- rowSums(OBJ@assays$RNA@counts[,cell_frags.binLis.copy[[tmp]][,"cell"]])
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


#' CallT2Dpeak_qvalue
#'
#' This function is to run regression (linear/logistic) based on the number of PCs, and characteristics of samples.
#' @param cellvsPeak.m.aggr, gene~traj_bin counts matrix
#' @param depths total fragments of cells within each traj_bin
#' @param index  avg pseudo-index of each bin
#' @param qcut   qvalue cutoff to measure the significance of glm(expr~index)
#' @param slopecut1 UP slope cutoff, slope from glm(expr~index), UP genes
#' @param slopecut2 DN slope cutoff, slope from glm(expr~index), DN genes
#' @param doscale if build glm model based on scaled_expression and scaled pseudo-inddex
#' @return The function return a list: "PCvariance", "PCanfpheno", "object.raw.withinfo", "model", "reg.plot.2d", "model.para"
#' @import Seurat ggplot2 Matrix RColorBrewer gridExtra pscl dplyr
#' @export
#' @examples betaT2D.diffGene.20bin.PCA <- CallT2Dpeak_qvalue(beta.RNA.PCA.20bin.ob$cellvsPeak.m.aggr, depths=beta.RNA.PCA.20bin.ob$depths, index=beta.RNA.PCA.20bin.ob$index, qcut=0.2,slopecut1=0.3,slopecut2=-0.3,doscale=T)


CallT2Dpeak_qvalue <- function(cellvsPeak.m.aggr=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$cellvsPeak.m.aggr, depths=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$depths, index=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$index,qcut=0.1,slopecut1=0.5,slopecut2=-0.5,doscale=T){
    cellvsPeak.m.aggr.norm <- t(t(cellvsPeak.m.aggr)/depths) # normalize by total cell depth per traj_bin, result is gene~traj_bin
    cellvsPeak.m.aggr.norm.scale <- apply(cellvsPeak.m.aggr.norm,1,function(x){(x-mean(x))/sd(x)}) %>% t # standardize, result is gene~traj_bin
    # cellvsPeak.m.aggr.scale<- apply(cellvsPeak.m.aggr,1,function(x){(x-mean(x))/sd(x)}) %>% t
    index.scale <- scale(index) # standadize avg_pseudo-index
    pseudoregress.all<-c()
    for (i in 1:nrow(cellvsPeak.m.aggr)){  #Loop gene by gene
        if(any(is.na(as.numeric(cellvsPeak.m.aggr.norm.scale[i,])))){
            pseudoregress.all<-rbind(pseudoregress.all,c(NA,NA,NA,NA,NA))
            next
        }
        if(doscale){
            gmodel1<-glm(as.numeric(cellvsPeak.m.aggr.norm.scale[i,] )~index.scale,family=gaussian) #build generalized linear models, traj_scale_expression~index_scale
        }else{
            gmodel1<-glm(as.numeric(cellvsPeak.m.aggr.norm[i,] )~index,family=gaussian) #build generalized linear models, traj_norm_expression~index
        }
        cur.slope1 <- summary(gmodel1)$coefficients[,1][2]   # to get slope of pseudoBMIindex
        cur.pvalues1 <- summary(gmodel1)$coefficients[,4][2]  # to get p value for pseudoBMIindex
        cur.corr <- cor(1e4*as.numeric(cellvsPeak.m.aggr.norm.scale[i,]),index.scale) # correlation between traj_expr and index
        cur.MaxMedian.FC <- max(cellvsPeak.m.aggr.norm[i,])/median(cellvsPeak.m.aggr.norm[i,])
        cur.Max <- 1e6*max(cellvsPeak.m.aggr.norm[i,])
        names(cur.slope1)<-"slope"
        names(cur.pvalues1)<-"pvalue"
        names(cur.corr)<-"cor"
        names(cur.MaxMedian.FC)<-"MaxMedianFC"
        names(cur.Max)<-"MaxRPKM"
        pseudoregress.all <- rbind(pseudoregress.all,c(cur.slope1,cur.pvalues1,cur.corr,cur.MaxMedian.FC,cur.Max))
    }
    row.names(pseudoregress.all) <- row.names(cellvsPeak.m.aggr)
    pseudoregress.all <- as.data.frame(pseudoregress.all)
    pseudoregress.all$qvalue <- qvalue(pseudoregress.all$pvalue)$qvalues # qvalue
    pseudoregress.all <- pseudoregress.all[complete.cases(pseudoregress.all),]
    UP <- subset(pseudoregress.all, qvalue<qcut & slope>slopecut1) %>% .[order(.$slope,decreasing=T),]
    DN <- subset(pseudoregress.all, qvalue<qcut & slope<slopecut2) %>% .[order(.$slope,decreasing=F),]
    UPDN.toplot<-t(cellvsPeak.m.aggr)/depths
    return(list(pseudoregress.all=pseudoregress.all,UP=UP,DN=DN,UPDN.toplot=UPDN.toplot))
}

#' CallT2Dpeak_qvalue_SpeedUP
#' This function is to run regression (linear/logistic) based on the number of PCs, and characteristics of samples.
#' @param cellvsPeak.m.aggr, gene~traj_bin counts matrix
#' @param depths total fragments of cells within each traj_bin
#' @param index  avg pseudo-index of each bin
#' @param qcut   qvalue cutoff to measure the significance of glm(expr~index)
#' @param slopecut1 UP slope cutoff, slope from glm(expr~index), UP genes
#' @param slopecut2 DN slope cutoff, slope from glm(expr~index), DN genes
#' @param doscale if build glm model based on scaled_expression and scaled pseudo-inddex
#' @return The function return a list: "PCvariance", "PCanfpheno", "object.raw.withinfo", "model", "reg.plot.2d", "model.para"
#' @import Seurat ggplot2 Matrix RColorBrewer gridExtra pscl dplyr
#' @export
#' @examples betaT2D.diffGene.20bin.PCA <- CallT2Dpeak_qvalue_SpeedUP(beta.RNA.PCA.20bin.ob$cellvsPeak.m.aggr, depths=beta.RNA.PCA.20bin.ob$depths, index=beta.RNA.PCA.20bin.ob$index, qcut=0.2,slopecut1=0.3,slopecut2=-0.3,doscale=T)


CallT2Dpeak_qvalue_SpeedUP <- function(cellvsPeak.m.aggr=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$cellvsPeak.m.aggr, depths=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$depths, index=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$index,qcut=0.1,slopecut1=0.5,slopecut2=-0.5,doscale=T){
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

#' SelectSigPCs
#'
#' This function is to run regression (linear/logistic) based on the number of PCs, and characteristics of samples.
#' @param BetaPCA, a dataframe, combined PCs and meta data
#' @param Sample sample/donor name
#' @param pheno  phenotype to compare, binary e.g. diseaseStat; continuous e.g. BMI
#' @param is_continuous T or F, if pheno is cpontinuous
#' @return a dataframe: PC_number, PC, PC_pval, ordered by PC_pval
#' @export
#' @examples sigPCs <- SelectSigPCs(BetaPCA, Sample, pheno, is_continuous=F)

SelectSigPCs <- function(BetaPCA, Sample, pheno, is_continuous=F){
    if(is_continuous==F){
        pheno1.cells <- rownames(BetaPCA)[which(BetaPCA[,pheno]==unique(BetaPCA[,pheno])[1])]
        pheno2.cells <- rownames(BetaPCA)[which(BetaPCA[,pheno]==unique(BetaPCA[,pheno])[2])]
        PC_pval <- c()
        for(PC_ind in 1:50){
            PC_pval <- c(PC_pval, wilcox.test(BetaPCA[pheno1.cells, PC_ind], BetaPCA[pheno2.cells, PC_ind])$p.value)
        }
    }
    if(is_continuous==T){
        sample_cell_lis <- list()
        for (d in levels(BetaPCA[,Sample])){
            sample_cell_lis[[d]] <- rownames(BetaPCA)[which(BetaPCA[,Sample]==d)]
        }
        PC_pval <- c()
        sample_pheno.df <- unique(BetaPCA[,c(Sample, pheno)])
        for(PC_ind in 1:50){
            sample_PCmedian_lis <- c()
            for (d in sample_pheno.df[,1]){
                sample_PCmedian_lis <- c(sample_PCmedian_lis, median(BetaPCA[sample_cell_lis[[d]], PC_ind]))
            }
            PC_pval <- c(PC_pval, wilcox.test(sample_pheno.df[, pheno], sample_PCmedian_lis)$p.value)
        }
    }
    PC_pval_df <- data.frame(PC_number=1:50, PC=paste("PC_",1:50,sep=''), pval=PC_pval)
    PC_pval_df.order <- PC_pval_df[order(PC_pval_df$pval, decreasing=F),]
    return(PC_pval_df.order[1:10,])
}


#' CallT2Dpeak_pvalueOneTail
#' @param cellvsPeak.m.aggr
#' @param depths
#' @param index
#' @param doscale
#' @param GlobalSlopes
#' @return
#' @export
#' @examples PCAInfo.20bin.donor.PCA <- CallT2Dpeak_pvalueOneTail(PCAInfo.20bin.donor.ob$cellvsPeak.m.aggr, PCAInfo.20bin.donor.ob$depths, PCAInfo.20bin.donor.ob$index,doscale=T, GlobalSlopes=PCAInfo.20bin.ob.PCA$pseudoregress.all[,1,drop=F])

CallT2Dpeak_pvalueOneTail <- function(cellvsPeak.m.aggr, depths, index, doscale=T, GlobalSlopes){
    cellvsPeak.m.aggr.norm <- t(t(cellvsPeak.m.aggr)/depths)
    cellvsPeak.m.aggr.norm.scale <- apply(cellvsPeak.m.aggr.norm,1,function(x){(x-mean(x))/sd(x)}) %>% t
    # cellvsPeak.m.aggr.scale<- apply(cellvsPeak.m.aggr,1,function(x){(x-mean(x))/sd(x)}) %>% t
    index.scale=scale(index)
    pseudoregress.all<-c()
    #Loop gene by gene
        startTime <- Sys.time()
    for(i in 1:nrow(cellvsPeak.m.aggr)){
        if(any(is.na(as.numeric(cellvsPeak.m.aggr.norm.scale[i,])))){
            pseudoregress.all<-rbind(pseudoregress.all,c(NA,NA,NA,NA,NA))
            next
        }
        if(doscale){
            gmodel1<-glm(as.numeric(cellvsPeak.m.aggr.norm.scale[i,] )~index.scale,family=gaussian) # glm(expr~index), cells from one donor
        }else{
            gmodel1<-glm(as.numeric(cellvsPeak.m.aggr.norm[i,] )~index,family=gaussian)
        }
        OneTail.p <- pt(summary(gmodel1)$coefficients[2,3], gmodel1$df.residual, lower.tail=GlobalSlopes[row.names(cellvsPeak.m.aggr.norm.scale)[i],]<0) # individual slope vs global slope
        cur.slope1 <- summary(gmodel1)$coefficients[,1][2]   # to get slope of pseudoBMIindex
        cur.pvalues1 <- OneTail.p  # to get p value for pseudoBMIindex, which is one tailed against global
        cur.corr <- cor(1e4*as.numeric(cellvsPeak.m.aggr.norm.scale[i,]),index.scale)
        cur.MaxMedian.FC <- max(cellvsPeak.m.aggr.norm[i,])/median(cellvsPeak.m.aggr.norm[i,])
        cur.Max <- 1e6*max(cellvsPeak.m.aggr.norm[i,])
        names(cur.slope1) <- "slope"
        names(cur.pvalues1) <- "pvalue.onetail"
        names(cur.corr) <- "cor"
        names(cur.MaxMedian.FC) <- "MaxMedianFC"
        names(cur.Max) <- "MaxRPKM"
        pseudoregress.all <- rbind(pseudoregress.all,c(cur.slope1,cur.pvalues1,cur.corr,cur.MaxMedian.FC,cur.Max))
    }
        endTime <- Sys.time()
    print(endTime - startTime)

    row.names(pseudoregress.all) <- row.names(cellvsPeak.m.aggr)
    pseudoregress.all <- as.data.frame(pseudoregress.all)
    return(pseudoregress.all)
}

#' scRNA.RePACT optimize speed
#'
#' This function is to run regression based on the number of PCs, and characteristics of samples.
#' @param OBJ, a scRNA-seq Seurat object
#' @param Sample, colnames of donor or sample infomation in OBJ@meta.data
#' @param pheno, the column name of the "characteristics to compare" in OBJ@meta.data, e.g. diseaseStatus or BMI or HBA1C
#' @param pheno_levels, specify the levels of pheno column. e.g. c("HT", "T2D"), specify for binary pheno
#' @param is_continuous, if pheno is continous variable. default is F. e.g. diseaseStatus is F, BMI is T
#' @param if_donorWise, if perform donor wise RePACT, default is F
#' @param binnumber, number of bins for cell grouping
#' @param PCrange. default is "", automatically take top 10 PCs that are significant with phenotype, otherwise specify 10 PCs, e.g. 1:10 or 2:11
#' @param RePACT_qvalCut, qvalue cutoff to determine the significant genes along the disease or other characteristics trajactory, certain cutoff for slope can be further applied
#' @param donorWise_qvalCut, qvalue cutoff to determine the significant genes varying within donors or across donors.
#' @return The function return a list of objects
#' @import Seurat plyr dplyr ggrepel ggplot2 plot3D pscl gridExtra qvalue rlist gplots
#' @export
#' @examples
#' T2D.scRNA.RePACT <- scRNA.RePACT(OBJ=scRNA.OBJ,Sample="Donor", pheno="diseaseStat", is_continuous=F, if_donorWise=F)
scRNA.RePACT <- function(OBJ, Sample, pheno, pheno_levels, is_continuous=F, if_donorWise=F, binnumber=20, PCrange="", RePACT_qvalCut=0.005, donorWise_qvalCut=0.01){
    require(Seurat)
    require(plyr)
    require(dplyr)
    require(ggrepel)
    require(pscl)
    require(qvalue)
    require(rlist)
    GetRePACTLinearmodel.cca<-function(ccaWithinfo=cca.L2.info, prefix="CC",pheno="Disease",CCrange=1:10){
        CCnames <- paste("ccaWithinfo$",prefix,"_",CCrange, sep = "")
        ccaWithinfo[,pheno] <- as.factor(ccaWithinfo[,pheno])
        ccaWithinfo[,pheno] <- as.numeric(levels(ccaWithinfo[,pheno]))[ccaWithinfo[,pheno]]
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
    # BetaPCA contains 50PCs and OBJ meta data
    BetaPCA <- OBJ@reductions$pca@cell.embeddings %>% Tomerge_v2(.,OBJ@meta.data)
    BetaPCA <- BetaPCA[,1:50] %>% apply(.,2,function(x){x/sqrt(sum(x^2))}) %>% cbind(.,BetaPCA[,51:ncol(BetaPCA)])
    BetaPCA[,Sample] <- factor(BetaPCA[,Sample] ,levels=unique(BetaPCA[,Sample]))
    if(is_continuous==F){
	BetaPCA[,pheno] <- as.character(BetaPCA[,pheno])
#	BetaPCA[,pheno] <- factor(BetaPCA[,pheno] ,levels=sort(unique(BetaPCA[,pheno])))	
	BetaPCA[,pheno] <- factor(BetaPCA[,pheno] ,levels=pheno_levels)
    }
    beta.rna.pca.withinfo.subs<-list()
    print("Calculating cell pseudo-index for 100 times")
    # Loop 100 times
    for(i in 1:100){
        beta.rna.pca.withinfo.sub<-c()
        # Loop through each sample/donor, sample 200 cells from each sample/donor, better require samples having more than 200 cells to have better significance
        for (d in levels(BetaPCA[,Sample])){
            tmp<- BetaPCA[which(BetaPCA[,Sample]==d),]
            if(nrow(tmp)<=200){
                beta.rna.pca.withinfo.sub <- rbind(beta.rna.pca.withinfo.sub,tmp)
            }else{
                beta.rna.pca.withinfo.sub <- rbind(beta.rna.pca.withinfo.sub,tmp[sample(1:nrow(tmp),200),])
            }
        }
        beta.rna.pca.withinfo.subs <- c(beta.rna.pca.withinfo.subs,list(beta.rna.pca.withinfo.sub))
    }
    pseudo.indexes<-list()
    if(PCrange==""){
        sigPCs <- SelectSigPCs(BetaPCA, Sample, pheno, is_continuous=F)
        sigPCs_number <- as.numeric(sigPCs[,1])
    }else{
        sigPCs_number <- PCrange
    }
    # for each loop
    for(i in 1:100){
        # build regression model: Disease/BMI ~ (PC1+PC2+...PC10), consider only including significant PCs in the future
        if(is_continuous==F){
            md <- GetRePACTmodel.cca(ccaWithinfo=beta.rna.pca.withinfo.subs[[i]],prefix="PC",pheno=pheno,CCrange=c(sigPCs_number))
        }else{
            md <- GetRePACTLinearmodel.cca(ccaWithinfo=beta.rna.pca.withinfo.subs[[i]],prefix="PC",pheno=pheno,CCrange=c(sigPCs_number))
        }
        trainingdata <- beta.rna.pca.withinfo.subs[[i]] # sampled 200 cells datainfo
        Restdata <- BetaPCA[setdiff(row.names(BetaPCA),row.names(beta.rna.pca.withinfo.subs[[i]])),] # the rest of datainfo, excluding sampled 200 cells datainfo
        Alldata <- BetaPCA
        beta.rna.pca.withinfo.subs[[i]]$pseudo.index <- apply(trainingdata[,sigPCs_number],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])}) # fitted disease/BMI index: intercept + sum(coef1*PC1, coef2*PC2, coef3*PC3,...,coef10*PC10)
        Restdata$pseudo.index <- apply(Restdata[,sigPCs_number],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})
        pseudo.indexes <- c(pseudo.indexes,list(apply(Alldata[,sigPCs_number],1,function(x){md$coefficients[[1]]+sum(x*md$coefficients[2:11])})))
    } # every cell has 100 pseudo-index due to 100 times sampling
    BetaPCA$pseudo.index.balanced <- do.call(cbind,pseudo.indexes) %>% rowMeans() # Get one average pseudo-index from 100 indexes for every cell
    #adjustrange = seq(0, 180, length.out = 13)
    # compare pseudo-index distribution between two categorical variables, e.g. disease vs non-disease,
    #ks.test(BetaPCA[which(BetaPCA[,pheno]==unique(BetaPCA[,pheno])[1]),]$pseudo.index.balanced,BetaPCA[which(BetaPCA[,pheno]==unique(BetaPCA[,pheno])[2]),]$pseudo.index.balanced)
    BetaPCA$rank <- rank(BetaPCA$pseudo.index.balanced) # rank cells based on their pseudo-index
    # bin the cells based on the pseudo-index rank, coonsidering the similar depth and cell number with each bin
    print("binning the cells along the trajectory")
    beta.RNA.PCA.20bin.ob <- MakeEvenBinBydepth_SpeedUP(OBJ=OBJ, data.info=BetaPCA[,51:ncol(BetaPCA)], binnumber=binnumber)
    # measure the gene trend across two conditions or a continuous progression by regressing pseudo-index~expression. Get slope and pvalue/qvalue as stats
    print("Getting statistics")
#    betaT2D.diffGene.20bin.PCA <- CallT2Dpeak_qvalue(beta.RNA.PCA.20bin.ob$cellvsPeak.m.aggr, beta.RNA.PCA.20bin.ob$depths, beta.RNA.PCA.20bin.ob$index, qcut=0.2,slopecut1=0.3,slopecut2=-0.3,doscale=T)
    betaT2D.diffGene.20bin.PCA <- CallT2Dpeak_qvalue_SpeedUP(beta.RNA.PCA.20bin.ob$cellvsPeak.m.aggr, beta.RNA.PCA.20bin.ob$depths, beta.RNA.PCA.20bin.ob$index, qcut=0.2,slopecut1=0.3,slopecut2=-0.3,doscale=T)

    bin20.g.up <- row.names(subset(betaT2D.diffGene.20bin.PCA$UP,qvalue<RePACT_qvalCut & MaxRPKM>quantile(betaT2D.diffGene.20bin.PCA$pseudoregress.all$MaxRPKM,0.25)))
    bin20.g.dn <- row.names(subset(betaT2D.diffGene.20bin.PCA$DN,qvalue<RePACT_qvalCut & MaxRPKM>quantile(betaT2D.diffGene.20bin.PCA$pseudoregress.all$MaxRPKM,0.25)))
    RePACT_call <- list(UP=bin20.g.up, DN=bin20.g.dn)
    # run RePACT within each sample/donor to explore intra-donor heterogeneity
    if(if_donorWise==T){
        print("Start donor-wise RePACT")
        PCAInfo <- BetaPCA
        RePACT.diff.genes.bydonor<-list()
        # for each donor
        for (curDonor in levels(PCAInfo[,Sample])){
            print(paste("processing",which(levels(PCAInfo[,Sample])==curDonor), "/",length(levels(PCAInfo[,Sample]))))
            PCAInfo.donor <- subset(PCAInfo, Sample==curDonor) # take out the cells from the donor
            PCAInfo.donor$rank <- rank(PCAInfo.donor$pseudo.index.balanced) # rank the cells based on the globally assigned pseudo-index
            OBJ.donor <- subset(OBJ, cells=rownames(PCAInfo.donor))
            PCAInfo.20bin.donor.ob <- MakeEvenBinBydepth_SpeedUP(OBJ=OBJ.donor, data.info=PCAInfo.donor[,51:ncol(PCAInfo.donor)], binnumber=binnumber)
            PCAInfo.20bin.donor.PCA <- CallT2Dpeak_pvalueOneTail(PCAInfo.20bin.donor.ob$cellvsPeak.m.aggr, PCAInfo.20bin.donor.ob$depths, PCAInfo.20bin.donor.ob$index,doscale=T, GlobalSlopes=betaT2D.diffGene.20bin.PCA$pseudoregress.all[,1,drop=F])
            RePACT.diff.genes.bydonor <- c(RePACT.diff.genes.bydonor,list(PCAInfo.20bin.donor.PCA))
        }
        print("Combining donor-wise results")
        names(RePACT.diff.genes.bydonor) <- levels(PCAInfo[,Sample])
        ps <- lapply(RePACT.diff.genes.bydonor,function(x){x$pvalue.onetail}) %>% do.call(cbind,.) # combine pvalue, gene~donor 
        row.names(ps) <- row.names(RePACT.diff.genes.bydonor[[1]])
        ps[is.na(ps)] <- 1
        #slopes <- lapply(RePACT.diff.genes.bydonor,function(x){x$pseudoregress.all$slope}) %>% do.call(cbind,.)
	#row.names(slopes) <- row.names(RePACT.diff.genes.bydonor[[1]]$pseudoregress.all)
	slopes <- lapply(RePACT.diff.genes.bydonor,function(x){x$slope}) %>% do.call(cbind,.) # combine slope, gene~donor
        row.names(slopes) <- row.names(RePACT.diff.genes.bydonor[[1]])
        FishersMethod.p <- apply(ps,1,function(x){-2*sum(log(x))}) %>% pchisq(.,2*11,lower.tail=FALSE)  ## k=11, df=2k
        FishersMethod.p <- FishersMethod.p[c(bin20.g.up,bin20.g.dn)] # subset RePACT genes
        FishersMethod.q <- qvalue(FishersMethod.p)$qvalues %>% .[order(.)] # calculate qvalues
        FishersMethod.q.df <- data.frame(gene=names(FishersMethod.q),qvalueInSample=FishersMethod.q) # <gene> <qvalue>
        Globaltag<-c()
        Globaltag[which(FishersMethod.q.df$gene %in% bin20.g.up)]<-"UP"
        Globaltag[which(FishersMethod.q.df$gene %in% bin20.g.dn)]<-"DN"
        Globaltag[which(!FishersMethod.q.df$gene %in% c(bin20.g.up,bin20.g.dn))]<-""
        FishersMethod.q.df$Globaltag <- Globaltag

#        GeneLabel <- c(subset(FishersMethod.q.df,Globaltag=="UP") %>% .[order(.$qvalueInSample),] %>% head(.,n=25) %>% row.names,subset(FishersMethod.q.df,Globaltag=="DN") %>% .[order(.$qvalueInSample),] %>% head(.,n=25) %>% row.names)
#        FishersMethod.q.df$GeneLabel <- ifelse(FishersMethod.q.df$gene %in% GeneLabel,row.names(FishersMethod.q.df),"")
#        FishersMethod.q.df <- Tomerge_v2(FishersMethod.q.df,PCAInfo.20bin.ob.PCA$pseudoregress.all[,"qvalue",drop=F])
        FishersMethod.q.df$Globaltag <- as.factor(FishersMethod.q.df$Globaltag)
        FishersMethod.q.dn.df <- subset(FishersMethod.q.df,Globaltag=="DN")
        FishersMethod.q.up.df <- subset(FishersMethod.q.df,Globaltag=="UP")
        FishersMethod.q.dn.df$rank <- rank(FishersMethod.q.dn.df$qvalueInSample)
        FishersMethod.q.up.df$rank <- rank(FishersMethod.q.up.df$qvalueInSample)
        FishersMethod.q.dn.df$InSampleTag <- ifelse(FishersMethod.q.dn.df$qvalueInSample<donorWise_qvalCut,"Intra-donor","Inter-donor")
        FishersMethod.q.up.df$InSampleTag <- ifelse(FishersMethod.q.up.df$qvalueInSample<donorWise_qvalCut,"Intra-donor","Inter-donor")
        DN.hetero <- names(FishersMethod.q[bin20.g.dn])[FishersMethod.q[bin20.g.dn]<=donorWise_qvalCut]
        DN.homo <- names(FishersMethod.q[bin20.g.dn])[FishersMethod.q[bin20.g.dn]>donorWise_qvalCut]
        UP.hetero <- names(FishersMethod.q[bin20.g.up])[FishersMethod.q[bin20.g.up]<=donorWise_qvalCut]
        UP.homo <- names(FishersMethod.q[bin20.g.up])[FishersMethod.q[bin20.g.up]>donorWise_qvalCut]
        RePACT_donorWise_intermediate <- list(FishersMethod.q.df=FishersMethod.q.df, FishersMethod.q.up.df=FishersMethod.q.up.df, FishersMethod.q.dn.df=FishersMethod.q.dn.df)
        RePACT_donorWise_call <- list(DN.interdonor=DN.hetero, DN.intradonor=DN.homo, UP.interdonor=UP.hetero, UP.intradonor=UP.homo)
        return(list(BetaPCA=BetaPCA,beta.RNA.PCA.20bin.ob=beta.RNA.PCA.20bin.ob, betaT2D.diffGene.20bin.PCA=betaT2D.diffGene.20bin.PCA, RePACT_call=RePACT_call, RePACT_donorWise_intermediate=RePACT_donorWise_intermediate, RePACT_donorWise_call=RePACT_donorWise_call))
    }else{
        return(list(BetaPCA=BetaPCA,beta.RNA.PCA.20bin.ob=beta.RNA.PCA.20bin.ob, betaT2D.diffGene.20bin.PCA=betaT2D.diffGene.20bin.PCA, RePACT_call=RePACT_call))
    }
}
#' GetRePACTmodel.cca
#'
#' This function is to run logistic regression based on the number of LSIs, and characteristics of samples.
#' @param ccaWithinfo, phenotable
#' @param PCrange, Specified PCs for regression
#' @param prefix
#' @param pheno, the column name of the "characteristics to compare" in the phenodic.use dataframe
#' @param CCrange
#' @return The function return logistic model
#' @import Seurat pscl
#' @export
#' @examples
#'
GetRePACTmodel.cca<-function(ccaWithinfo=cca.L2.info, prefix="CC",pheno="Disease",CCrange=1:10){
  require(pscl)
  CCnames <- paste("ccaWithinfo$",prefix,"_",CCrange, sep = "")
  ccaWithinfo[,pheno]<-as.factor(ccaWithinfo[,pheno])
  form <- formula(paste("ccaWithinfo[,pheno]", paste(CCnames, collapse = "+"), sep = "~"))
  model <- glm(form, family = "binomial")
# p <- ggplot(ccaWithinfo) + aes_string("CC_1", "CC_2", color = pheno) + geom_point() + scale_color_brewer(palette = "Set2") + geom_abline(slope = model$coefficients[3]/model$coefficients[2])
# model.para <- summary(model)
# lrt.pvalue <- 1 - pchisq(summary(model)$null.deviance - summary(model)$deviance, summary(model)$df.null - summary(model)$df.residual)
  return(model)
}

#' standard_fun
#'
#' This function is to run basic processing for Seurat OBJ
#' @param OBJ, Seurat OBJ
#' @return The function return a Seurat OBJ
#' @import Seurat plyr dplyr ggrepel ggplot2 plot3D pscl gridExtra qvalue rlist gplots
#' @export
#' @examples

standard_fun <- function(OBJ){
        require(Seurat)
        OBJ <- NormalizeData(OBJ)
        OBJ <- FindVariableFeatures(OBJ, selection.method = "vst", nfeatures = 2000)
        OBJ <- ScaleData(OBJ, features = rownames(OBJ))
        OBJ <- RunPCA(OBJ, features = VariableFeatures(object = OBJ))
        DefaultAssay(OBJ) <- "RNA"
	OBJ <- RunUMAP(OBJ, dims = 1:10)
        return(OBJ)
}

#' plot_scRNA_RePACT_heatmap
#'
#' This function is to take a vetor of genes and plot heatmaps from RePACT result
#' @param geneset, a vector of genes
#' @param RePACT_OBJ, RePACT reesult OBJ
#' @return The function return a plot
#' @import Seurat plyr reshape2 dplyr ggrepel ggplot2 plot3D pscl gridExtra qvalue rlist gplots
#' @export
#' @examples

plot_scRNA_RePACT_heatmap <- function(geneset, RePACT_OBJ, ifShowGene=T){
   require(ggrepel)
   require(ggplot2)
   require(plot3D)
   require(gridExtra)
   require(gplots)
   require(reshape2)
    normalize_01<-function(vector){
            normed<-(vector-min(vector))/(max(vector)-min(vector))
            return(normed)
    }
    T2D.bindata <- RePACT_OBJ$betaT2D.diffGene.20bin.PCA$UPDN.toplot[, intersect(geneset, colnames(RePACT_OBJ$betaT2D.diffGene.20bin.PCA$UPDN.toplot))]
    T2D.bindata <- data.frame(apply(T2D.bindata[, -ncol(T2D.bindata)-1, drop = F],2, normalize_01), bin = 1:nrow(T2D.bindata))
    T2D.bindata.m <- reshape2::melt(T2D.bindata, id.vars = "bin")
    T2D.bindata.m$bin <- factor(T2D.bindata.m$bin,levels=1:nrow(T2D.bindata))
    T2D.bindata.m$variable <- factor(T2D.bindata.m$variable, levels=geneset)
    T2D.bindata.m <- T2D.bindata.m[which(!is.na(T2D.bindata.m$variable)),]
    colnames(T2D.bindata.m)[3] <- "normExpr"
    if(ifShowGene==T){
	p1 <- ggplot(T2D.bindata.m) + aes(bin, variable, fill = normExpr) + geom_tile() + scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 0.6) +
            labs(y = "Variable genes")
    }else{
	p1 <- ggplot(T2D.bindata.m) + aes(bin, variable, fill = normExpr) + geom_tile() + scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 0.6) + 
             labs(y = "Variable genes")+theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
             theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ theme(legend.position = "none")
    }
    return(p1)
}



