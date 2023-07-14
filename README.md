# RePACT 
## (Regressing Principle components for the Assembly of Continuous Trajectory)
We developed RePACT as a general method to sensitively identify disease relevant gene signatures / ATAC-peaks using single cell data.The goal of RePACT is to find the best trajectory to rank single cells (e.g. in our case, the β cells) reflecting the change of disease status. <br /> 
It can identify:<br /> 
- Variable genes/ATAC-peaks across cells from multiple samples
- Intra-donor, and inter-donor variable genes/ATAC-peaks 

Please cite: <br /> 
- [Single-Cell Heterogeneity Analysis and CRISPR Screen Identify Key β-Cell-Specific Disease Genes. Zhou & Chen et al,2019](https://doi.org/10.1016/j.celrep.2019.02.043) 
- Single cell multi-omics analysis reveals diabetes-associated β-cell heterogeneity in human islets (under revision)

## Install
```
library(devtools)
devtools::install_github("chenweng1991/RePACT") #install RePACT
library("RePACT")
```

## Prepare a Seurat object for RePACT
The Seurat object should only contain cells from one cell type, and rerun PCA/LSI after subsetting.<br /> 
An example:
```
# scRNA
OBJ <- NormalizeData(OBJ)
OBJ <- FindVariableFeatures(OBJ, selection.method = "vst", nfeatures = 2000)
OBJ <- ScaleData(OBJ, features = rownames(OBJ))
OBJ <- RunPCA(OBJ, features = VariableFeatures(object = OBJ))
DefaultAssay(OBJ) <- "RNA"

# snATAC
OBJ <- RunTFIDF(OBJ)
OBJ <- FindTopFeatures(OBJ, min.cutoff = "q5")
OBJ <- RunSVD(OBJ)
DefaultAssay(OBJ) <- "ATAC"
```
The OBJ@meta.data should contain a "Sample" column and a column to compare using RePACT, e.g. diseaseStat (binary, "HT"/"T2D")
```
# scRNA meta example (cell name omitted):
   nCount_RNA nFeature_RNA percent.mt Sample diseaseStat
1        1079          573          0     H5           H
2        2481         1209          0     H6           H
3        5265         2315          0   T2D2         T2D
4        4202         2000          0     H3           H
5        1462          987          0     H3           H

# snATAC meta example (cell name omitted):
   nCount_ATAC nFeature_ATAC Fragments PeakRatio Sample diseaseStat
1         1093          1093      2190 0.4059361   T2D3         T2D
2          865           865      1555 0.4398714    HT8          HT
3         2153          2153      3838 0.4317353    HT3          HT
4         3456          3456      6628 0.3919734   T2D5         T2D
5          788           788      1196 0.5610368    HT8          HT
```

## RePACT on scRNA-seq
Examples to run RePACT on the downloaded Seurat object.
```
scRNA.OBJ <- readRDS("data/Old9.Beta.scRNA.rds") # Zhou & Chen et al,2019
T2D.scRNA.RePACT <- scRNA.RePACT(OBJ=scRNA.OBJ,Sample="Sample", pheno="diseaseStat", is_continuous=F, if_donorWise=T)
```
Output
```
T2D.scRNA.RePACT$BetaPCA: cell ~ (PC + metadata + pseudoindex)
T2D.scRNA.RePACT$RePACT_call: UP and DN genes
T2D.scRNA.RePACT$beta.RNA.PCA.20bin.ob: intermediate files
T2D.scRNA.RePACT$RePACT_donorWise_intermediate: fisher exact test applied to UP and DN genes within or across donors/samples
T2D.scRNA.RePACT$betaT2D.diffGene.20bin.PCA: 20bins on trajectory ~ genes
T2D.scRNA.RePACT$RePACT_donorWise_call: significant genes changed within or across samples
```
Plot RePACT genes
```
# 3D scatterplot on three specified PCs
scatter3D(T2D.scRNA.RePACT$BetaPCA[,"PC_1"], T2D.scRNA.RePACT$BetaPCA[,"PC_2"],T2D.scRNA.RePACT$BetaPCA[,"PC_3"],ticktype = "detailed", pch = 20, theta = 150, phi = 180, colvar = ifelse(T2D.scRNA.RePACT$BetaPCA[,"diseaseStat"]==unique(T2D.scRNA.RePACT$BetaPCA[,"diseaseStat"])[1],1,0), bty = "b2", cex = 0.6, col = alpha.col(col = c("blue", "brown"), 0.6))
# pseudoindex distribution violin for each donor
ggplot(T2D.scRNA.RePACT$BetaPCA)+aes(Sample,pseudo.index.balanced,fill=diseaseStat)+geom_violin()+geom_boxplot(width=0.2,outlier.shape = NA,notch=F,coef = 0,fill="grey25",color="grey75")+coord_flip()+theme_classic()+scale_fill_manual(values=c("blue", "brown"))+theme_bw()+theme(legend.position="none")
# pseudoindex distribution from healthy and T2D cells
ggplot(T2D.scRNA.RePACT$BetaPCA)+aes(pseudo.index.balanced,fill=diseaseStat)+geom_density(alpha=0.75)+scale_fill_manual(values=c("blue","brown"))+theme_classic()+theme(legend.position="none")
# heatmap to show gene trend along T2D trajectory
geneset <- c(T2D.scRNA.RePACT$RePACT_call$UP[1:30],T2D.scRNA.RePACT$RePACT_call$DN[1:30]) # select top 30 genes from UP and DN
T2D.bindata <- T2D.scRNA.RePACT$betaT2D.diffGene.20bin.PCA$UPDN.toplot[, geneset]
T2D.bindata <- data.frame(apply(T2D.bindata[, -ncol(T2D.bindata), drop = F],2, normalize_01), bin = 1:nrow(T2D.bindata))
T2D.bindata.m <- reshape2::melt(T2D.bindata, id.vars = "bin")
T2D.bindata.m$bin <- factor(T2D.bindata.m$bin,levels=1:nrow(T2D.bindata))
T2D.bindata.m$variable <- factor(T2D.bindata.m$variable, levels=geneset)
ggplot(T2D.bindata.m) + aes(bin, variable, fill = value) + geom_tile() + scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 0.6) + labs(y = "Variable genes")
```
<img src="https://raw.githubusercontent.com/chenweng1991/RePACT/RePACT.organized/image/RePACT.violinheat.png" alt="drawing" width="600"/>

Plot donorWise RePACT genes
```
# ranked scatterplot to show intra-donor DN genes
ggplot(T2D.scRNA.RePACT$RePACT_donorWise_intermediate$FishersMethod.q.dn.df)+aes(rank,-log10(qvalueInSample),label=GeneLabel,shape=InSampleTag)+geom_point(color="red")+geom_text_repel(color="red")+geom_hline(yintercept=2,linetype=2)+theme_classic()+scale_shape_manual(values=c(20,2))+ggtitle("DN")
# ranked scatterplot to show intra-donor UP genes
ggplot(T2D.scRNA.RePACT$RePACT_donorWise_intermediate$FishersMethod.q.up.df)+aes(rank,-log10(qvalueInSample),label=GeneLabel,shape=InSampleTag)+geom_point(color="darkgreen")+geom_text_repel(color="darkgreen")+geom_hline(yintercept=2,linetype=2)+theme_classic()+scale_shape_manual(values=c(20,2))+ggtitle("UP")
# ranked scatterplot to show intra-donor genes, global_qvalue ~ intradonor_qvalue
ggplot(T2D.scRNA.RePACT$RePACT_donorWise_intermediate$FishersMethod.q.df)+aes(-log10(qvalue),-log10(qvalueInSample),color=Globaltag,label=GeneLabel)+geom_point(size=0.5)+geom_text_repel(aes(color=Globaltag))+scale_color_manual(values=c("red","darkgreen"))+geom_hline(yintercept=2,linetype=2)+geom_vline(xintercept=-log10(0.005),linetype=2)+theme_classic()+xlab("-log10(Global qvalue)")
```

![](https://github.com/chenweng1991/RePACT/blob/master/image/scRNA.RePACT.donorWise.PNG)

## RePACT on snATAC-seq



## RePACT on snATAC-seq to identify intra-dnor and inter-donor ATAC-peaks
```
OBJ <- readRDS("data/Beta.snATAC.rds") 
T2D.snATAC.RePACT <- snATAC.RePACT(OBJ=OBJ, Sample="Sample", pheno="diseaseStat", is_continuous=F, if_donorWise=T, RePACT_qvalCut=0.01, donorWise_qvalCut=0.01)
```
Output
```
T2D.snATAC.RePACT$LSIInfo: cell ~ (LSI + metadata + pseudoindex)
T2D.snATAC.RePACT$LSIInfo.20bin.ob: intermediate files for binning the cells 
T2D.snATAC.RePACT$LSIInfo.20bin.ob.LSI: RePACT peaks 
T2D.snATAC.RePACT$RePACT_donorWise_intermediate: fisher exact test applied to UP and DN peaks within or across donors/samples
T2D.snATAC.RePACT$RePACT_donorWise_call: donorwise RePACT peaks
```
Plot RePACT peaks
```
## Plot 3D scatterplot LSI_1 LSI_2 LSI_3
scatter3D(T2D.snATAC.RePACT$LSIInfo[,"LSI_1"], T2D.snATAC.RePACT$LSIInfo[,"LSI_2"], T2D.snATAC.RePACT$LSIInfo[,"LSI_3"],ticktype = "detailed", pch = 20, theta = 90, phi = 30, colvar = ifelse(T2D.snATAC.RePACT$LSIInfo[,'diseaseStat']=='H',1,0), bty = "b2", cex = 0.3, col = alpha.col(col = c("steelblue", "red"), 0.6))
## Plot violin pseudoindex
ggplot(T2D.snATAC.RePACT$LSIInfo)+aes(Sample,pseudo.index.balanced,fill=diseaseStat)+geom_violin()+geom_boxplot(width=0.2,outlier.shape = NA,notch=F,coef = 0,fill="grey25",color="grey75")+coord_flip()+theme_classic()+scale_fill_manual(values=c("steelblue", "red"))+theme_bw()+theme(legend.position="none")
## Plot pseudoindex density between HT and T2D
ggplot(T2D.snATAC.RePACT$LSIInfo)+aes(pseudo.index.balanced,fill=diseaseStat)+geom_density(alpha=0.75)+scale_fill_manual(values=c("steelblue","red"))+theme_classic()+theme(legend.position="none")
```
![](https://github.com/chenweng1991/RePACT/blob/master/image/T2D.snATAC.RePACT.PNG)
Plot donorwise RePACT peaks
```
# Plot intra-donor dn peaks
ggplot(T2D.snATAC.RePACT$RePACT_donorWise_intermediate$FishersMethod.q.dn.df)+aes(rank,-log10(qvalueInSample),label=tag2,color=tag1)+geom_point()+geom_text_repel()+geom_hline(yintercept=2)+theme_classic()+ggtitle("DN peaks")
# Plot intra-donor up peaks
ggplot(T2D.snATAC.RePACT$RePACT_donorWise_intermediate$FishersMethod.q.up.df)+aes(rank,-log10(qvalueInSample),label=tag2,color=tag1)+geom_point()+geom_text_repel()+geom_hline(yintercept=2)+theme_classic()+ggtitle("UP peaks")
```
