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
## scRNA RePACT
```
scRNA.OBJ <- readRDS("Old9.Beta.scRNA.rds") # load your data
T2D.scRNA.RePACT <- scRNA.RePACT(scRNA.OBJ,
                                 Sample="Sample",
                                 pheno="diseaseStat",
                                 pheno_levels=c("H","T2D"),
                                 is_continuous=F,
                                 if_donorWise=F,
                                 binnumber=20, PCrange="", RePACT_qvalCut=0.005, donorWise_qvalCut=0.01)
```
scRNA RePACT PLOTS:
3D scatterplot
```
scatter3D(T2D.scRNA.RePACT$BetaPCA[,"PC_1"], T2D.scRNA.RePACT$BetaPCA[,"PC_2"],T2D.scRNA.RePACT$BetaPCA[,"PC_3"],
          ticktype = "detailed", pch = 20, theta = 150, phi = 180, 
          colvar = ifelse(T2D.scRNA.RePACT$BetaPCA[,"diseaseStat"]==unique(T2D.scRNA.RePACT$BetaPCA[,"diseaseStat"])[1],1,0),
          bty = "b2", cex = 0.6, col = alpha.col(col = c("blue", "brown"), 0.6))
```
Pseudoindex distribution violin for each donor
```
ggplot(T2D.scRNA.RePACT$BetaPCA %>% .[complete.cases(.),])+aes(Sample,pseudo.index.balanced,fill=diseaseStat)+
       geom_violin()+
       geom_boxplot(width=0.2,outlier.shape = NA,notch=F,coef = 0,fill="grey25",color="grey75")+
       coord_flip()+theme_classic()+scale_fill_manual(values=c("blue", "brown"))+
       theme_bw()+theme(legend.position="none")
```
pseudoindex distribution from healthy and T2D cells
```
ggplot(T2D.scRNA.RePACT$BetaPCA)+aes(pseudo.index.balanced,fill=diseaseStat)+geom_density(alpha=0.75)+
       scale_fill_manual(values=c("blue","brown"))+theme_classic()+theme(legend.position="none")
```
heatmap to show gene trend along T2D trajectory
```
plot_scRNA_RePACT_heatmap(geneset=c(T2D.scRNA.RePACT$RePACT_call$UP[1:30], T2D.scRNA.RePACT$RePACT_call$DN[1:30]),
                          RePACT_OBJ=T2D.scRNA.RePACT, ifShowGene=T)
```
DonorWise scRNA RePACT to identify intra-donor heterogeneity
```
T2D.scRNA.RePACT <- scRNA.RePACT(OBJ=scRNA.OBJ,
                                 Sample="Sample",
                                 pheno="diseaseStat",
                                 pheno_levels=c("H","T2D"),
                                 is_continuous=F,
                                 if_donorWise=T,
                                 binnumber=20, PCrange="", RePACT_qvalCut=0.005, donorWise_qvalCut=0.01)

lengths(T2D.scRNA.RePACT$RePACT_donorWise_call)
```
## snATAC RePACT
```
snATAC.OBJ <- readRDS("Beta.snATAC.rds") # load your data
T2D.snATAC.RePACT <- snATAC.RePACT(OBJ=snATAC.OBJ,
                                   Sample="Sample",
                                   pheno="diseaseStat",
                                   pheno_levels=c("HT", "T2D"),
                                   is_continuous=F,
                                   if_donorWise=F,
                                   binnumber=20, LSIrange="", RePACT_qvalCut=0.01, donorWise_qvalCut=0.01)
```
Plot 3D scatterplot LSI_1 LSI_2 LSI_3
```
scatter3D(T2D.snATAC.RePACT$LSIInfo[,"LSI_1"], T2D.snATAC.RePACT$LSIInfo[,"LSI_2"], T2D.snATAC.RePACT$LSIInfo[,"LSI_3"],
          ticktype = "detailed", pch = 20, theta = 90, phi = 30, colvar = ifelse(T2D.snATAC.RePACT$LSIInfo[,'diseaseStat']=='HT',1,0),
          bty = "b2", cex = 0.3, col = alpha.col(col = c("steelblue", "red"), 0.6))
```
Plot violin pseudoindex
```
ggplot(T2D.snATAC.RePACT$LSIInfo)+aes(Sample,pseudo.index.balanced,fill=diseaseStat)+
       geom_violin()+geom_boxplot(width=0.2,outlier.shape = NA,notch=F,coef = 0,fill="grey25",color="grey75")+
       coord_flip()+theme_classic()+scale_fill_manual(values=c("steelblue", "red"))+
       theme_bw()+theme(legend.position="none")
```
Plot pseudoindex density between HT and T2D
```
ggplot(T2D.snATAC.RePACT$LSIInfo)+aes(pseudo.index.balanced,fill=diseaseStat)+
      geom_density(alpha=0.75)+scale_fill_manual(values=c("steelblue","red"))+
      theme_classic()+theme(legend.position="none")
```
heatmap to show peak trend along T2D trajectory
```
plot_scATAC_RePACT_heatmap(peakSet=c(T2D.snATAC.RePACT$RePACT_call$UP, T2D.snATAC.RePACT$RePACT_call$DN), T2D.snATAC.RePACT)
```
DonorWise snATAC RePACT to identify intra-donor heterogeneity
```
T2D.snATAC.RePACT <- snATAC.RePACT(OBJ=snATAC.OBJ,
                                   Sample="Sample",
                                   pheno="diseaseStat",
                                   pheno_levels=c("HT", "T2D"),
                                   is_continuous=F,
                                   if_donorWise=T,
                                   binnumber=20, LSIrange="", RePACT_qvalCut=0.01, donorWise_qvalCut=0.01)

lengths(T2D.snATAC.RePACT$RePACT_donorWise_call)
```
