# RePACT 
## (Regressing Principle components for the Assembly of Continuous Trajectory)
We developed RePACT as a general method to sensitively identify disease relevant gene signatures / ATAC-peaks using single cell data.The goal of RePACT is to find the best trajectory to rank single cells (e.g. in our case, the β cells) reflecting the change of disease status. <br /> 
It can identify:<br /> 
- Variable genes/ATAC-peaks across cells from multiple samples
- Intra-donor, and inter-donor variable genes/ATAC-peaks 

Please cite: <br /> 
- [Single-Cell Heterogeneity Analysis and CRISPR Screen Identify Key β-Cell-Specific Disease Genes. Zhou & Chen et al,2019](https://doi.org/10.1016/j.celrep.2019.02.043) 
- Single cell multi-omics analysis reveals diabetes-associated β-cell heterogeneity in human islets (In preparation)

## Install
```
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
OBJ <- ScaleData(OBJ, features = rownames(pbmc))
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
scRNA.OBJ <- readRDS("../Old9.Beta.scRNA.rds")
T2D.scRNA.RePACT <- repact.scRNA(scRNA.OBJ, PCrange=1:20,pheno="diseaseStat",is_continuous="F",norm_index=T,SlopeCut=0.05,output_name2="T2D_Beta.scRNA.RePACT")

p1 <- ggplot(T2D.scRNA.RePACT$RepACT.obj$PCanfpheno %>% .[complete.cases(.),])+aes_string("Sample","pseudo.index",group="Sample",fill=pheno)+geom_violin()+coord_flip()+scale_fill_manual(values=c("blue", "red"))+theme_classic()
p2 <- Do_heatmap(T2D.scRNA.RePACT$RepACT.2nd.ob$bin.data,df1=T2D.scRNA.RePACT$RepACT.2nd.ob$BINlinear.result.summarized$UP,df2=T2D.scRNA.RePACT$RepACT.2nd.ob$BINlinear.result.summarized$DOWN,rankname="rank",top_gene_num=20)
grid.arrange(p1,p2,ncol=2)
```

The function will return a list for plots and downstream analysis, also the output files will be written out to the disk.
```
T2D_Beta.scRNA.RePACTDOWN_gene.csv
T2D_Beta.scRNA.RePACT.genes.FCs.csv
T2D_Beta.scRNA.RePACT.genes.Full_list.csv
T2D_Beta.scRNA.RePACTUP_gene.csv
```
<img src="https://raw.githubusercontent.com/chenweng1991/RePACT/RePACT.organized/image/RePACT.violinheat.png" alt="drawing" width="600"/>

## RePACT on scRNA-seq to identify intra-dnor and inter-donor gene signitures 
```
scRNA.OBJ <- readRDS("../Old9.Beta.scRNA.rds")
donorWise.T2D.RePACT <- DonorWise.scRNA.RePACT(OBJ=scRNA.OBJ, pheno="diseaseStat")
# Take output List for key plots
pdf("T2D.scRNA.RePACT.donorWise.pdf",18,5)
p1 <- ggplot(donorWise.T2D.RePACT$FishersMethod.q.dn.df)+aes(rank,-log10(qvalueInSample),label=GeneLabel,shape=InSampleTag)+geom_point(color="red")+geom_text_repel(color="red")+geom_hline(yintercept=2,linetype=2)+theme_classic()+scale_shape_manual(values=c(20,2))+ggtitle("DN")
p2 <- ggplot(donorWise.T2D.RePACT$FishersMethod.q.up.df)+aes(rank,-log10(qvalueInSample),label=GeneLabel,shape=InSampleTag)+geom_point(color="darkgreen")+geom_text_repel(color="darkgreen")+geom_hline(yintercept=2,linetype=2)+theme_classic()+scale_shape_manual(values=c(20,2))+ggtitle("UP")
p3 <- ggplot(donorWise.T2D.RePACT$FishersMethod.q.df)+aes(-log10(qvalue),-log10(qvalueInSample),color=Globaltag,label=GeneLabel)+geom_point(size=0.5)+geom_text_repel(aes(color=Globaltag))+scale_color_manual(values=c("red","darkgreen"))+geom_hline(yintercept=2,linetype=2)+geom_vline(xintercept=-log10(0.005),linetype=2)+theme_classic()+xlab("-log10(Global qvalue)")
grid.arrange(grobs=list(p1,p2,p3),ncol=3)
gene.grobs <- list()
for(gene in c(donorWise.T2D.RePACT$FishersMethod.q.dn.df$gene[1:3], donorWise.T2D.RePACT$FishersMethod.q.up.df$gene[1:3])){
    gene.grobs[[gene]] <- InSample.compare.nb(OBJ=scRNA.OBJ, pseudo=donorWise.T2D.RePACT$PCAInfo[,"pseudo.index.balanced",drop=F], gene=gene, GetDatatoplot=T)[[3]] %>% ggplot(.)+geom_violin(aes(Q,pseudo.index.balanced,fill=Ave.expr),alpha=0.9) + facet_grid(cols = vars(Sample),space="free")+scale_fill_gradient(low="white",high="red")+theme_classic()+ggtitle(gene)
}
grid.arrange(grobs=gene.grobs[1:3],ncol=3)
grid.arrange(grobs=gene.grobs[4:6],ncol=3)
dev.off()
```
![](https://github.com/chenweng1991/RePACT/blob/master/image/scRNA.RePACT.donorWise.PNG)
## RePACT on snATAC-seq
For example, we ran RePACT on beta cell snATAC-seq Seurat obj(donwload here).
```
snATAC.OBJ <- readRDS("../Beta.snATAC.rds")
T2D.snATAC.RePACT <- repact_logistic.snATAC(OBJ=scRNA.OBJ, LSIrange=1:20, pheno="diseaseStat", outputname)
pdf("T2D.snATAC.RePACT.pdf")
scatter3D(T2D.snATAC.RePACT$LSI.withinfo[,LSI_Top[1]], T2D.snATAC.RePACT$LSI.withinfo[,LSI_Top[2]], T2D.snATAC.RePACT$LSI.withinfo[,LSI_Top[3]],ticktype = "detailed", pch = 20, theta = 90, phi = 30, colvar = ifelse(LSI.withinfo[,pheno]==pheno.2,1,0), bty = "b2", cex = 0.3, col = alpha.col(col = c("steelblue", "red"), 0.6))
ggplot(T2D.snATAC.RePACT$LSI.withinfo)+aes(pseudo.index.balanced,fill=get(pheno))+geom_density()+scale_fill_manual(values=c("steelblue","red"))+theme_classic()+theme(legend.position="none")
ggplot(T2D.snATAC.RePACT$LSI.withinfo)+aes(Sample,pseudo.index.balanced,fill=get(pheno))+geom_violin()+geom_boxplot(width=0.2,outlier.shape = NA,notch=F,coef = 0,fill="grey25",color="grey75")+coord_flip()+theme_classic()+scale_fill_manual(values=c("steelblue", "red"))+theme_bw()+theme(legend.position="none")
# p4 <- ggplot(T2D.snATAC.RePACT$Evenbin.donorContribute)+aes(Sample,value,fill=disease)+geom_bar(stat="identity",color="black")+facet_grid(~evenfragbin)+theme(axis.text=element_blank(),axis.ticks=element_blank())+scale_fill_manual(values=c("steelblue","red"))+theme_bw()
apply(T2D.snATAC.RePACT$tmpT2D.diffPeaks.20bin.LSI$UPDN.toplot[,c(rownames(T2D.snATAC.RePACT$tmpT2D.diffPeaks.20bin.LSI$UP), rownames(tmpT2D.diffPeaks.20bin.LSI$DN))],2,function(x){scale(x)}) %>% melt() %>% ggplot()+aes(Var1,Var2,fill=value)+geom_tile()+scale_fill_gradient2(low="steelblue",mid="white",high="red")+theme_classic()+theme(axis.text=element_blank())+ggtitle("LSI-20bins")+xlab("20 trajectory bins") + ylab("ATAC peaks")
dev.off()
```

![](https://github.com/chenweng1991/RePACT/blob/master/image/T2D.snATAC.RePACT.PNG)

## RePACT on snATAC-seq to identify intra-dnor and inter-donor ATAC-peaks
```
snATAC.OBJ <- readRDS("../Beta.snATAC.rds")
donorWise.T2D.snATAC.RePACT <- DonorWise.snATAC.RePACT(OBJ=snATAC.OBJ, pheno="diseaseStat")
p1 <- ggplot(donorWise.T2D.snATAC.RePACT$FishersMethod.q.dn.df)+aes(rank,-log10(qvalueInSample),label=tag2,color=tag1)+geom_point()+geom_text_repel()+geom_hline(yintercept=2)+theme_classic()+ggtitle("DN peaks")
p2 <- ggplot(donorWise.T2D.snATAC.RePACT$FishersMethod.q.up.df)+aes(rank,-log10(qvalueInSample),label=tag2,color=tag1)+geom_point()+geom_text_repel()+geom_hline(yintercept=2)+theme_classic()+ggtitle("UP peaks")
p1 + p2

```
