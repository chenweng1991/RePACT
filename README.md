# EZsinglecell
---


## Install
```
# Install development version from GitHub:
# install.packages("devtools")

install("/mnt/NFS/homeGene/JinLab/cxw486/lib/dropseqlib/analysisscript/RePACT")


library("RePACT")


setwd("/mnt/NFS/homeGene/JinLab/cxw486/lib/dropseqlib/analysisscript/RePACT")
setwd("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/DGEanalysis/Islet412+511+919+T2D1+S4/workplaceLink/2017.9.21.revise")

```
---

## Introduction
We developed RePACT (Regressing Principle components for the Assembly of Continuous Trajectory) as a general method to sensitively identify disease relevant gene signatures using single cell data. The key step is to find the best trajectory to rank single cells (e.g. β cells) reflecting the change of disease status. In this study, we used RePACT to study obesity (denoted by a continuous BMI variable) and T2D (denoted by a dichotomous variable T2D).

Before running a main RePACT analysis, We highly recommend to perform cell clustering to generate a clear picture of cell type composition in the data. With cell type informrmation annotated, we can then focus on ONE specific cell type across different donors for a "clean" disease trajectory analysis.

We integrate some major functions from [Seurat](https://satijalab.org/seurat/) into our pipeline for basic dimension reduction and clustering analysis.

Here, we share our raw data as well as some necessary intermediate data to demonstrate the usage of RePACT toolkit.

Please follow the pipeline below that examplify the RePACT analysis to generate major results in our manuscript ([Zhou & Chen et al,2019](https://doi.org/10.1016/j.celrep.2019.02.043))

---

## Preliminary clustering analysis

We provide the digital gene expression(dge) matrix data from each donor. There are 9 dge matrix data in total.
**H1, H2, H3, H4, H5, H6** are healthy donors
**T2D1, T2D2 and T2D3** are T2D donors


##### 1. We share all dge raw data here

```
data(H1.D.clean.dge,H2.D.clean.dge,H3.D.clean.dge,H4.D.clean.dge,H5.D.clean.dge,H6.D.clean.dge,T2D1.D.clean.dge,T2D1.D.clean.dge,T2D2.D.clean.dge,T2D3.D.clean.dge)

#Data look like this
H1.D.clean.dge[1:5,1:5]
#         TGTGAGCTGAGA TTGATCTGCCCA GCTAACCTCTCN GTCTAATCCCGT TCTCACCCTTCN
#A1BG                3            2            1            0            2
#A1BG-AS1            0            0            0            0            1
#A1CF                0            0            1            0            0
#A2M                 0            0            0            0            0
#A2M-AS1             0            0            0            0            0
```


##### 2. Normalization and clustering for multiple samples
One line commend **_`docluster.multi()`_** for basic dimension reduction and clustering analysis.It takes multiple dge matrix data as input, and output a Seurat style metadata object. In this example, I took H1 and H2 as smallest data as a short example
```
H1H2.ob<-docluster.multi(Number=500,txcutoff=500,sets=list(H1.D.clean.dge,H2.D.clean.dge),nms=c("H1","H2"))
```
##### 3. Visulization
With the analyzed object above, we provide a commend **_`Fullplot_v2()`_** to generate major informative figures. **_`Fullplot_v2()`_** will generate a PDF file containing
- A. tSNE plot colored by default cluster
- B. heatmap showing signature genes for each cluster
- C. tSNE plot colored by samples (as 2 donors in this example)
- D. PCA plot colored by default cluster
- E. PCA plot colored by donors (as 2 donors in this example)
- F. Heatmap to show Genes that are driving PC1-PC4, respectively
- G. tSNE staining to show marker genes

```
H1H2.fullplot<-Fullplot_v2(H1H2.ob,"./PDF/example.fullplot.pdf",signiture=c("INS", "GCG", "SST", "PPY", "KRT19", "COL1A2"),doreturn=T)
```
The above two line will generate example PDFs for [H1H2](https://github.com/chenweng1991/RePACT/blob/RePACT.organized/PDF/example.fullplot.pdf)
> Example figure , left panel is the tsne illustration of cell clustering. It is colored by unsupervised cell cluster. Right panel is the same tSNE plot with interested genes highlighted . Color intensity reflect expression level in Z score.
![](https://raw.githubusercontent.com/chenweng1991/RePACT/RePACT.organized/image/Fig1.png)

##### 4. Cell type specific secondary clustering
For a robust RePACT, we suggest to redo clustering with a focus on one specific cell type across different Samples.

To extract sample information coupled with each single cell, which is the important info in RePACT to bridge between single cell and potential interesting phenotype
```
Sample.dict<-H1H2.ob@data.info[,"Sample",drop=F]
```
We show an example of secondary clustering on insulin-positive cells, as is in the case above we choose cluster 1 Note, Sample.dict extracted should be input as adjoint info to maintain the sample info for each single cell
```
H1H2.C1.raw<-as.matrix(H1H2.ob@raw.data[,row.names(subset(H1H2.ob@data.info,res.0.6==1))])
H1H2.C1.ob<-docluster.single(500,H1H2.C1.raw,dict=Sample.dict)
```
As above, you may want to visulize it by Fullplot_v2
```
H1H2.C1.fullplot<-Fullplot_v2(H1H2.C1.ob,"./PDF/H1H2.C1.fullplot.pdf",signiture=c("INS", "GCG", "SST", "PPY", "KRT19", "COL1A2"),doreturn=T,cell.use=100)
```

## RePACT main section
---
#### 1. Things to prepare before started
- First of all, the major input of RePACT pipeline is an preliminarily analyzed object from **_`docluster.multi()`_** as was shown in the previous section.  We suggest to use cell type specific data.
- A table with phenotype information that would be used to build the model.
In our case, the phenotype table look like below. When making your customized phenotype table, please follow this format, especially the first column, which is the sample name and must be the same as that was used in **_`docluster.multi()`_**, denoted in parameter `nms`. In our study, we used *BMI* as well as *Disease* to build the model respectively.
```
data(phenotable)
> phenotable
  Sample Gender Age   BMI HbA1c Disease
1     H1      M  27 20.60 5.40%       H
2     H2      M  21 22.80 5.20%       H
3     H3      F  38 34.40 5.00%       H
4     H4      M  52 22.00 5.60%       H
5     H5      M  28 30.80 4.90%       H
6     H6      M  44 34.60 5.40%       H
7   T2D1      M  58 39.30 8.90%       T
8   T2D2      M  61 28.10 5.20%       T
9   T2D3      M  51 35.59 7.10%       T
```
Below, we are using our full data for the RePACT demonstration because a good number of cells and samples with variable phenotypes are important to give enough power to uild a decent model. We share our full non-stressed beta cell dataset after filtering,
```
data(Allbeta.dge,sample.dict.full)
Beta.HSnegonly.ob<-docluster.single(500,Allbeta.dge,dict=Sample.dict.full)
```
#### 2. Build RePACT model




1.  Do clustering for one dge sample
```
docluster(dgepreprocess(s7.RockII_1.dge,500,norowname=T),
GetinformativeGene(dgepreprocess(s7.RockII_1.dge,500,norowname=T),500),
"s7.RockII_1",reso=0.6)
```

2.  Do clustering for combined multiple Samples, Number is informative gene number.
```
docluster.multi(Number=500,txcutoff=500,sets=list(s7.RockII_1=s7.RockII_1.dge,s7.B=s7.B.dge),nms=c("s7.RockII_1","s7.B"),israw=T)
```

### DE analysis
1. Make data pair that can be used for DE analysis and bubble plot
```
ROCKvsnorock.endo.paired<-datapair.mk(list(S7rock=S7rock_1.ob,S7=S7.ob),cols=c("Sample","Sample.2nd"),pick.list=list(c("s7.RockII_1"),c("s7.B_1")),normalizecellsize=F,randomizecelloirder=T)
```
2. Run negtive binomial based differentially expression gene calling
```
library(scran)  # Use computeSumFactors to compute size factor
ROCKvsnorock.endo.tri.dummy<-DE.gettripple(ROCKvsnorock.endo.paired,cpcol="name")
ROCKvsnorock.endo.de<-DoDE(ROCKvsnorock.endo.tri.dummy,"name",onlyoneSample=T,cpus=16)
```

### Ploting

1. Fullplot_v2,  plot basic clustering results for seurat object
```
mytopgenes<-Gettopgenes(object, number)  # Sometimes, precalculated cluster gene signatures can save some time
Fullplot_v2(S7rock_1.ob,"S7rock_1.ob.pdf",topgene=NULL,resolusion="res.0.6",signiture=c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A","DNAJB1","GHRL"))
```
2.  GettsnesignatureSuper, plot tsne or PCA stained by gene relative expression
```
GettsnesignatureSuper(object, object.all, signiture = c("INS", "RBP4","FFAR4", "ID1", "ID2", "ID3", "DNAJB1"), doPCA = F, dotSNE = T,usePC34 = F, extratitle = "", toreturn = F, dotsize = 0.3,buttomgrey = T, nolegend = T, highcolor = "red")
```

3. Generalbubbleplot,  to plot bubble plot for data pair, df is a dataframe whose row name is gene name, there is an extra columns named "tag", that is used for facet if usefacet is set as T.
```
Generalbubbleplot(ROCKvsnorock.non.paired,cpcol="name",genelist=df,donormalscale=F,usefacet=F)
```

### RePACT main functions

0. Loading packages
```
library(pscl)  # For Repact analysis to calculate logistic regression pseudo-p value
library(plot3D) # To plot 3D
```
1. Prepareforpseudoregress.g  This function is to perform the initial regression to prepare the trajectory study
```
T2D.tjct.ob<-Prepareforpseudoregress.g(Beta.HSnegonly.ob,PCrange=1:10,phenodic.use=phenodic,pheno="disease",linear=F)
```

2. Toplot3Dtjct   his function is to make 3D example plot for regression analysis
```
Toplot3Dtjct(T2D.tjct.ob,PCrange=c(1,3,4),pheno="disease",linear=F,multiplotname="test.pdf",titlename="tittle")
```


3. Tjct.core.gen  This function is to compute significant trajectory genes by linear regression
```
T2D.tjct.2nd.ob<-Tjct.core.gen(T2D.tjct.ob)
```

4. Tjct.core.plot  This function is to generate major plots for Repact analysis
```
Tjct.core.plot(T2D.tjct.ob,T2D.tjct.2nd.ob,pheno="Disease",f1.name="T2D.tjct.10d.violin.pdf",f2.name="T2D.tjct.his.pdf",f3.name="T2D.tjct.trj.heatmap.pdf",f3.height=14,f3.tittle="cell type:Changing genes on phenotype trajectory\ntop6%",table1.name="T2D.tjct.traj.up.genes-q0.05Full.csv",table2.name="T2D.tjct.traj.dowb.genes-q0.05Full.csv",rankcut=0.04)
```
### Data included
- Insulin regulator gene by Crispr-screening
  - Crisp.t1
  - Crisp.t2
- DB/OB related GWAS data
  - GWASdata
- Surfaceome data
  - Surfaceome.data
- Prepare cell cycle data
  - G1.S
  - S
  - G2.M
  - M
  - M.G1
  - all.cellcycle
- Transcription factor
  - TFfromDBD(Old version)
  - TFvector(version 18.12.18)
- GSEA msigdb
  - msig.db
- Beta/alpha/delta/pp/episilon specific geneset
  - all.Beta.maker
  - allAlpha.marker
  - delta.markers
  - pp.markers
  - epsilon.markers
- Diabetes and obesity trajectory genes
  - Beta.BMI.trj.genes
  - Beta.T2D.trj.genes
- Druggable gene
  - gene_drug.db
- Hocomoco genes
  - Hocomocogenes




  ### To do

  incorporate the  library
  library(Seurat, lib.loc="/mnt/NFS/homeGene/JinLab/cxw486/R/x86_64-redhat-linux-gnu-library/3.3")
