# RePACT
---


## Install
```
# Install development version from GitHub:
# install.packages("devtools")

install("/mnt/NFS/homeGene/JinLab/cxw486/lib/dropseqlib/analysisscript/RePACT")


library("RePACT")

devtools::install_github("chenweng1991/RePACT")

library(Seurat, lib.loc="/mnt/NFS/homeGene/JinLab/cxw486/R/x86_64-redhat-linux-gnu-library/3.3")
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


#### 1. We share all dge raw data here

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


#### 2. Normalization and clustering for multiple samples
One line commend **_`docluster.multi()`_** for basic dimension reduction and clustering analysis.It takes multiple dge matrix data as input, and output a Seurat style metadata object. In this example, I took H1 and H2 as smallest data as a short example
```
H1H2.ob<-docluster.multi(Number=500,txcutoff=500,sets=list(H1.D.clean.dge,H2.D.clean.dge),nms=c("H1","H2"))
```
#### 3. Visulization
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

#### 4. Cell type specific secondary clustering
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
Below, we are using our beta cell only dataset for the RePACT demonstration because a good number of cells and samples with variable phenotypes are essential to give enough power to build a decent model. We share our full non-stressed beta cell dataset after filtering here.
```
data(Allbeta.dge,sample.dict.full)
Beta.HSnegonly.ob<-docluster.single(500,Allbeta.dge,dict=sample.dict.full)
Fullplot_v2(Beta.HSnegonly.ob,"./PDF/Beta.HSnegonly.fullplot.pdf",signiture=NULL,doreturn=T,cell.use=100)
```
The analyzed only beta cells look like [this](https://github.com/chenweng1991/RePACT/blob/RePACT.organized/PDF/Beta.HSnegonly.fullplot.pdf)

#### 2. Build RePACT model
The goal is to create a model that fit the relationship between the phenotype and single cell transcriptome. We provide
two regression models. For binary phenotype, we fit the model by logistic regression, for example disease status, samples are either healthy or diabetic, denoted as H and T (`linear=F`). For continuous phenotype for example,BMI, samples have different index values. We fit the model using linear regression(`linear=T`). The parameter `PCrange` determine how many PCs to use for the phenotype prediction. The default PC numbers are PC1 through PC10
```
T2D.tjct.ob<-Prepareforpseudoregress.g(Beta.HSnegonly.ob,PCrange=1:10,phenodic.use=phenotable,pheno="Disease",linear=F)
BMI.tjct.ob<-Prepareforpseudoregress.g(Beta.HSnegonly.ob,PCrange=1:10,phenodic.use=phenotable,pheno="BMI",linear=T)
```
I will use `T2D.tjct.ob` as an example to explain the RePACT model object,.
- Use `BMI.tjct.ob$model` to check out the model and significance of each PC as a phenotype predictor
```
BMI.tjct.ob$model.para
# Call:
# lm(formula = form)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -14.4640  -2.2229   0.3886   2.6119  15.3527
#
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
# (Intercept)     30.12522    0.06325 476.280  < 2e-16 ***
# PCandPheno$PC1   0.49493    0.02053  24.108  < 2e-16 ***
# PCandPheno$PC2   1.42528    0.02686  53.064  < 2e-16 ***
# PCandPheno$PC3   0.39394    0.03201  12.309  < 2e-16 ***
# PCandPheno$PC4   0.49144    0.03320  14.803  < 2e-16 ***
# PCandPheno$PC5  -0.38479    0.03393 -11.341  < 2e-16 ***
# PCandPheno$PC6   0.18351    0.03783   4.850 1.28e-06 ***
# PCandPheno$PC7  -0.27799    0.03870  -7.183 8.13e-13 ***
# PCandPheno$PC8   0.46673    0.04217  11.068  < 2e-16 ***
# PCandPheno$PC9   0.12591    0.04318   2.916  0.00357 **
# PCandPheno$PC10  0.75956    0.04427  17.157  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 3.95 on 3889 degrees of freedom
# Multiple R-squared:  0.5306,	Adjusted R-squared:  0.5294
# F-statistic: 439.7 on 10 and 3889 DF,  p-value: < 2.2e-16
```
- For a partial illustration `BMI.tjct.ob$reg.plot.2d` shows the corelation only between PC1,PC2 and the phenotype (BMI). But keep in mind the actual model is on 10 dimensions
![](https://raw.githubusercontent.com/chenweng1991/RePACT/RePACT.organized/image/BMI.2dplot.png)

*Optional* For a 3D illustration of phenotype against any three PCs, use commend `Toplot3Dtjct`
```
Toplot3Dtjct(T2D.tjct.ob,PCrange=c(1,3,4),pheno="Disease",linear=F,multiplotname="PDF/test.pdf",titlename="tittle")
```

#### 3. Reconstruct disease trajectory
If RePACT model tells that single cell transcriptome variation,on PC space, is highly correlated with the phenotype you are interested, the next step is to calculate the phenotype associated pseudostates index, by which I will be able to identify genes that are differentially expressed across the disease trajectory.

```
T2D.tjct.2nd.ob<-Tjct.core.gen(T2D.tjct.ob,binnumber=20)
```
I will explain the most important components in this trajectory object using `T2D.tjct.2nd.ob` as an example
- T2D.tjct.2nd.ob$raw.bin[[1]] is a raw digital expression matrix with a column `pseudo.index`, which is the calculated index value for each single cell on the disease associated trajectory. `tag` is the bin number after put cells into several bins on that trajectory(along pseudo.index). T2D.tjct.2nd.ob$raw.bin[[2]] shows how the bins are defined
- T2D.tjct.2nd.ob$BINlinear.result.summarized$UP is a dataframe showing the genes identified as upregulated along the disease associated trjectory. T2D.tjct.2nd.ob$BINlinear.result.summarized$DOWN are the genes identified as down-regulated.

#### 4. Generate and visulize intepretable RePACT results

```
Tjct.core.plot(T2D.tjct.ob,T2D.tjct.2nd.ob,pheno="Disease",f1.name="T2D.tjct.10d.violin.pdf",f2.name="T2D.tjct.his.pdf",f3.name="T2D.tjct.trj.heatmap.pdf",f3.height=10,f3.tittle="Beta cell:Variable genes on T2D trajectory\ntop1%\nLow T2D ---------------> High T2D",table1.name="T2D.tjct.traj.up.genes-q0.05Full.csv",table2.name="T2D.tjct.traj.dowb.genes-q0.05Full.csv",rankcut=0.01,colorset="Set1")
```
This command will generate 3 figures and two tables
- [T2D.tjct.10d.violin.pdf](https://github.com/chenweng1991/RePACT/blob/RePACT.organized/PDF/T2D.tjct.10d.violin.pdf)  The violin plot like was shown in our manuscript, shows the single cell distribution along the pseudo-index
- [T2D.tjct.his.pdf](https://github.com/chenweng1991/RePACT/blob/RePACT.organized/PDF/T2D.tjct.his.pdf)  This is also showing distribution along pseudo-index, but by histogram
- [T2D.tjct.trj.heatmap.pdf](https://github.com/chenweng1991/RePACT/blob/RePACT.organized/PDF/T2D.tjct.trj.heatmap.pdf) The heatmap to show top X% genes as `rankcut` decided that is highly variable along the trajectory.
- [T2D.tjct.traj.dowb.genes-q0.05Full.csv](https://github.com/chenweng1991/RePACT/blob/RePACT.organized/PDF/T2D.tjct.traj.dowb.genes-q0.05Full.csv)  Gene list of significantly down-regulated along trajectory
-  [T2D.tjct.traj.up.genes-q0.05Full.csv](https://github.com/chenweng1991/RePACT/blob/RePACT.organized/PDF/T2D.tjct.traj.up.genes-q0.05Full.csv)  Gene list of significantly up-regulated along trajectory
>T2D.tjct.10d.violin and T2D.tjct.trj.heatmap.pdf
![](https://raw.githubusercontent.com/chenweng1991/RePACT/RePACT.organized/image/RePACT.violinheat.png))






  ### To do

  incorporate the  library
  library(Seurat, lib.loc="/mnt/NFS/homeGene/JinLab/cxw486/R/x86_64-redhat-linux-gnu-library/3.3")
