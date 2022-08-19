# RePACT 
## (Regressing Principle components for the Assembly of Continuous Trajectory)

## Introduction

We developed RePACT as a general method to sensitively identify disease relevant gene signatures / ATAC-peaks using single cell data.The goal of RePACT is to find the best trajectory to rank single cells (e.g. in our case, the β cells) reflecting the change of disease status. In this study, we used RePACT to study  type II diabetes(T2D) (denoted by a dichotomous variable T2D).Please follow the pipeline below that serve as a tutorial for RePACT analysis. It generates major results in our manuscript ([Zhou & Chen et al,2019](https://doi.org/10.1016/j.celrep.2019.02.043)).

## Install
```
devtools::install_github("chenweng1991/RePACT") #install RePACT
library("RePACT")
```
## RePACT on scRNA-seq
It requires a Seurat object for a specific cell type you're interested, and the meta.data in the object should contain "Sample" column and another column to compare using RePACT. 
For example, we ran RePACT on beta cell scRNA-seq Seurat obj(donwload here), part of meta.data are shown below. Note, the column containg donor ID has to be named as "Sample".

```
  Sample Gender Age   BMI HbA1c diseaseStat
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
Examples to run RePACT on the downloaded Seurat object.
```
scRNA.OBJ <- readRDS("../Old9.Beta.scRNA.rds")
repact.scRNA(scRNA.OBJ, PCrange=1:20,pheno="diseaseStat",is_continuous="F",norm_index=T,SlopeCut=0.05,output_name2="T2D_Beta.scRNA.RePACT")
```
Output files will be written out to the disk, they are intermediate rds files, genes on T2D trajectory, plots.
```
T2D_Beta.scRNA.RePACT.RePACT.pdf
T2D_Beta.scRNA.RePACTDOWN_gene.csv
T2D_Beta.scRNA.RePACT.genes.FCs.csv
T2D_Beta.scRNA.RePACT.genes.Full_list.csv
T2D_Beta.scRNA.RePACT.RepactOBJ.1.rds
T2D_Beta.scRNA.RePACT.RepactOBJ.2.rds
T2D_Beta.scRNA.RePACTUP_gene.csv
```
![](https://raw.githubusercontent.com/chenweng1991/RePACT/RePACT.organized/image/RePACT.violinheat.png))

## RePACT on snATAC-seq
It requires a Seurat object for a specific cell type you're interested, the raw data should be "cell~peak" format, and the meta.data in the object should contain "Sample" column and another column to compare using RePACT. 

For example, we ran RePACT on beta cell snATAC-seq Seurat obj(donwload here).
```
snATAC.OBJ <- readRDS("../Beta.snATAC.rds")
repact_logistic.snATAC(snATAC.OBJ, LSIrange=1:20, pheno="diseaseStat", outputname="T2D_Beta.snATAC.RePACT")
```
Output files will be written out to the disk
```
T2D_Beta.snATAC.RePACT.CARePACT.pdf
T2D_Beta.snATAC.RePACT.CARePACT.T2D_diffPeaks.20bin.LSI.rds
T2D_Beta.snATAC.RePACT.T2D.diffPeaks.repact.bulk.rds
```
![](https://github.com/chenweng1991/RePACT/blob/81317e850892b415a554c97433646b8abfbde9f0/image/RePACT.snATAC.graphic.PNG)

