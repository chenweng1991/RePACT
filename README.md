# EZsinglecell



### Install
```
# Install development version from GitHub:
# install.packages("devtools")

install("/mnt/NFS/homeGene/JinLab/cxw486/lib/dropseqlib/analysisscript/RePACT")


library("RePACT")


setwd("/mnt/NFS/homeGene/JinLab/cxw486/lib/dropseqlib/analysisscript/RePACT")
setwd("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/DGEanalysis/Islet412+511+919+T2D1+S4/workplaceLink/2017.9.21.revise")

```
---

### Introduction
We developed RePACT (Regressing Principle components for the Assembly of Continuous Trajectory) as a general method to sensitively identify disease relevant gene signatures using single cell data. The key step is to find the best trajectory to rank single cells (e.g. β cells) reflecting the change of disease status. In this study, we used RePACT to study obesity (denoted by a continuous BMI variable) and T2D (denoted by a dichotomous variable T2D).

Before running a main RePACT analysis, We highly recommend to perform cell clustering to generate a clear picture of cell type composition in the data. With cell type informrmation annotated, we can then focus on ONE specific cell type across different donors for a "clean" disease trajectory analysis.

We integrate some major functions from [Seurat](https://satijalab.org/seurat/) into our pipeline for basic dimension reduction and clustering analysis.

Here, we share our raw data as well as some necessary intermediate data to demonstrate the usage of RePACT toolkit.

Please follow the pipeline below that examplify the RePACT analysis to generate major results in our manuscript ([Zhou & Chen et al,2019](https://doi.org/10.1016/j.celrep.2019.02.043))


### Clustering analysis



We provide the digital gene expression(dge) matrix data from each donor. There are 9 dge matrix data in total.
H1, H2, H3, H4, H5, H6 are healthy donors
T2D1, T2D2 and T2D3 are T2D donors

To gain all dge raw data
```
data(H1.D.clean.dge,H2.D.clean.dge,H3.D.clean.dge,H4.D.clean.dge,H5.D.clean.dge,H6.D.clean.dge,T2D1.D.clean.dge,T2D1.D.clean.dge,T2D2.D.clean.dge,T2D3.D.clean.dge)

#Data look like this
H1.D.clean.dge[1:5,1:5]
         TGTGAGCTGAGA TTGATCTGCCCA GCTAACCTCTCN GTCTAATCCCGT TCTCACCCTTCN
A1BG                3            2            1            0            2
A1BG-AS1            0            0            0            0            1
A1CF                0            0            1            0            0
A2M                 0            0            0            0            0
A2M-AS1             0            0            0            0            0
```



One line commend **_docluster.multi()_** for basic dimension reduction and clustering analysis.It takes multiple dge matrix data as input, and output a Seurat style metadata object. In this example, I took H1 and H2 as smallest data as a fast example
```
H1H2.ob<-docluster.multi(Number=500,txcutoff=500,sets=list(H1.D.clean.dge,H2.D.clean.dge),nms=c("H1","H2"))
```
Alternatively, the full analyzed object that include 9 donors can be directly loaded
```
data(SeuratALL.filtered.0.6)
```
With the analyzed object above, we provide a commend **_Fullplot_v2_** to generate major informative figures. **_Fullplot_v2_** will generate a PDF file containing
A. tSNE plot colored by default cluster
B. heatmap showing signature genes for each cluster
C. tSNE plot colored by samples (as 2 donors in this example)
D. PCA plot colored by default cluster
E. PCA plot colored by donors (as 2 donors in this example)
F. Heatmap to show Genes that are driving PC1-PC4, respectively
G. tsne staining to show marker genes
```
H1H2.fullplot<-Fullplot_v2(H1H2.ob,"./PDF/example.fullplot.pdf",signiture=c("INS", "GCG", "SST", "PPY", "KRT19", "COL1A2"),doreturn=T)
```
The above two line will generate example PDFs for [H1H2]()

<!-- tiff("./image/Fig1A.tiff",res=300, width = 6, height = 6, units = 'in')
H1H2.fullplot[[1]]
dev.off()

newloist<-list()
for (i in 1:6){
newloist<-c(newloist,list(H1H2.fullplot[[9]][[i]]+theme(axis.text=element_blank(),axis.title=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())))
}




tiff("./image/Fig1G.tiff",res=300, width = 6, height = 6, units = 'in')
grid.arrange(grobs=newloist,nrow=2)
dev.off() -->

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
