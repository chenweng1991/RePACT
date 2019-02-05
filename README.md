# EZsinglecell



### Install
```
# Install development version from GitHub:
# install.packages("devtools")
#devtools::install_github("chenweng1991/RePACT")
```
---

### Clustering

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
