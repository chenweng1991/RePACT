# Most frequently used analysis and plotting for single cells RNAseq

### Loading packages
```
library(Seurat)  #‘1.4.0.16’
library(grid)
library(gridExtra)
library(gtable)
library(dplyr)
library(RColorBrewer)
library(EZsinglecell)
```

### Loading example dge data
```
example.dge<-read.table("./EZsinglecell/test.data/Example.dge.txt.gz",header=TRUE,stringsAsFactors=FALSE)
example.dge[,1]<-gsub("-",".",example.dge[,1])
Example.ob<-readRDS("./EZsinglecell/test.data//N706_Stage6.obj")
```

### Basic clustering plots
```
Fullplot_v2(Example.ob,"test.pdf",p.tsnesample=F,p.pcasample=F)
```
