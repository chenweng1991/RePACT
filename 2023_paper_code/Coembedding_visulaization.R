  Islet12ALL.scATAC.FILTER<-readRDS("Islet12ALL.scATAC.FILTER")
  library(e1071)
  Islet.scATAC.metadata.HiQ<-subset(Islet.scATAC.metadata,PeakRatio>0.25 & Fragments>2000)
  Islet12ALL.scATAC.FILTER<-RunscATACclustering.Filter(Mergedset=Islet.scATAC.dac.biggenes.Mtx[,Islet.scATAC.metadata.HiQ$Name],atacmeta=Islet.scATAC.metadata.HiQ,ref=Islet12.scRNA.seurat3.filtered,cl="Cell_type")
  saveRDS(Islet12ALL.scATAC.FILTER,"/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/Islet12ALL.scATAC.FILTER.HiQ")
  # saveRDS(Islet12ALL.scATAC.FILTER,"/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/Islet12ALL.scATAC.FILTER")
  Islet12ALL.scATAC.FILTER<-readRDS("/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/Islet12ALL.scATAC.FILTER")
  ## Step 1 Split the rna and atac
  cca.rna<-Islet12ALL.scATAC.FILTER$cca.1000.L2normed[!grepl("TA:Z",row.names(Islet12ALL.scATAC.FILTER$cca.1000.L2normed)),]
  cca.atac<-Islet12ALL.scATAC.FILTER$cca.1000.L2normed[grepl("TA:Z",row.names(Islet12ALL.scATAC.FILTER$cca.1000.L2normed)),]
  cca.rna.celltype<-Tomerge_v2(cca.rna,Islet12.scRNA.seurat3.filtered@meta.data[,"Cell_type",drop=F])
  ##region(Optional) To test SVM, Multinomial logistic regression, and random forest model accuracy, optional
  accuracy.SVM<-c()
  for (i in 1:50){
    print(i)
  s<-sample(nrow(cca.rna.celltype),6000)
  cca.rna.celltype_train<-cca.rna.celltype[s,]
  cca.rna.celltype_test<-cca.rna.celltype[-s,]
  cca.rna.celltype_train$Cell_type<-as.factor(cca.rna.celltype_train$Cell_type)
  svmfit<-svm(Cell_type~.,data=cca.rna.celltype_train, kernel="radial",cost=10,scale=T)
  # tuned<-tune(svm,y ~.,data=cca.rna.celltype_train, kernel="linear",ranges=list(cost=c(0.001,0.01,0.1,1,10,100)))
  # summary(tuned)
  p<-predict(svmfit,cca.rna.celltype_test[,1:30],type="class")
  accuracy.SVM<-c(accuracy.SVM,length(which(cca.rna.celltype_test$Cell_type==p))/nrow(cca.rna.celltype_test))
  }
  #################Test random forest ########################################
  library(randomForest)
  accuracy.rf<-c()
  for (i in 1:50){
    print (i)
  s<-sample(nrow(cca.rna.celltype),6000)
  cca.rna.celltype_train<-cca.rna.celltype[s,]
  cca.rna.celltype_test<-cca.rna.celltype[-s,]
  cca.rna.celltype_train$Cell_type<-as.factor(cca.rna.celltype_train$Cell_type)
  rf <- randomForest(
    Cell_type ~ .,
    data=cca.rna.celltype_train
  )
  p = predict(rf, newdata=cca.rna.celltype_test[,1:30])
  accuracy.rf<-c(accuracy.rf,length(which(cca.rna.celltype_test$Cell_type==p))/nrow(cca.rna.celltype_test))
  }
  ############## Test Mutilnomial logistic regression########################
  library(nnet)
  accuracy.mlr<-c()
  for (i in 1:50){
    print (i)
  s<-sample(nrow(cca.rna.celltype),6000)
  cca.rna.celltype_train<-cca.rna.celltype[s,]
  cca.rna.celltype_test<-cca.rna.celltype[-s,]
  cca.rna.celltype_train$Cell_type<-as.factor(cca.rna.celltype_train$Cell_type)
  mlr <- multinom(as.formula(paste("Cell_type ~ ",paste(paste("CC",1:30,sep="_"),collapse="+"))), data = cca.rna.celltype_train)
  p = predict(mlr, newdata=cca.rna.celltype_test[,1:30])
  accuracy.mlr<-c(accuracy.mlr,length(which(cca.rna.celltype_test$Cell_type==p))/nrow(cca.rna.celltype_test))
  }
  data.frame(SVM=accuracy.SVM,LogisticRegression=accuracy.mlr,RandomForest=accuracy.rf) %>% melt %>% ggplot()+aes(variable,value)+geom_boxplot()+geom_jitter()+theme(axis.text=element_text(size=24,color="black"))+ggtitle("Prediction model: 6000(Out of 16,000 cells) as trainingg set, The rest cells for test")
  ##endregion Step 2  To test SVM, Multinomial logistic regression, and random forest model accuracy, optional
  ## Step2: Predict celltypes using the trained model SVM
  # s<-sample(nrow(cca.rna.celltype),6000)
  cca.rna.celltype_train<-cca.rna.celltype#[s,]
  svmfit<-svm(Cell_type~.,data=cca.rna.celltype_train, kernel="radial",cost=10,scale=T)
  p1<-predict(svmfit,cca.atac,type="class",decision.value=T) ##Actually
  # saveRDS(p1,"/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/p1.svm")
  p1<-readRDS("/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/p1.svm")
  ## Step 3 Visulize on UMAP
  myconfig<-umap.defaults
  myconfig$a<-5
  myconfig$b<-0.5
  atacrna.umap<-umap(Islet12ALL.scATAC.FILTER$cca.1000.L2normed,config=myconfig)
  Islet12.umap.meta<-data.frame(atacrna.umap$layout)
  names(Islet12.umap.meta)<-c("UMAP1","UMAP2")
  Islet12.umap.meta$tech<-ifelse(grepl("TA:Z",row.names(Islet12.umap.meta)),"ATAC","RNA")
  Islet12.umap.meta<-Tomerge_v2(Islet12.umap.meta,Islet12.scRNA.seurat3.filtered@meta.data[,"Cell_type",drop=F],leavex=T)
  Islet12.umap.meta<-Tomerge_v2(Islet12.umap.meta,data.frame(p1))
  plot1<-Islet12.umap.meta[order(Islet12.umap.meta$tech,na.last=F),]%>% ggplot(.)+aes(UMAP1,UMAP2,color=Cell_type)+geom_point(size=0.5)+ guides(colour = guide_legend(override.aes = list(size=10)))+ggtitle("RNA CellType Label")
  plot2<-subset(Islet12.umap.meta,tech=="ATAC") %>% ggplot(size=0.5)+aes(UMAP1,UMAP2,color=p1)+geom_point(size=0.5)+ guides(colour = guide_legend(override.aes = list(size=10)))+ggtitle("Predicted ATAC CellType Label")
  SVM.KNN<-Tomerge_v2(subset(Islet12.umap.meta,tech=="ATAC"),Islet12ALL.scATAC.FILTER$DoubletsAnalysis[,c("First","Second","ratio","Fragments")])
  # saveRDS(SVM.KNN,"/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/SVM.KNN")
  Doublets1<- subset(SVM.KNN,ratio<=2) %>% row.names
  Doublets2<-subset(SVM.KNN,First %in% c("Beta","Alpha","PSC") & Second %in% c("Beta","Alpha","PSC")) %>% subset(.,ratio<=4) %>% row.names
  Doublets3<-subset(SVM.KNN,First %in% c("PSC")) %>% subset(.,ratio<=100000) %>% row.names
  Doublets5<-subset(SVM.KNN,First %in% c("Duct")) %>% subset(.,ratio<=100000) %>% row.names
  Doublets4<-subset(SVM.KNN,First %in% c("Beta","Alpha","Delta","PP") & Second %in% c("PSC")) %>% subset(.,ratio<=16) %>% row.names
  Doublets5<-subset(SVM.KNN,First %in% c("Alpha") & Second %in% c("PSC")) %>% subset(.,ratio<=100000) %>% row.names
  Doublets6<-subset(SVM.KNN,First %in% c("PP") & Second %in% c("Alpha")) %>% subset(.,ratio<=4) %>% row.names
  Doublets7<-subset(SVM.KNN,First %in% c("Beta","Alpha","Delta") & Second %in% c("Duct")) %>% subset(.,ratio<=16) %>% row.names
  Doublets8<-subset(SVM.KNN,First %in% c("Beta","Alpha") & Second %in% c("Duct")) %>% subset(.,ratio<=100000) %>% row.names
  # Doublets3<-subset(SVM.KNN,First %in% c("Beta","Alpha") & Second=="PSC") %>% subset(.,ratio<=100) %>% row.names
  # Doublets4<-subset(SVM.KNN,First %in% c("Beta","Alpha") & Second=="PSC") %>% subset(.,ratio<=100) %>% row.names
  SVM.KNN.filtered<-SVM.KNN[!row.names(SVM.KNN) %in% unique(c(Doublets1,Doublets2,Doublets3,Doublets4,Doublets5,Doublets6,Doublets7,Doublets8)),]
  SVM.KNN.filtered<-SVM.KNN.filtered[as.character(SVM.KNN.filtered$p1)==as.character(SVM.KNN.filtered$First),]

PSC.TominusAlpha<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="PSC",]),] %>% .[order(abs(.[,"Alpha/PSC"])),] %>% .[1:400,"Alpha/PSC"] %>% names  # Saved
Alpha.TominusPSC<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="Alpha",]),] %>% .[order(abs(.[,"Alpha/PSC"])),] %>% .[1:500,"Alpha/PSC"] %>% names # saved
SVM.KNN.filtered<-SVM.KNN.filtered[!row.names(SVM.KNN.filtered) %in% c(Alpha.TominusPSC,PSC.TominusAlpha,PP.TominusAlpha),]
readRDS("")
SVM.KNN.filtered<-readRDS("/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/SVM.KNN.filtered")
PP.TominusAlpha<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="PP",]),] %>% .[order(abs(.[,"Alpha/PP"])),] %>% .[1:50,"Alpha/PP"] %>% names  # Not saved yet
Duct.TominusAlpha<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="Duct",]),] %>% .[order(abs(.[,"Alpha/Duct"])),] %>% .[1:800,"Alpha/Duct"] %>% names  #300
Duct.TominusBeta<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="Duct",]),] %>% .[order(abs(.[,"Beta/Duct"])),] %>% .[1:900,"Beta/Duct"] %>% names  #250
Beta.TominusPSC<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="Beta",]),] %>% .[order(abs(.[,"Beta/PSC"])),] %>% .[1:250,"Beta/PSC"] %>% names
Alpha.TominusPSC<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="Alpha",]),] %>% .[order(abs(.[,"Alpha/PSC"])),] %>% .[1:250,"Alpha/PSC"] %>% names
Alpha.TominusDuct<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="Alpha",]),] %>% .[order(abs(.[,"Alpha/Duct"])),] %>% .[1:250,"Alpha/Duct"] %>% names
Delta.TominusDuct<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="Delta",]),] %>% .[order(abs(.[,"Duct/Delta"])),] %>% .[1:40,"Duct/Delta"] %>% names
PSC.TominusAlpha<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="PSC",]),] %>% .[order(abs(.[,"Alpha/PSC"])),] %>% .[1:380,"Alpha/PSC"] %>% names
# PSC.TominusBeta<-attr(p1,"decision.values")[row.names(SVM.KNN.filtered[SVM.KNN.filtered$p1=="PSC",]),] %>% .[order(abs(.[,"Beta/PSC"])),] %>% .[1:50,"Beta/PSC"] %>% names
SVM.KNN.filtered.final<-SVM.KNN.filtered[!row.names(SVM.KNN.filtered) %in% c(PP.TominusAlpha,Duct.TominusAlpha,Duct.TominusBeta,Beta.TominusPSC,Alpha.TominusPSC,Alpha.TominusDuct,Delta.TominusDuct,PSC.TominusAlpha),]   #,PSC.TominusBeta
  # Doublets4<-subset(SVM.KNN.filtered,p1 %in% c("Alpha","Beta","PSC","Duct") &Fragments<2000) %>% row.names
  # SVM.KNN.filtered<-SVM.KNN.filtered[!row.names(SVM.KNN.filtered) %in% Doublets4,]

## Visulaize
  ##region 2020-8-28-b Visulize the clustering umap
  Islet.scATAC.dac.biggenes.Mtx<-readRDS("Islet.scATAC.dac.biggenes.Mtx")
  Islet12.scRNA.seurat3.filtered<-readRDS("Islet12.scRNA.seurat3.filtered")
  Islet.scATAC.metadata<-read.table("Islet.scATAC.metadata")
  row.names(Islet.scATAC.metadata)<-Islet.scATAC.metadata$Name
  SVM.KNN.filtered<-readRDS("SVM.KNN.filtered")
  CellToremove.1<-row.names(subset(SVM.KNN.filtered,p1 %in% c("Alpha","Beta","PSC","Duct") & Fragments<3000))
  CellToremove.2<-row.names(rbind(subset(SVM.KNN.filtered,p1=="Alpha" &Second=="PP"),subset(SVM.KNN.filtered,p1=="Alpha" &Second=="Delta"),subset(SVM.KNN.filtered,p1=="Delta" &Second=="PP"),subset(SVM.KNN.filtered,p1=="PP" &Second=="Delta")))
  CellToremove.3<-subset(SVM.KNN.filtered,p1=="Alpha" &ratio<10) %>% row.names
  SVM.KNN.filtered.Visulize<-SVM.KNN.filtered[!row.names(SVM.KNN.filtered) %in% c(CellToremove.1,CellToremove.2,CellToremove.3),]

Islet12.scRNA.seurat3.filtered.sub<-readRDS("Islet12.scRNA.seurat3.filtered.sub")
Islet12ALL.coembed.Visulize<-RunscATACclustering.v2(Mergedset=Islet.scATAC.dac.biggenes.Mtx[,row.names(SVM.KNN.filtered.Visulize)],atacmeta=subset(Islet.scATAC.metadata,Name %in% row.names(SVM.KNN.filtered.Visulize)),ref=Islet12.scRNA.seurat3.filtered.sub,kw=25)
Islet12ALL.coembed.Visulize[[1]] <- FindNeighbors(Islet12ALL.coembed.Visulize[[1]], dims = 1:30)
Islet12ALL.coembed.Visulize[[1]] <- FindClusters(Islet12ALL.coembed.Visulize[[1]], resolution = 0.2)

Islet12ALL.coembed.Visulize.meta<-Tomerge_v2(Islet12ALL.coembed.Visulize[[1]]@meta.data,SVM.KNN.filtered.Visulize[,"p1",drop=F]) %>% Tomerge_v2(.,Islet12ALL.coembed.Visulize[[1]]@reductions$umap@cell.embeddings)
Islet12ALL.coembed.Visulize.meta.atac<-subset(Islet12ALL.coembed.Visulize.meta,tech=="ATAC")
Islet12ALL.coembed.Visulize.meta.atac$Clustering.celltype<-mapvalues(Islet12ALL.coembed.Visulize.meta.atac$integrated_snn_res.0.2,from=c(0,8,11,1,6,4,5,2,3,7,9),to=c("Beta","Beta","Beta","Alpha","Alpha","PP|Delta","Delta","Duct","PSC","Acinar","Endothelial"))
Islet12ALL.coembed.Visulize.meta.atac.filtered<-Islet12ALL.coembed.Visulize.meta.atac[apply(Islet12ALL.coembed.Visulize.meta.atac,1,function(x){grepl(x[16],x[13])}),]
Final.filter<-c(row.names(subset(Islet12ALL.coembed.Visulize.meta.atac.filtered,p1=="Delta" & UMAP_1 > 2.5)),row.names(subset(Islet12ALL.coembed.Visulize.meta.atac.filtered,p1=="PP" & UMAP_1 < 0)))
Islet12ALL.coembed.Visulize.meta.atac.filtered<-Islet12ALL.coembed.Visulize.meta.atac.filtered[!row.names(Islet12ALL.coembed.Visulize.meta.atac.filtered) %in% Final.filter,]
p1data<-Tomerge_v2(Islet12ALL.coembed.Visulize[[1]]@reductions$umap@cell.embeddings,Islet12ALL.coembed.Visulize[[1]]@meta.data)
CellToremove<-row.names(subset(p1data,tech=="ATAC"))[!row.names(subset(p1data,tech=="ATAC")) %in% row.names(Islet12ALL.coembed.Visulize.meta.atac.filtered)]
p1data<-p1data[!row.names(p1data) %in% CellToremove,]
p1<-ggplot(data=subset(p1data,tech!="ATAC"),aes(UMAP_1,UMAP_2,color=celltype))+geom_point(size=1)+theme_classic()+ guides(color = guide_legend(override.aes = list(size=5)))+scale_color_brewer(palette="Set1")+geom_point(data=subset(p1data,tech=="ATAC"),aes(UMAP_1,UMAP_2),color="azure3",size=1.2,alpha=0.25)+theme(axis.line=element_line(size=0.5),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+ guides(color = guide_legend(override.aes = list(size=5)))
p1data<-Tomerge_v2(p1data,Islet12ALL.coembed.Visulize.meta.atac.filtered[,"p1",drop=F],leavex=T)
p1.highlightATAC<-ggplot(data=subset(p1data,tech!="RNA"),aes(UMAP_1,UMAP_2,color=p1))+geom_point(size=1)+theme_classic()+ guides(color = guide_legend(override.aes = list(size=5)))+scale_color_brewer(palette="Set1")+geom_point(data=subset(p1data,tech=="RNA"),aes(UMAP_1,UMAP_2),color="azure3",size=0.6,alpha=0.2)+theme(axis.line=element_line(size=0.5),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+ guides(color = guide_legend(override.aes = list(size=5)))
p2<-ggplot(Islet12ALL.coembed.Visulize.meta.atac.filtered)+aes(UMAP_1,UMAP_2,color=p1)+geom_point(size=1.2)+theme_classic()+theme(axis.line=element_line(size=0.5),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+ guides(color = guide_legend(override.aes = list(size=5)))+scale_color_brewer(palette="Set1")
setwd("~/Chip-seq/ATAC/scATAC/IsletscATAC/Fig1.Clustering")
png("Coembeding.ALL.png",height=1650,width=2000,res=300)
print(p1)
dev.off()

png("Coembeding.ALL.HLatac.png",height=1650,width=2000,res=300)
print(p1.highlightATAC)
dev.off()

png("Coembeding.ATAC.png",height=1650,width=2000,res=300)
print(p2)
dev.off()
##endregion

  ##region 2020-9-13-a Visulize TF motifs on the Umap
library(Matrix.utils)
library(chromVAR)
library(SummarizedExperiment)
library(Matrix)
library(motifmatchr)
library(universalmotif)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
Islet.scATAC.metadata<-read.table("Islet.scATAC.metadata")
row.names(Islet.scATAC.metadata)<-Islet.scATAC.metadata$Name
hocomoco.motifsList<-readRDS("hocomoco.motifsList")
Islet12ALL.coembed.Visulize.meta.atac.filtered<-readRDS("Islet12ALL.coembed.Visulize.meta.atac.filtered")
SVM.KNN.filtered.final.print<-read.table("SVM.KNN.filtered.dic.Final.tab")
names(SVM.KNN.filtered.final.print)<-c("Name","CellType")
row.names(SVM.KNN.filtered.final.print)<-SVM.KNN.filtered.final.print$Name
Islet.scATAC.dac.IsletALLPeak.reproducible<-read.table("Islet.scATAC.dac.IsletALLPeak.reproducible")
names(Islet.scATAC.dac.IsletALLPeak.reproducible)<-c("Cell","Peak","Cts")
ALLPeak.reproducible.Mtx<-dMcast(Islet.scATAC.dac.IsletALLPeak.reproducible,Peak~Cell)
colnames(ALLPeak.reproducible.Mtx)<-gsub("Cell","",colnames(ALLPeak.reproducible.Mtx))
ALLPeak.reproducible.filtered.Mtx<-ALLPeak.reproducible.Mtx[,colnames(ALLPeak.reproducible.Mtx) %in% row.names(Islet12ALL.coembed.Visulize.meta.atac.filtered)]
ALLPeak.reproducible.filtered.Mtx<-ALLPeak.reproducible.filtered.Mtx[rowSums(ALLPeak.reproducible.filtered.Mtx)!=0,]

Chr<-strsplit(row.names(ALLPeak.reproducible.filtered.Mtx),"_") %>% lapply(.,function(x){x[1]}) %>% unlist
Start<-strsplit(row.names(ALLPeak.reproducible.filtered.Mtx),"_") %>% lapply(.,function(x){x[2]}) %>% unlist %>% as.numeric
End<-strsplit(row.names(ALLPeak.reproducible.filtered.Mtx),"_") %>% lapply(.,function(x){x[3]}) %>% unlist %>% as.numeric
peaks<-GRanges(Chr,strand="*",ranges=IRanges(start=Start,end=End))
fragment_counts <- SummarizedExperiment(assays = list(counts = ALLPeak.reproducible.filtered.Mtx),rowRanges = peaks)
# fragment_counts@colData<-Islet.scATAC.metadata[SVM.KNN.filtered.final.print$Name,] %>% DataFrame
fragment_counts@colData<-Islet.scATAC.metadata[row.names(Islet12ALL.coembed.Visulize.meta.atac.filtered),] %>% DataFrame
names(fragment_counts@colData)[2]<-"depth"
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg19)
counts_filtered <- filterPeaks(fragment_counts)
row.names(counts_filtered@colData)<-counts_filtered@colData$Name
motif_ix <- matchMotifs(hocomoco.motifsList, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg19)
bg <- getBackgroundPeaks(object = counts_filtered)
dev <- computeDeviations(object = counts_filtered, annotations = motif_ix,background_peaks = bg)
MotifCellTypeList<-list()
MotifCellTypeList.summary<-list()
for (ref in c("Alpha","Acinar","Beta","Delta","Duct","PP","PSC"))
{
print(ref)
pvalues<-c()
slopes<-c()
for(i in 1:nrow(dev@assays@data$z)){
data<-Tomerge_v2(data.frame(z=dev@assays@data$z[i,]),SVM.KNN.filtered.final.print)
data$CellType<-relevel(data$CellType,ref=ref)
md<-lm(z~CellType,data=data)
slopes<-rbind(slopes,md$coefficients[2:length(md$coefficients)])
pvalues<-rbind(pvalues,summary(md)$coefficient[,4][2:length(md$coefficients)])
}
row.names(pvalues)<-row.names(dev@assays@data$z)
row.names(slopes)<-row.names(dev@assays@data$z)
MotifCellTypeList<-c(MotifCellTypeList,list(list(slopes=slopes,pvalues=pvalues)))
summary<-data.frame(pvalues=apply(pvalues,1,function(x){exp(mean(log(x)))}),slopes=rowMeans(slopes))
summary<-summary[order(summary$slopes),]
MotifCellTypeList.summary<-c(MotifCellTypeList.summary,list(summary))
}
names(MotifCellTypeList)<-c("Alpha","Acinar","Beta","Delta","Duct","PP","PSC")
names(MotifCellTypeList.summary)<-c("Alpha","Acinar","Beta","Delta","Duct","PP","PSC")

tocheck.alpha<-MotifCellTypeList.summary$Alpha %>% head(.,n=10) %>% row.names
tocheck.beta<-MotifCellTypeList.summary$Beta %>% head(.,n=10) %>% row.names
tocheck.delta<-MotifCellTypeList.summary$Delta %>% head(.,n=10) %>% row.names
tocheck.PP<-MotifCellTypeList.summary$PP %>% head(.,n=10) %>% row.names
tocheck.Duct<-MotifCellTypeList.summary$Duct %>% head(.,n=10) %>% row.names
tocheck.acinar<-MotifCellTypeList.summary$Acinar %>% head(.,n=10) %>% row.names
tocheck.PSC<-MotifCellTypeList.summary$PSC %>% head(.,n=10) %>% row.names
allpss<-list()
for(tocheck in list(tocheck.alpha,tocheck.beta,tocheck.delta,tocheck.PP,tocheck.Duct,tocheck.acinar,tocheck.PSC)){
motiftoplot<-dev@assays@data$z[tocheck,] %>% t %>% as.data.frame %>% Tomerge_v2(.,counts_filtered@colData)
motiftoplot<-Tomerge_v2(Islet12ALL.coembed.Visulize.meta.atac.filtered[,c("UMAP_1","UMAP_2")],motiftoplot)
ps<-list()
for(mtf in names(motiftoplot)[3:12]){
p<-ggplot(motiftoplot)+aes_string("UMAP_1","UMAP_2",color=mtf)+geom_point(size=0.1)+scale_color_gradient(low="lightsteelblue1",high="darkred",limit=c(-1,2),oob=scales::squish)+theme_classic()
ps<-c(ps,list(p))
}
allpss<-c(allpss,list(ps))
}
allpss<-readRDS("/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/allpss")
motiftoplot<-dev@assays@data$z[c("NEUROD1.MA1109.1","PDX1_HUMAN.H11MO.1.A","IRF1_HUMAN.H11MO.0.A","FOXA2_HUMAN.H11MO.0.A","EHF_HUMAN.H11MO.0.B","ELF3_HUMAN.H11MO.0.A","RUNX3_HUMAN.H11MO.0.A"),] %>% t %>% as.data.frame %>% Tomerge_v2(.,counts_filtered@colData)
motiftoplot<-Tomerge_v2(Islet12ALL.coembed.Visulize.meta.atac.filtered[,c("UMAP_1","UMAP_2")],motiftoplot)

p.neurod1<-ggplot(motiftoplot)+aes_string("UMAP_1","UMAP_2",color="NEUROD1.MA1109.1")+geom_point(size=1.8,alpha=0.5)+scale_color_gradient(low="lightsteelblue1",high="darkred",limit=c(-1.5,2),oob=scales::squish)+theme_void()
p.pdx1<-ggplot(motiftoplot)+aes_string("UMAP_1","UMAP_2",color="PDX1_HUMAN.H11MO.1.A")+geom_point(size=1.8,alpha=0.5)+scale_color_gradient(low="lightsteelblue1",high="darkred",limit=c(-0.5,2.8),oob=scales::squish)+theme_void()
p.EHF<-ggplot(motiftoplot)+aes_string("UMAP_1","UMAP_2",color="EHF_HUMAN.H11MO.0.B")+geom_point(size=1.8,alpha=0.5)+scale_color_gradient(low="lightsteelblue1",high="darkred",limit=c(0.5,2),oob=scales::squish)+theme_void()

p.RUNX<-ggplot(motiftoplot)+aes_string("UMAP_1","UMAP_2",color="RUNX3_HUMAN.H11MO.0.A")+geom_point(size=1.8,alpha=0.5)+scale_color_gradient(low="lightsteelblue1",high="darkred",limit=c(0.5,2),oob=scales::squish)+theme_void()
setwd("~/Chip-seq/ATAC/scATAC/IsletscATAC/Fig1.Clustering")
png("TF.umap.NEUROD.png",height=1650,width=2300,res=300)
print(p.neurod1)
dev.off()
png("TF.umap.PDX1.png",height=1650,width=2300,res=300)
print(p.pdx1)
dev.off()
png("TF.umap.EHF.png",height=1650,width=2300,res=300)
print(p.EHF)
dev.off()
png("TF.umap.RUNX.png",height=1650,width=2300,res=300)
print(p.RUNX)
dev.off()

  ##region 2020-8-28-c Visulize single-cell track examples
setwd("~/Chip-seq/ATAC/scATAC/IsletscATAC")
Islet.scATAC.dac.Bin5K.Mtx<-readRDS("Islet.scATAC.dac.Bin5K.Mtx")
SVM.KNN.filtered.final<-readRDS("SVM.KNN.filtered.final")
GenomicRange.Index<-read.table("hg19.genome_split_5000")
Loci.df<-data.frame(Celltype=c("Beta","Alpha","Delta","PP","Duct","PSC","Acinar"),Loci=c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A"),Chr=c("chr11","chr2","chr3","chr17","chr17","chr7","chr2"),Starts=c(2143641,162976712,187361895,41988356,39667409,94015400,79334255),Ends=c(2209812,163029039,187416856,42053853,39695067,94065961,79363874))
Loci.bed<-MakeBinBed(df=Loci.df,res=100)
write.table(Loci.bed,"Islet12_7celltype.marker.100bp.bed",row.names=F,col.names=F,sep="\t",quote=F)
## To generate signal in bins in bash
"
bedtools intersect -a Islet12_7celltype.marker.100bp.bed -b Islet12.ALL.Monoclonal.BC.sort.bed -wa -wb -loj > Islet12_7celltype.marker.100bp.monoclonal.beds
"
## Readin the overlapped bedfile and modify
Islet12_7celltype.marker.100bp.monoclonal.beds<-read.table("Islet12_7celltype.marker.100bp.monoclonal.beds")
Islet12_7celltype.marker.100bp.monoclonal<-data.frame(Pixel=Islet12_7celltype.marker.100bp.monoclonal.beds$V4,Cell=Islet12_7celltype.marker.100bp.monoclonal.beds$V8,Signal=1)
Islet12_7celltype.marker.100bp.monoclonal.Mtx<-dMcast(Islet12_7celltype.marker.100bp.monoclonal,Pixel~Cell)
colnames(Islet12_7celltype.marker.100bp.monoclonal.Mtx)<-gsub("Cell","",colnames(Islet12_7celltype.marker.100bp.monoclonal.Mtx))
## To prepare for plotting
colors<-c(rgb(51,153,102,maxColorValue=255),rgb(0,102,204,maxColorValue=255),rgb(153,51,255,maxColorValue=255),rgb(204,102,0,maxColorValue=255),rgb(255,102,0,maxColorValue=255),rgb(255,153,204,maxColorValue=255),rgb(255,0,0,maxColorValue=255))
Celltypes<-c("Beta","Alpha","Delta","PP","Duct","PSC","Acinar")
ps<-list()
for(i in 1:7){
Celltype.focus=Celltypes[i]
TopN=200
TopCells<-subset(SVM.KNN.filtered.final,p1==Celltype.focus) %>% .[order(.$Fragments,decreasing=T),] %>% head(.,n=TopN) %>% row.names
Selected.MTX<-Islet12_7celltype.marker.100bp.monoclonal.Mtx[,colnames(Islet12_7celltype.marker.100bp.monoclonal.Mtx) %in% TopCells]
Selected.MTX.long<-melt(as.matrix(Selected.MTX))
Selected.MTX.long<-strsplit(as.character(Selected.MTX.long$Var1),"_") %>% do.call(rbind,.)  %>% cbind(Selected.MTX.long,.)
names(Selected.MTX.long)<-c("Binname","Cell","Signal","Loci","BinN")
Selected.MTX.long$BinN<-as.numeric(as.character(Selected.MTX.long$BinN))
Selected.MTX.long$Signal[Selected.MTX.long$Signal>1]<-1
Selected.MTX.long$Signal<-as.factor(Selected.MTX.long$Signal)
Selected.MTX.long$Loci<-factor(Selected.MTX.long$Loci,levels=c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A"))
p<-ggplot(Selected.MTX.long)+aes(BinN,Cell,fill=Signal)+geom_tile()+facet_grid(~Loci,scales="free")+scale_fill_manual(values=c("white",colors[i]))+theme_void()+theme(panel.border=element_rect(size=0.2,color="grey90",fill=NA),legend.position="null")
ps<-c(ps,list(p))
}
## To plot
setwd("~/Chip-seq/ATAC/scATAC/IsletscATAC/Fig1.Clustering")
png("Beta.scTrack.png",height=600,width=5000,res=300)
print(ps[[1]])
dev.off()
png("Alpha.scTrack.png",height=600,width=5000,res=300)
print(ps[[2]])
dev.off()
png("Delta.scTrack.png",height=600,width=5000,res=300)
print(ps[[3]])
dev.off()
png("PP.scTrack.png",height=600,width=5000,res=300)
print(ps[[4]])
dev.off()
png("Duct.scTrack.png",height=600,width=5000,res=300)
print(ps[[5]])
dev.off()
png("PSC.scTrack.png",height=600,width=5000,res=300)
print(ps[[6]])
dev.off()
png("Acinar.scTrack.png",height=600,width=5000,res=300)
print(ps[[7]])
dev.off()
##endregion

  ##region 2020-8-28-c Visulize signature gene RNA expression
Islet12ALL.coembed.Visulize<-readRDS("/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/Islet12ALL.coembed.Visulize")
Islet12.scRNA.seurat3.filtered<-readRDS("/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/Islet12.scRNA.seurat3.filtered")

RNASig.tiplot<-t(Islet12.scRNA.seurat3.filtered@assays$RNA@scale.data[c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A","FLT1","ANGPT2"),]) %>% Tomerge_v2(.,Islet12ALL.coembed.Visulize[[1]]@reductions$umap@cell.embeddings) %>% .[complete.cases(.),]
setwd("~/Chip-seq/ATAC/scATAC/IsletscATAC/Fig1.Clustering")
for(gene in c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A","FLT1")){
png(paste(gene,"RNA.png",sep="."),height=2000,width=2100,res=300)
p<-ggplot(RNASig.tiplot)+aes_string("UMAP_1","UMAP_2",color=gene)+geom_point()+scale_color_gradient(low="azure3",high="red")+theme_void()
print(p)
dev.off()
}
##endregion
