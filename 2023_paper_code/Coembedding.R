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
