  ##region  Prepare scRNAseq data
  Islet.scATAC.metadata<-read.table("Islet.scATAC.metadata")
  row.names(Islet.scATAC.metadata)<-Islet.scATAC.metadata$Name
  Islet.scRNA.seurat3<-readRDS("All_islet_30.OBJ3.subset.with_cell_type.rds")
  Islet.scRNA.seurat3@meta.data$Cell_type<-mapvalues(Islet.scRNA.seurat3@meta.data$Cell_type,from=c("Acinar","Alpha","Alpha_proliferating","Beta","Delta","Duct","Endothelial","Epsilon","PP","PSC_active","PSC_proliferating","PSC_quiescent"),to=c("Acinar","Alpha","Alpha","Beta","Delta","Duct","Endothelial","Epsilon","PP","PSC","PSC","PSC"))
  Secondary_Filtering.panel<-Islet.scRNA.seurat3@assays$RNA@counts[c("INS","GCG","PPY","SST","COL1A2","KRT19"),] %>% t() %>% as.data.frame() %>% Tomerge_v2(.,Islet.scRNA.seurat3@meta.data[,"Cell_type",drop=F])
  INScut<-30  #10%
  GCGcut<-20 #10%
  PPYcut<-33 #10%
  SSTcut<-57 #10%
  COLcut<-2 #10%
  KRTcut<-1 #10%
  tofilter.names.alpha<-row.names(subset(Secondary_Filtering.panel,Cell_type=="Alpha"  & (GCG<GCGcut |PPY>PPYcut | SST >SSTcut | INS>INScut |COL1A2>COLcut| KRT19>KRTcut )))
  tofilter.names.beta<-row.names(subset(Secondary_Filtering.panel,Cell_type=="Beta"  & (INS<INScut |PPY>PPYcut | SST >SSTcut | GCG>GCGcut |COL1A2>COLcut| KRT19>KRTcut)))
  tofilter.names.delta<-row.names(subset(Secondary_Filtering.panel,Cell_type=="Delta"  & (SST<SSTcut |PPY>PPYcut | INS >INScut | GCG>GCGcut |COL1A2>COLcut| KRT19>KRTcut)))
  tofilter.names.pp<-row.names(subset(Secondary_Filtering.panel,Cell_type=="PP"  & (PPY<PPYcut | INS>INScut | SST >SSTcut | GCG>GCGcut |COL1A2>COLcut| KRT19>KRTcut)))
  tofilter.names.duct<-row.names(subset(Secondary_Filtering.panel,Cell_type=="Duct" & (KRT19<KRTcut| INS>INScut | SST >SSTcut | GCG>GCGcut | PPY >PPYcut | COL1A2>COLcut)))
  tofilter.names.psc<-row.names(subset(Secondary_Filtering.panel,Cell_type=="PSC" & (COL1A2< COLcut |INS>INScut | SST >SSTcut | GCG>GCGcut | PPY >PPYcut | KRT19>KRTcut)))
  tofilter.names.acinarPSC<-row.names(subset(Secondary_Filtering.panel,Cell_type %in% c("Acinar","Endothelial")  & (INS>INScut | SST >SSTcut | GCG>GCGcut | PPY >PPYcut)))
  tofilter.names<-c(tofilter.names.alpha,tofilter.names.beta,tofilter.names.delta,tofilter.names.pp,tofilter.names.duct,tofilter.names.psc,tofilter.names.acinarPSC)

  ##region Plot for supplementary fig1
  tofilter.names.INS<-row.names(subset(Secondary_Filtering.panel, Cell_type!="Beta" &INS>INScut))
  tofilter.names.GCG<-row.names(subset(Secondary_Filtering.panel, Cell_type!="Alpha" &GCG>GCGcut))
  tofilter.names.PPY<-row.names(subset(Secondary_Filtering.panel, Cell_type!="PP" &PPY>PPYcut))
  tofilter.names.SST<-row.names(subset(Secondary_Filtering.panel, Cell_type!="Delta" &SST>SSTcut))
  tofilter.names.toplot<-list(tofilter.names.INS,tofilter.names.GCG,tofilter.names.SST,tofilter.names.PPY)
ps<-list()
for(i in 1:4){
gene<-c("INS","GCG","SST","PPY")[i]
FT<-tofilter.names.toplot[[i]]
datatoplot<-ifelse(row.names(Secondary_Filtering.panel) %in% FT,"Remove","") %>% cbind(Secondary_Filtering.panel,tag=.)
datatoplot$Cell_type<-factor(datatoplot$Cell_type,levels=c("Beta","Alpha","Delta","PP","Duct","Acinar","PSC","Endothelial"))
datatoplot<-datatoplot[complete.cases(datatoplot),]
p<-ggplot(datatoplot)+aes_string("Cell_type",gene,color="tag",size="tag")+geom_jitter()+scale_color_manual(values=c("black","red"))+scale_size_manual(values=c(0.1,0.6))#+theme_void()
ps<-c(ps,list(p))
}

  png("Secondary.Filter.png",height=2000,width=2500,res=300)
  grid.arrange(grobs=ps)
  dev.off()

  Islet.scRNA.seurat3.meta.filtered<-Islet.scRNA.seurat3@meta.data[!row.names(Islet.scRNA.seurat3@meta.data) %in% tofilter.names,] %>% subset(.,nCount_RNA>=1000)
  Islet.Endo.scRNA.seurat3.meta.filtered<-subset(Islet.scRNA.seurat3.meta.filtered,Cell_type %in% c("Alpha","Beta","Delta","PP"))
  # Clustering for all cells
  Islet12.scRNA.seurat3.filtered<-docluster(Islet.scRNA.seurat3@assays$RNA@counts[,row.names(Islet.scRNA.seurat3.meta.filtered)],meta=Islet.scRNA.seurat3.meta.filtered,reso=0.6,nGene=2000,mincell=3,minfeature=200)
  # Clustering for endocrine cells
  names(Islet.Endo.scRNA.seurat3.meta.filtered)[5]<-"Donor"
  Islet12.Endo.scRNA.seurat3.filtered<-docluster(Islet.scRNA.seurat3@assays$RNA@counts[,row.names(Islet.Endo.scRNA.seurat3.meta.filtered)],meta=Islet.Endo.scRNA.seurat3.meta.filtered[,c("Donor","Cell_type")],reso=0.6,nGene=1000,mincell=3,minfeature=200)
  FindNeighbors(Islet12.scRNA.seurat3.filtered)
  Islet12.Endo.scRNA.seurat3.filtered <- FindNeighbors(Islet12.Endo.scRNA.seurat3.filtered, dims = 1:20)
  Islet12.Endo.scRNA.seurat3.filtered <- FindClusters(Islet12.Endo.scRNA.seurat3.filtered, resolution = 0.2)
  Islet12.Endo.scRNA.seurat3.filtered <- RunUMAP(Islet12.Endo.scRNA.seurat3.filtered, dims = 1:20)
  # Third filtering for endo
  ThirdFilter<-subset(Islet12.Endo.scRNA.seurat3.filtered@meta.data,(Cell_type=="Alpha" & !RNA_snn_res.0.2 %in% c(1,2,3)) |(Cell_type=="Beta" & !RNA_snn_res.0.2 %in% c(0)) | (Cell_type=="PP" & !RNA_snn_res.0.2 %in% c(4))|(Cell_type=="Delta" & !RNA_snn_res.0.2 %in% c(5))) %>% row.names
  CellsToKeep<-row.names(Islet12.Endo.scRNA.seurat3.filtered@meta.data)[!row.names(Islet12.Endo.scRNA.seurat3.filtered@meta.data) %in% ThirdFilter]
  Islet12.Endo.scRNA.seurat3.filtered<-subset(Islet12.Endo.scRNA.seurat3.filtered,cells=CellsToKeep)
  # Third filtering also filter for all cells
  CellsToKeep<-row.names(Islet12.scRNA.seurat3.filtered@meta.data)[!row.names(Islet12.scRNA.seurat3.filtered@meta.data) %in% ThirdFilter]
  Islet12.scRNA.seurat3.filtered<-subset(Islet12.scRNA.seurat3.filtered,cells=CellsToKeep)
  #
  Islet12.scRNA.seurat3.filtered@meta.data$Donor<-dic[Islet12.scRNA.seurat3.filtered@meta.data$Sample,]
  Islet12.Endo.scRNA.seurat3.filtered@meta.data$DonorName<-dic[Islet12.Endo.scRNA.seurat3.filtered@meta.data$Donor,]
  # Tovisulize the clustering
  Tomerge_v2(Islet12.scRNA.seurat3.filtered@reductions$umap@cell.embeddings,Islet12.scRNA.seurat3.filtered@meta.data) %>% ggplot(.)+aes(UMAP_1,UMAP_2,color=Cell_type)+geom_point()
  Tomerge_v2(Islet12.Endo.scRNA.seurat3.filtered@reductions$umap@cell.embeddings,Islet12.Endo.scRNA.seurat3.filtered@meta.data) %>% ggplot(.)+aes(UMAP_1,UMAP_2,color=Cell_type)+geom_point()
  Islet12.scRNA.seurat3.filtered<-readRDS("/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/Islet12.scRNA.seurat3.filtered")
  Islet12.Endo.scRNA.seurat3.filtered<-readRDS("/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/RDS/Islet12.Endo.scRNA.seurat3.filtered")
  CelltypeCopo<-as.data.frame(table(Islet12.scRNA.seurat3.filtered@meta.data$Cell_type))
  CelltypeCopo$ratio<-CelltypeCopo$Freq/sum(CelltypeCopo$Freq)
