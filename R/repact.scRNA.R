#' repact.scRNA.R
#'
#' This function is to run repact on a Seurat object.
#' @param OBJ, A scRNA Seurat (V3/V4), with PCs and clustering info, OBJ@meta.data should contain "Sample"(donor) column and the column of the chracteristics to compare, e.g. if_T2D, BMI
#' @param PCrange, Specified PCs for regression
#' @param pheno, the column name of the "characteristics to compare" in the OBJ@meta.data
#' @param is_continuous, if the "characteristics to compare" is continuous
#' @param norm_index, if the pseudo index should be normalized
#' @param SlopeCut, the cutoff of slope to define differentially genes across pseudo index (n=binnumber)
#' @param output_name2, a path/outputname to write out the results to disk, including .rds, .csv, .pdf
#' @return .rds, .csv, .pdf to disk
#' @import Seurat ggplot2 Matrix RColorBrewer gridExtra pscl dplyr
#' @export
#' @examples
#' repact.scRNA.R(OBJ, PCrange=1:20,pheno="diseaseStat",is_continuous="F",norm_index=T,SlopeCut=0.05,output_name2="T2D_Beta.scRNA.RePACT")

repact.scRNA <- function(OBJ,PCrange=1:20,pheno,is_continuous,norm_index=T,SlopeCut=0,output_name2){
  require(Seurat)
  require(ggplot2)
  require(Matrix)
  require(RColorBrewer)
  require(dplyr)
  require(pscl)
  require(gridExtra)
  phenotable <- OBJ@meta.data
  sub_phenotable = phenotable[complete.cases(phenotable[,pheno]),]
  if(is_continuous == "T"){
    RepACT.obj<-Prepareforpseudoregress.g(OBJ,PCrange=PCrange,phenodic.use=sub_phenotable,pheno=pheno,linear=T)
  }else{
    RepACT.obj<-Prepareforpseudoregress.g(OBJ,PCrange=PCrange,phenodic.use=sub_phenotable,pheno=pheno,linear=F)
  }
  RepACT.2nd.ob<-Tjct.core.gen(RepACT.obj,binnumber=20,norm_index=T,SlopeCut=SlopeCut)
  saveRDS(RepACT.obj, paste(output_name2,"RepactOBJ.1.rds",sep='.'))
  saveRDS(RepACT.2nd.ob, paste(output_name2,"RepactOBJ.2.rds",sep='.'))
  qval=qvalue(RepACT.2nd.ob$BINlinear.result[,2])
  whole_lis=data.frame(RepACT.2nd.ob$BINlinear.result,qval$qvalues)
  write.table(whole_lis,paste(output_name2,"genes.Full_list.csv",sep='.'),quote=FALSE,sep=',')
  write.table(RepACT.2nd.ob$FC.result,paste(output_name2,"genes.FCs.csv",sep='.'),quote=FALSE,sep=',')
  RepACT.obj$PCanfpheno <- RepACT.obj$PCanfpheno[order(RepACT.obj$PCanfpheno[,pheno]),]
  pdf(paste(output_name2,"RePACT.pdf",sep='.'))
  if(is_continuous == "T"){
    if(length(RepACT.obj$model$coefficients)>3){
      print(ggplot(RepACT.obj$PCanfpheno)+aes_string(colnames(RepACT.obj$PCanfpheno)[3],colnames(RepACT.obj$PCanfpheno)[4],color=pheno)+geom_point(size=0.3) + scale_color_gradient(low="blue",high="red")+geom_abline(slope=RepACT.obj$model$coefficients[3]/RepACT.obj$model$coefficients[2]))
    }
    if(length(RepACT.obj$model$coefficients)>5){
      print(ggplot(RepACT.obj$PCanfpheno)+aes_string(colnames(RepACT.obj$PCanfpheno)[5],colnames(RepACT.obj$PCanfpheno)[6],color=pheno)+geom_point(size=0.3) + scale_color_gradient(low="blue",high="red")+geom_abline(slope=RepACT.obj$model$coefficients[5]/RepACT.obj$model$coefficients[4]))
    }
    if(length(RepACT.obj$model$coefficients)>7){
      print(ggplot(RepACT.obj$PCanfpheno)+aes_string(colnames(RepACT.obj$PCanfpheno)[7],colnames(RepACT.obj$PCanfpheno)[8],color=pheno)+geom_point(size=0.3) + scale_color_gradient(low="blue",high="red")+geom_abline(slope=RepACT.obj$model$coefficients[7]/RepACT.obj$model$coefficients[6]))
    }
  }else{
    if(length(RepACT.obj$model$coefficients)>3){
      print(ggplot(RepACT.obj$PCanfpheno)+aes_string(colnames(RepACT.obj$PCanfpheno)[3],colnames(RepACT.obj$PCanfpheno)[4],color=pheno)+geom_point(size=0.3) + guides(color=guide_legend(override.aes=list(size=2)))+scale_color_brewer(palette="Set2")+geom_abline(slope=RepACT.obj$model$coefficients[3]/RepACT.obj$model$coefficients[2]))
    }
    if(length(RepACT.obj$model$coefficients)>5){
      print(ggplot(RepACT.obj$PCanfpheno)+aes_string(colnames(RepACT.obj$PCanfpheno)[5],colnames(RepACT.obj$PCanfpheno)[6],color=pheno)+geom_point(size=0.3) + guides(color=guide_legend(override.aes=list(size=2)))+scale_color_brewer(palette="Set2")+geom_abline(slope=RepACT.obj$model$coefficients[5]/RepACT.obj$model$coefficients[4]))
    }
    if(length(RepACT.obj$model$coefficients)>7){
      print(ggplot(RepACT.obj$PCanfpheno)+aes_string(colnames(RepACT.obj$PCanfpheno)[7],colnames(RepACT.obj$PCanfpheno)[8],color=pheno)+geom_point(size=0.3) + guides(color=guide_legend(override.aes=list(size=2)))+scale_color_brewer(palette="Set2")+geom_abline(slope=RepACT.obj$model$coefficients[7]/RepACT.obj$model$coefficients[6]))
    }
  }
#   if(length(RepACT.obj$model$coefficients)>3){
#     print(ggplot(RepACT.obj$PCanfpheno) + aes_string(colnames(RepACT.obj$PCanfpheno)[3],colnames(RepACT.obj$PCanfpheno)[4],label = "experiment") + geom_point(size=0.1,aes(colour = experiment),show.legend = TRUE) + guides(color=guide_legend(override.aes=list(size=2))))
#   }
#   if(length(RepACT.obj$model$coefficients)>5){
#     print(ggplot(RepACT.obj$PCanfpheno) + aes_string(colnames(RepACT.obj$PCanfpheno)[5],colnames(RepACT.obj$PCanfpheno)[6],label = "experiment") + geom_point(size=0.1,aes(colour = experiment),show.legend = TRUE) + guides(color=guide_legend(override.aes=list(size=2))))
#   }
#   if(length(RepACT.obj$model$coefficients)>7){
#     print(ggplot(RepACT.obj$PCanfpheno) + aes_string(colnames(RepACT.obj$PCanfpheno)[7],colnames(RepACT.obj$PCanfpheno)[8],label = "experiment") + geom_point(size=0.1,aes(colour = experiment),show.legend = TRUE) + guides(color=guide_legend(override.aes=list(size=2))))
#   }
  Tjct.core.plot.ss(RepACT.obj,RepACT.2nd.ob,pheno=pheno,phenodic.use=sub_phenotable,top_gene_num=15,output_name=output_name2,norm_index=T)
  dev.off()
}

