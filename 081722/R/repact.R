#' repact
#'
#' This function is to run repact on a Seurat object.
#' @param OBJ, A scRNA Seurat (V3/V4), with PCs and clustering info
#' @param pheno, the column name of the "characteristics to compare" in the phenodic.use dataframe
#' @param is_continuous, if the "characteristics to compare" is continuous
#' @param output_name2, a path/outputname to write out the results to disk, including .rds, .csv, .pdf
#' @param phenotable, a dataframe contains Sample column(it refers to donor in this study), characteristics to compare(it referes to if the donor/cells are from healthy or T2D, or it can be continuous, e.g. BMI)
#' @param norm_index, if the pseudo index should be normalized
#' @param PCrange, Specified PCs for regression
#' @param SlopeCut, the cutoff of slope to define differentially genes across pseudo index (n=binnumber)
#' @return .rds, .csv, .pdf to disk
#' @import Seurat ggplot2 Matrix RColorBrewer gridExtra pscl
#' @export
#' @examples
#' repact(OBJ,"diseaseStat","F","diseaseStat",phenotable,norm_index=T,SlopeCut=0.05)

repact <- function(OBJ,pheno,is_continuous,output_name2,phenotable,norm_index=T,PCrange=1:20,SlopeCut=0){
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

