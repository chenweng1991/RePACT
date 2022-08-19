#' CallT2Dpeak_qvalue
#'
#' This function is to run regression (linear/logistic) based on the number of PCs, and characteristics of samples.
#' @param cellvsPeak.m.aggr, A scRNA Seurat (V3/V4), with PCs and clustering info
#' @param depths
#' @param index
#' @param qcut
#' @param slopecut1
#' @param slopecut2
#' @param doscale
#' @return The function return a list: "PCvariance", "PCanfpheno", "object.raw.withinfo", "model", "reg.plot.2d", "model.para"
#' @import Seurat ggplot2 Matrix RColorBrewer gridExtra pscl
#' @export
#' @examples
#' RepACT.obj<-Prepareforpseudoregress.g(OBJ,PCrange=PCrange,phenodic.use=sub_phenotable,pheno=pheno,linear=T)

CallT2Dpeak_qvalue <- function(cellvsPeak.m.aggr=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$cellvsPeak.m.aggr, depths=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$depths, index=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$index,qcut=0.1,slopecut1=0.5,slopecut2=-0.5,doscale=T){
cellvsPeak.m.aggr.norm<-t(t(cellvsPeak.m.aggr)/depths)
cellvsPeak.m.aggr.norm.scale<- apply(cellvsPeak.m.aggr.norm,1,function(x){(x-mean(x))/sd(x)}) %>% t
# cellvsPeak.m.aggr.scale<- apply(cellvsPeak.m.aggr,1,function(x){(x-mean(x))/sd(x)}) %>% t
index.scale=scale(index)
pseudoregress.all<-c()
for (i in 1:nrow(cellvsPeak.m.aggr))     #Loop gene by gene
{
        if(any(is.na(as.numeric(cellvsPeak.m.aggr.norm.scale[i,])))){
                pseudoregress.all<-rbind(pseudoregress.all,c(NA,NA,NA,NA,NA))
                next
        }
        if(doscale){
                gmodel1<-glm(as.numeric(cellvsPeak.m.aggr.norm.scale[i,] )~index.scale,family=gaussian)
        }       else{
 gmodel1<-glm(as.numeric(cellvsPeak.m.aggr.norm[i,] )~index,family=gaussian)
        }
 cur.slope1<-summary(gmodel1)$coefficients[,1][2]   # to get slope of pseudoBMIindex
 cur.pvalues1<-summary(gmodel1)$coefficients[,4][2]  # to get p value for pseudoBMIindex
 cur.corr<-cor(1e4*as.numeric(cellvsPeak.m.aggr.norm.scale[i,]),index.scale)
 cur.MaxMedian.FC<-max(cellvsPeak.m.aggr.norm[i,])/median(cellvsPeak.m.aggr.norm[i,])
 cur.Max<-1e6*max(cellvsPeak.m.aggr.norm[i,])
 names(cur.slope1)<-"slope"
 names(cur.pvalues1)<-"pvalue"
 names(cur.corr)<-"cor"
 names(cur.MaxMedian.FC)<-"MaxMedianFC"
 names(cur.Max)<-"MaxRPKM"
#   names(meanUMI)<-"meanUMI"
 pseudoregress.all<-rbind(pseudoregress.all,c(cur.slope1,cur.pvalues1,cur.corr,cur.MaxMedian.FC,cur.Max))
 print(i)
}
row.names(pseudoregress.all)<-row.names(cellvsPeak.m.aggr)
pseudoregress.all<-as.data.frame(pseudoregress.all)
pseudoregress.all$qvalue<-qvalue(pseudoregress.all$pvalue)$qvalues

UP<-subset(pseudoregress.all,qvalue< qcut &slope >slopecut1 )  %>% .[order(.$slope,decreasing=T),]
DN<-subset(pseudoregress.all,qvalue< qcut &slope <slopecut2 )  %>% .[order(.$slope,decreasing=F),]
UPDN.toplot<-t(cellvsPeak.m.aggr)/depths
return(list(pseudoregress.all=pseudoregress.all,UP=UP,DN=DN,UPDN.toplot=UPDN.toplot))
}
