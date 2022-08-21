#' CallT2Dpeak_pvalueOneTail
#'
#' 
#' @param cellvsPeak.m.aggr
#' @param depths
#' @param index
#' @param doscale
#' @param GlobalSlopes
#' @return 
#' @import 
#' @export
#' @examples

CallT2Dpeak_pvalueOneTail <- function(cellvsPeak.m.aggr=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$cellvsPeak.m.aggr,depths=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$depths,index=ATAConCCA.betaT2D.tjct.3nd.10bin.ob$index,doscale=T,GlobalSlopes){
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
 OneTail.p<-pt(summary(gmodel1)$coefficients[2,3],gmodel1$df.residual,lower.tail=GlobalSlopes[row.names(cellvsPeak.m.aggr.norm.scale)[i],]<0)
 cur.slope1<-summary(gmodel1)$coefficients[,1][2]   # to get slope of pseudoBMIindex
 cur.pvalues1<-OneTail.p  # to get p value for pseudoBMIindex, which is one tailed against global
 cur.corr<-cor(1e4*as.numeric(cellvsPeak.m.aggr.norm.scale[i,]),index.scale)
 cur.MaxMedian.FC<-max(cellvsPeak.m.aggr.norm[i,])/median(cellvsPeak.m.aggr.norm[i,])
 cur.Max<-1e6*max(cellvsPeak.m.aggr.norm[i,])
 names(cur.slope1)<-"slope"
 names(cur.pvalues1)<-"pvalue.onetail"
 names(cur.corr)<-"cor"
 names(cur.MaxMedian.FC)<-"MaxMedianFC"
 names(cur.Max)<-"MaxRPKM"
#   names(meanUMI)<-"meanUMI"
 pseudoregress.all<-rbind(pseudoregress.all,c(cur.slope1,cur.pvalues1,cur.corr,cur.MaxMedian.FC,cur.Max))
 #print(i)
}
row.names(pseudoregress.all)<-row.names(cellvsPeak.m.aggr)
pseudoregress.all<-as.data.frame(pseudoregress.all)
return(pseudoregress.all)
}

