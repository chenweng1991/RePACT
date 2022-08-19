#' MakeEvenBinBydepth
#'
#' This function is to run regression (linear/logistic) based on the number of PCs, and characteristics of samples.
#' @param cellvsPeak.m, 
#' @param data.info
#' @param binnumber
#' @return The function return a list: "PCvariance", "PCanfpheno", "object.raw.withinfo", "model", "reg.plot.2d", "model.para"
#' @import Seurat ggplot2 Matrix RColorBrewer gridExtra pscl
#' @export
#' @examples
#' RepACT.obj<-Prepareforpseudoregress.g(OBJ,PCrange=PCrange,phenodic.use=sub_phenotable,pheno=pheno,linear=T)

MakeEvenBinBydepth<-function(cellvsPeak.m=BetacellvsPeak.m.filtered,data.info=BetaPeak.data.info,binnumber=20){
cellvsPeak.m<-cellvsPeak.m[row.names(data.info[order(data.info$rank),]),]
fragment.cutoff<-as.integer(sum(cellvsPeak.m)/binnumber)
step<-as.integer(nrow(data.info)/binnumber)
lastpointer<-1
pointer<-step
while(fragment.cutoff-sum(cellvsPeak.m[lastpointer:pointer,])>1000){
  pointer<-pointer+1
  print(pointer)
}
while(fragment.cutoff-sum(cellvsPeak.m[lastpointer:pointer,])< -5000){
  pointer<-pointer-1
  print(pointer)
}
print(paste("Current pointer at",pointer))
data.info.withbin<-c()
print(paste("total cells:",max(data.info$rank)))

for (i in 1:binnumber){
print(paste(i,lastpointer,"<","ranks","<=",pointer,"total fragments:",sum(cellvsPeak.m[lastpointer:pointer,])))
tmp<-subset(data.info,rank<=pointer &rank>lastpointer) %>% cbind(.,evenfragbin=i)
data.info.withbin<-rbind(data.info.withbin,tmp)
lastpointer<-pointer
pointer<-pointer+step
if(pointer>max(data.info$rank)){
  pointer<-max(data.info$rank)
}
while(fragment.cutoff-sum(cellvsPeak.m[lastpointer:pointer,])>5000){
  pointer<-pointer+1
  if(pointer>max(data.info$rank)){
  break
  }
  print(pointer)
}
if(pointer>max(data.info$rank)){
  break
}
while(fragment.cutoff-sum(cellvsPeak.m[lastpointer:pointer,])< -5000){
  pointer<-pointer-1
  print(pointer)
}
}
Peak.m<-cellvsPeak.m[row.names(data.info[order(data.info$rank),]),]
fragment.cutoff<-as.integer(sum(cellvsPeak.m)/binnumber)
step<-as.integer(nrow(data.info)/binnumber)
lastpointer<-1
pointer<-step
while(fragment.cutoff-sum(cellvsPeak.m[lastpointer:pointer,])>1000){
  pointer<-pointer+1
  print(pointer)
}
while(fragment.cutoff-sum(cellvsPeak.m[lastpointer:pointer,])< -5000){
  pointer<-pointer-1
  print(pointer)
}
print(paste("Current pointer at",pointer))
data.info.withbin<-c()
print(paste("total cells:",max(data.info$rank)))

for (i in 1:binnumber){
print(paste(i,lastpointer,"<","ranks","<=",pointer,"total fragments:",sum(cellvsPeak.m[lastpointer:pointer,])))
tmp<-subset(data.info,rank<=pointer &rank>lastpointer) %>% cbind(.,evenfragbin=i)
data.info.withbin<-rbind(data.info.withbin,tmp)
lastpointer<-pointer
pointer<-pointer+step
if(pointer>max(data.info$rank)){
  pointer<-max(data.info$rank)
}
while(fragment.cutoff-sum(cellvsPeak.m[lastpointer:pointer,])>5000){
  pointer<-pointer+1
  if(pointer>max(data.info$rank)){
  break
  }
  print(pointer)
}
if(pointer>max(data.info$rank)){
  break
}
while(fragment.cutoff-sum(cellvsPeak.m[lastpointer:pointer,])< -5000){
  pointer<-pointer-1
  print(pointer)
}
}

data.info.withbin<-as.data.frame(data.info.withbin)
data.info.withbin$evenfragbin<-as.factor(data.info.withbin$evenfragbin)
cellvsPeak.m.bin<-Tomerge_v2(as.data.frame(cellvsPeak.m),as.data.frame(data.info.withbin[,"evenfragbin",drop=F]))  # Each row a cell, each column a peak, add one more column of bin info
cellvsPeak.m.aggr<-data.frame(name=colnames(cellvsPeak.m.bin)[1:(length(cellvsPeak.m.bin)-1)])
for(i in 1:max(as.numeric(as.character(data.info.withbin$evenfragbin)))){
tmp<-subset(cellvsPeak.m.bin,evenfragbin==i)
tmp.aggregate<-colSums(tmp[,1:(ncol(tmp)-1)])
cellvsPeak.m.aggr<-cbind(cellvsPeak.m.aggr,tmp.aggregate)
names(cellvsPeak.m.aggr)[i+1]<-paste("traj",i,sep="")
print(i)
}
cellvsPeak.m.aggr[,2:ncol(cellvsPeak.m.aggr)] %>% as.data.frame()  -> cellvsPeak.m.aggr
index<-ddply(as.data.frame(data.info.withbin),.(evenfragbin),summarise,mean(pseudo.index.balanced))[,2]
depths<-colSums(cellvsPeak.m.aggr)
return(list(data.info.withbin=data.info.withbin,cellvsPeak.m.aggr=cellvsPeak.m.aggr,index=index,depths=depths))
}
