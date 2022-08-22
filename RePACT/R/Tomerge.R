Tomerge <- function(A,B){
mergeAB<-merge(A,B,by="row.names",all=TRUE)
row.names(mergeAB)<-mergeAB[,1]
mergeAB<-mergeAB[,-1]
return(mergeAB)
}

