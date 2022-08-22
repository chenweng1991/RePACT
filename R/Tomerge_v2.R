
Tomerge_v2 <- 
function(A,B,leavex=T,leavey=F){
        mergeAB<-merge(A,B,by="row.names",all.x=leavex,all.y=leavey)
        row.names(mergeAB)<-mergeAB[,1]
        mergeAB<-mergeAB[,-1]
        mergeAB<-mergeAB[rownames(A),]
        return(mergeAB)
}

