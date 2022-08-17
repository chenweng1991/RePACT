#' Tjct.core.gen
#'
#' This function is to spearate cells into the specified number of bins based on even distance on PC space.
#' @param object, output of Prepareforpseudoregress.g function
#' @param binnumber, the number of bins separated based on even distance on PC space
#' @param qcut, the cutoff of qvalue to define differentially genes across pseudo index (n=binnumber)
#' @param norm_index, if the pseudo index should be normalized
#' @param SlopeCut, the cutoff of slope to define differentially genes across pseudo index (n=binnumber)
#' @return a list contains geneUMI.dic, raw.bin, bin.data, binregress, BINlinear.result, BINlinear.result.summarized, parameter, FC.result
#' @import Seurat ggplot2 Matrix RColorBrewer gridExtra pscl qvalue
#' @export
#' @examples
#' Tjct.core.gen(RepACT.obj,binnumber=20,norm_index=T,SlopeCut=SlopeCut)

Tjct.core.gen <- function(object=NULL, binnumber=20, qcut=0.05, norm_index=F, SlopeCut=0){
        mydplyr<-function(df,by="tag",thefunction=mean)
        {
                newdf<-c()
                for (i in as.character(unique(df[,by])))
                {
                        dftm<-df[which(as.character(df[,by])==i),]
                        newdf<-rbind(newdf,apply(dftm,2,thefunction))
                }
                return(newdf)
        }
        DopseuBINregress<-function(rawdata,regressby="pseudo.index",start=1,end=ncol(rawdata)-4){
                require(qvalue)
                pseudoregress.all<-c()
                for (i in start:end)     #Loop gene by gene
                {
                        gmodel1<-glm(rawdata[,i]~rawdata[,regressby],offset=log10(totalreads),data=rawdata,family=gaussian)
                        cur.slope1<-summary(gmodel1)$coefficients[,1][2]   # to get slope of pseudoBMIindex
                        cur.pvalues1<-summary(gmodel1)$coefficients[,4][2]  # to get p value for pseudoBMIindex
                        names(cur.slope1)<-paste("slope.",regressby,sep="")
                        names(cur.pvalues1)<-paste("p.",regressby,sep="")
                        #names(meanUMI)<-"meanUMI"
                        pseudoregress.all<-rbind(pseudoregress.all,c(cur.slope1,cur.pvalues1))
                }
                pseudoregress.all<-as.data.frame(pseudoregress.all)
                row.names(pseudoregress.all)<-names(rawdata)[start:end]
                return(pseudoregress.all)
        }
        givebintag<-function(df,ordername="pseudoBMIindex",bin=10,genenumbercut=50){
                leftendpoint<-quantile(df[,ordername],genenumbercut/nrow(df))
                rightendpoint<-quantile(df[,ordername],1-genenumbercut/nrow(df))
                binunit<-(rightendpoint-leftendpoint)/(bin-2)
                cutoff1<-leftendpoint
                cutoff2<-cutoff1+binunit
                roadmark<-c(cutoff1,cutoff2)
                newdf<-c()
                newdf<-rbind(newdf,cbind(df[which(df[,ordername]<=leftendpoint),],tag=1))
                newdf<-rbind(newdf,cbind(df[which(df[,ordername]>rightendpoint),],tag=bin))
                for (i in 2:(bin-1)){
                        #to skip the bin with 0 cells (generally since too small cells in total but still want to get result)
                        if(nrow(df[which(df[,ordername]>cutoff1 & df[,ordername]<=cutoff2),]) != 0){
                                newdf<-rbind(newdf,cbind(df[which(df[,ordername]>cutoff1 & df[,ordername]<=cutoff2),],tag=i))
                        }
                        cutoff1<-cutoff2
                        cutoff2<-cutoff2+binunit
                        roadmark<-c(roadmark,cutoff2)
                }
                roadmark<-roadmark[-length(roadmark)]
                return(list(newdf,roadmark))
        }
        GetBirank.g<-function(df,cutoff=0.05,Tscut=100,rankby="slope",SlopeCut){  # rankby == "slope" or "qvalue"
                require(qvalue)
                tellquality<-function(vector){   # T2D  pos=3,  BMI pos=4
                        if (vector["qvalue"]>cutoff){
                                return("NS")
                        }
                        else if (vector["slope.pseudo.index"] > SlopeCut){
                                return("UP")
                        }
                        else if (vector["slope.pseudo.index"] < -SlopeCut){
                                return("DOWN")
                        }else{
                                return("NS")
                        }
                }
                df2<-cbind(df,qvalue=qvalue(df[,"p.pseudo.index"])$qvalues)
                ##construct a null data frame
                tmp_df <- df2
                tmp_df$feature = 0
                tmp_df$rank = 0
                null_data <- subset(tmp_df,Transcripts<0)
                ##
                df3<-subset(df2,Transcripts>Tscut & qvalue<cutoff)
                if(nrow(df3) == 0){
                        UP <- null_data
                        DOWN <- null_data
                }else{
                        df4<-cbind(df3,feature=apply(df3,1,tellquality))
                        UP<-subset(df4,feature=="UP")
                        DOWN<-subset(df4,feature=="DOWN")
                        NS<-subset(df4,feature=="NS")
                        if(nrow(UP) == 0 & nrow(DOWN) != 0){
                                if (rankby=="slope"){
                                        DOWN.rank<-data.frame(rank=rank(-abs(DOWN$slope.pseudo.index)),row.names=row.names(DOWN))
                                }else if (rankby=="qvalue"){
                                        DOWN.rank<-data.frame(rank=rank(-abs(DOWN$qvalue)),row.names=row.names(DOWN))
                                }
                                DOWN.rank<-DOWN.rank/max(DOWN.rank)
                                df4<-Tomerge_v2(df4,DOWN.rank)
                                DOWN<-subset(df4,feature=="DOWN")
                                DOWN<-DOWN[order(DOWN$rank),]
                                UP <- null_data
                        }else if(nrow(UP) != 0 & nrow(DOWN) == 0){
                                if (rankby=="slope"){
                                        UP.rank<-data.frame(rank=rank(-abs(UP$slope.pseudo.index)),row.names=row.names(UP))
                                }else if (rankby=="qvalue"){
                                        UP.rank<-data.frame(rank=rank(-abs(UP$qvalue)),row.names=row.names(UP))
                                }
                                UP.rank<-UP.rank/max(UP.rank)
                                df4<-Tomerge_v2(df4,UP.rank)
                                UP<-subset(df4,feature=="UP")
                                UP<-UP[order(UP$rank),]
                                DOWN <- null_data
                        }else{
                                if (rankby=="slope")
                                {
                                        UP.rank<-data.frame(rank=rank(-abs(UP$slope.pseudo.index)),row.names=row.names(UP))
                                        DOWN.rank<-data.frame(rank=rank(-abs(DOWN$slope.pseudo.index)),row.names=row.names(DOWN))
                                }else if (rankby=="qvalue")
                                {
                                        UP.rank<-data.frame(rank=rank(-abs(UP$qvalue)),row.names=row.names(UP))
                                        DOWN.rank<-data.frame(rank=rank(-abs(DOWN$qvalue)),row.names=row.names(DOWN))
                                }

                                UP.rank<-UP.rank/max(UP.rank)
                                DOWN.rank<-DOWN.rank/max(DOWN.rank)
                                rank.dict<-rbind(UP.rank,DOWN.rank)
                                df4<-Tomerge_v2(df4,rank.dict)
                                UP<-subset(df4,feature=="UP")
                                UP<-UP[order(UP$rank),]

                                DOWN<-subset(df4,feature=="DOWN")
                                DOWN<-DOWN[order(DOWN$rank),]
                        }
                }
                return(list(UP=UP,DOWN=DOWN))
        }
        ##claculate FC between two bins, and mean values
        Gene.bin.FC <- function(rawdata,start=1,end=ncol(rawdata)-4){
                fc.all=c()
                for(i in start:end){
                        fc.each_two_bin=c()
                        for(n in 2:nrow(rawdata)){
                                fc.bin <- rawdata[n,i]/rawdata[(n-1),i]
                                fc.each_two_bin <- c(fc.each_two_bin, fc.bin)
                        }
                        fc.each_two_bin.2 = fc.each_two_bin
                        fc.each_two_bin.2[is.infinite(fc.each_two_bin.2)]=NA
                        fc.each_two_bin.2 <- na.omit(fc.each_two_bin.2)
                        if(sum(fc.each_two_bin.2 == 0)){
                                fc.each_two_bin.2 <- fc.each_two_bin.2[-which(fc.each_two_bin.2 == 0)]
                        }
                        fc.each_two_bin <- c(fc.each_two_bin,mean(fc.each_two_bin.2),gm_mean(fc.each_two_bin.2))
                        fc.all <- rbind(fc.all,fc.each_two_bin)
                }
                fc.all<-as.data.frame(fc.all)
                rownames(fc.all)<-names(rawdata)[start:end]
                colnames(fc.all)<-c(1:(nrow(rawdata)-1),"Arithmetic_mean","Geometric_mean")
                return(fc.all)
        }

#---------------
        if(norm_index){
                ###############scale pseudo index to 0~1
                scaled.index <- (object$object.raw.withinfo$pseudo.index - mean(object$object.raw.withinfo$pseudo.index))/sd(object$object.raw.withinfo$pseudo.index)
                object$object.raw.withinfo$pseudo.index = scaled.index
                object$PCanfpheno$pseudo.index = scaled.index
        }
        geneUMI.dic<- data.frame(Transcripts=apply(object$object.raw.withinfo,2,sum)[-((ncol(object$object.raw.withinfo)-2):ncol(object$object.raw.withinfo))])
        # To add Bin tag very everey single cell   #  The result is list with two elements., First is the raw data with info;  Second is the position of making bin
        raw.bin<-givebintag(object$object.raw.withinfo,ordername="pseudo.index",genenumbercut=round(nrow(object$object.raw.withinfo)/binnumber),bin=binnumber)
        print("rawBin has been made")
        # To calculate mean number of every bin for every gene
        bin.data<-as.data.frame(mydplyr(raw.bin[[1]]))
        # Do linear regression for BInTjct.core(alpha.BMI.tjct.ob
        binregress<-DopseuBINregress(bin.data,regressby="pseudo.index")
        # To get a plotable result
        print("regression finished")
        BINlinear.result<-Tomerge(binregress,geneUMI.dic)
        BINlinear.result.summarized<-GetBirank.g(BINlinear.result,rankby="slope",cutoff=qcut,SlopeCut=SlopeCut)
        # To get FCs between bins
        FC.result <- Gene.bin.FC(bin.data)
        # to generate plots
        print("calculation finished")
        parameter<-paste("Take",binnumber,"bins", "the q value to identify sig genes is 1<=",qcut)
        print("in total identified genes are UP=")
        print(row.names(BINlinear.result.summarized$UP))
        print("in total identified genes are DOWN=")
        print(row.names(BINlinear.result.summarized$DOWN))
        return(list(geneUMI.dic=geneUMI.dic,raw.bin=raw.bin,bin.data=bin.data,binregress=binregress,BINlinear.result=BINlinear.result,BINlinear.result.summarized=BINlinear.result.summarized,parameter=parameter,FC.result=FC.result))
}
