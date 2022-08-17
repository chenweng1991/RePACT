Tjct.core.plot.ss <- function(object=NULL,secondobj=NULL,phenodic.use=NULL,pheno=NULL,top_gene_num=10,output_name,norm_index=F){
    do.return=T
    colorset="Set1"
    f3.height=12
    ## Function 1 Do_heatmap
    Do_heatmap<-function(bindata,df1,df2,rankname,top_gene_num=10,hardadd=NULL,last=T,insertinto=0,doreturn=T){
        normalize_01<-function(vector){
                normed<-(vector-min(vector))/(max(vector)-min(vector))
                return(normed)
        }
        require(ggplot2)
        require(gridExtra)
        top_gene_num_DOWN=top_gene_num
        top_gene_num_UP=top_gene_num
        if(nrow(df1)<top_gene_num){
                top_gene_num_UP=nrow(df1)
        }
        if(nrow(df2)<top_gene_num){
                top_gene_num_DOWN=nrow(df2)
        }
        DEgenelist<-c(row.names(df1)[1:top_gene_num_UP],row.names(df2)[1:top_gene_num_DOWN])
        if(last==T){
                DEgenelist.L<-DEgenelist[1:(length(DEgenelist)-insertinto)]
                DEgenelist.R<-setdiff(DEgenelist,DEgenelist.L)
                DEgenelist<-c(DEgenelist.L,hardadd,DEgenelist.R)
        }else{
                DEgenelist.R<-DEgenelist[(insertinto+1):length(DEgenelist)]
                DEgenelist.L<-setdiff(DEgenelist,DEgenelist.R)
                DEgenelist<-c(DEgenelist.L,hardadd,DEgenelist.R)
        }
        DEgenelist<-na.omit(DEgenelist)
        bindata<-bindata[,c(DEgenelist,"tag")]
        bindata<-data.frame(apply(bindata[,-ncol(bindata),drop=F],2,normalize_01),bin=bindata[,ncol(bindata)])
        bindata.m<-reshape2::melt(bindata,id.vars="bin")
        p<-ggplot(bindata.m)+aes(bin,variable,fill=value)+geom_tile()+scale_fill_gradient2(low="white",high="red",mid="orange",midpoint=0.6)+labs(y="Variable genes")
        return(p)
    }
    ## detect input of Tjct.core.plot function
    object$PCanfpheno$Sample <- factor(object$PCanfpheno$Sample, levels=sort(unique(as.character(phenodic.use$Sample))))
    if(norm_index){
            object$PCanfpheno[,"pseudo.index"] <- secondobj$raw.bin[[1]][rownames(object$PCanfpheno),][,"pseudo.index"]
    }
    if(is.null(pheno)){
            print(head(object$PCanfpheno))
            print("please enter a column name for histograme fill")
    }else if(is.null(object)){
            print("please enter parameter like below")
            #print("object=tjct.ob,binnumber=20,qcut=0.05,f3.height=12)
    }else if (is.null(secondobj)){
            print ("please enter secoindary object, make sure it is consistant with the primary trajectory object")
    }else{
            bin.data<-secondobj$bin.data
            BINlinear.result.summarized<-secondobj$BINlinear.result.summarized
            raw.bin<-secondobj$raw.bin
            if(!is.factor(object$PCanfpheno[,pheno])){
                    print("1")
                    p<-ggplot(object$PCanfpheno)+aes_string("Sample","pseudo.index",group="Sample",fill=pheno)+geom_violin()+coord_flip()+scale_fill_gradient(low="white",high="red")+geom_hline(yintercept=raw.bin[[2]],linetype=5,size=0.25)
                    PCandPheno <- object$PCanfpheno %>% .[complete.cases(.),]
                    PCandPheno$Sample <- factor(PCandPheno$Sample, levels=phenodic.use[order(phenodic.use[,pheno]),]$Sample)
                    pp<-ggplot(PCandPheno)+aes_string("Sample","pseudo.index",group="Sample",fill=pheno)+geom_violin()+coord_flip()+scale_fill_gradient(low="white",high="red")+geom_hline(yintercept=raw.bin[[2]],linetype=5,size=0.25)
            }else{
                    print("2")
                    p <- ggplot(object$PCanfpheno)+aes_string("Sample","pseudo.index",group="Sample",fill=pheno)+geom_violin()+coord_flip()+geom_hline(yintercept=raw.bin[[2]],linetype=5,size=0.25)+scale_fill_manual(values=c("blue","red"))
                    PCandPheno <- object$PCanfpheno %>% .[complete.cases(.),]
                    order_list = c()
                    for(ii in unique(PCandPheno[,pheno])){
                            sub_set = PCandPheno[which(PCandPheno[,pheno]==ii),]
                            order_list2 = c()
                            for(samples in unique(sub_set$Sample)){
                                    sub_set2 = sub_set[which(sub_set$Sample==samples),]
                                    avg_subset2 = mean(sub_set2$pseudo.index)
                                    order_list2 = c(order_list2,avg_subset2)
                            }
                            order_list2 = data.frame(order_list2)
                            rownames(order_list2) = unique(sub_set$Sample)
                            order_list2 <- order_list2[order(order_list2),1,drop=F]
                            order_list = c(order_list,rownames(order_list2))
                    }
                    PCandPheno$Sample <- factor(PCandPheno$Sample, levels=order_list)
                    pp<-ggplot(PCandPheno)+aes_string("Sample","pseudo.index",group="Sample",fill=pheno)+geom_violin()+coord_flip()+geom_hline(yintercept=raw.bin[[2]],linetype=5,size=0.25)+scale_fill_manual(values=c("blue","red"))
            }

            write.csv(BINlinear.result.summarized$UP,paste(output_name,"UP_gene.csv",sep=''))
            write.csv(BINlinear.result.summarized$DOWN,paste(output_name,"DOWN_gene.csv",sep=''))
            p1<-ggplot(object$PCanfpheno)+aes(pseudo.index)+geom_histogram(fill="orange")+geom_vline(xintercept=raw.bin[[2]],linetype=5,size=0.25)
            p2<-ggplot(object$PCanfpheno)+aes(pseudo.index,residues,color=Sample)+geom_point(size=0.3)+geom_vline(xintercept=raw.bin[[2]],linetype=5,size=0.25)
            p3<-ggplot(object$PCanfpheno)+aes(pseudo.index,fill=Sample)+geom_histogram(position="fill")+geom_vline(xintercept=raw.bin[[2]],linetype=5,size=0.25)+scale_fill_brewer(palette=colorset)
            p4<-ggplot(object$PCanfpheno)+aes(pseudo.index,fill=Sample)+geom_histogram(position="stack")+geom_vline(xintercept=raw.bin[[2]],linetype=5,size=0.25)+scale_fill_brewer(palette=colorset)
            if(nrow(BINlinear.result.summarized$UP) != 0 | nrow(BINlinear.result.summarized$DOWN) != 0){
                    p5<-Do_heatmap(bin.data,df1=BINlinear.result.summarized$UP,df2=BINlinear.result.summarized$DOWN,rankname="rank",top_gene_num=top_gene_num)
            }
            print(p)
            print(pp)
            if(nrow(BINlinear.result.summarized$UP) != 0 | nrow(BINlinear.result.summarized$DOWN) != 0){
                    print(p5)
            }
    }
}

