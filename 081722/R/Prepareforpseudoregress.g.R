Prepareforpseudoregress.g<-function(object=NULL,PCrange=1:10,phenodic.use=NULL,pheno,dotsize=0.3,linear=T){   
    # if downsample==T, I'll do downsampling for donor3'
    require(ggplot2)
    require(RColorBrewer)
    require(gridExtra)
    require(pscl)
    phenodic.use = phenodic.use[complete.cases(phenodic.use[,pheno]),]
    if (is.null(object)){
        print("The parameter to enter is like : object=NULL,PCrange=1:40,phenodic.use=phenodic,pheno,linear=T,fam='binomial'")
        print("object(This is a Seutrat object),PCrange=1:40(decide which PCs to use for regression),phenodic.use=phenodic(This is a phenotype table loaded from outside),pheno(t
his is the column name that is used for regresion analysis),linear=T(if it is F, a non-linear relationship such as T2D/nonT2D that could be logistic using binomial),fam='binomial'(only
used when linear=F)")
    }else{
        ###refer to method 'PCA' of Seurat V1.04
        re.data.use=object@assays$RNA@scale.data
        re.pc.genes=object@assays$RNA@var.features
        re.pc.genes = unique(re.pc.genes[re.pc.genes%in%rownames(re.data.use)])
        re.pc.genes.var = apply(re.data.use[re.pc.genes,],1,var)
        re.pc.genes.use=re.pc.genes[re.pc.genes.var>0]; re.pc.genes.use=re.pc.genes.use[!is.na(re.pc.genes.use)]
        re.pc.data = re.data.use[re.pc.genes.use,]
        re.pca.obj = prcomp(t(re.pc.data))
        ###
        PCvariance<-summary(re.pca.obj)$importance[,40]
        # KeyPCinformation<-data.frame(name=row.names(GetPCAcelldata_v2(object)),GetPCAcelldata_v2(object)[,c(PCrange,which(colnames(GetPCAcelldata_v2(object))=="Sample"),which(colnames(GetPCAcelldata_v2(object))=="orig.ident"))])
        KeyPCinformation<-merge(object@reductions$pca@cell.embeddings[,PCrange], phenodic.use,by=0,all.y=TRUE)
        rownames(KeyPCinformation) <- KeyPCinformation$Row.names
        KeyPCinformation <- KeyPCinformation[,-1]
        # if(!is.null(phenodic.use)){
        #         PCandPheno<-merge(KeyPCinformation,phenodic.use,by="Sample")
        # }else{
        #         PCandPheno<-KeyPCinformation
        # }
        PCandPheno<-KeyPCinformation
        PCnames<-paste("PCandPheno$PC_",PCrange,sep="")
        form<-formula(paste("PCandPheno[,pheno]",paste(PCnames,collapse="+"),sep="~"))
        if(linear==T){
            model<-lm(form)
            PCandPheno[,pheno]<-as.factor(PCandPheno[,pheno])
            PCandPheno[,pheno]=as.numeric(levels(PCandPheno[,pheno]))[PCandPheno[,pheno]]
            p<-ggplot(PCandPheno)+aes_string("PC_1","PC_2",color=pheno)+geom_point(size=dotsize)+scale_color_gradient(low="white",high="red")+geom_abline(slope=model$coefficients[3]/model$coefficients[2])
            PCandPheno<-cbind(PCandPheno,pseudo.index=model$fitted.values,residues=model$residuals)
            PCandPheno[,pheno]<-as.numeric(as.character(PCandPheno[,pheno]))
            model.para<-summary(model)
            model.para<-list(modelsummary=model.para)
        }else{
            PCandPheno[, pheno]<-factor(PCandPheno[, pheno])
            model<-glm(form,family="binomial")
            p<-ggplot(PCandPheno)+aes_string("PC_1","PC_2",color=pheno)+geom_point(size=dotsize) + guides(color=guide_legend(override.aes=list(size=2)))+scale_color_brewer(palette="Set2")+geom_abline(slope=model$coefficients[3]/model$coefficients[2])
            PCandPheno<-unique(cbind(PCandPheno,pseudo.index=model$linear.predictors,residues=model$residuals))
            model.para<-summary(model)
            lrt.pvalue<-1-pchisq(summary(model)$null.deviance-summary(model)$deviance,summary(model)$df.null-summary(model)$df.residual)
            pR2<-pR2(model)[4]
            model.para<-list(modelsummary=model.para,lrt.pvalue=lrt.pvalue,pR2=pR2)
        }
        # row.names(PCandPheno)<-PCandPheno$name
        PCandPheno.names<-row.names(PCandPheno)
        object.raw<-t(as.matrix(object@assays$RNA@counts))[PCandPheno.names,]
        object.raw<-cbind(object.raw,totalreads=apply(object.raw,1,sum))
        object.raw.withinfo<-Tomerge_v2(object.raw,PCandPheno[,c(ncol(PCandPheno)-1,ncol(PCandPheno))])
        return(list(PCvariance=PCvariance,PCanfpheno=PCandPheno,object.raw.withinfo=object.raw.withinfo,model=model,reg.plot.2d=p,model.para=model.para))
    }
}

