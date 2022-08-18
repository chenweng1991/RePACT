#' Prepareforpseudoregress.g
#'
#' This function is to run regression (linear/logistic) based on the number of PCs, and characteristics of samples.
#' @param object, A scRNA Seurat (V3/V4), with PCs and clustering info
#' @param PCrange, Specified PCs for regression
#' @param phenodic.use, a dataframe contains Sample column(it refers to donor in this study), characteristics to compare(it referes to if the donor/cells are from healthy or T2D, or it can be continuous, e.g. BMI)
#' @param pheno, the column name of the "characteristics to compare" in the phenodic.use dataframe
#' @param dotsize, the size of dots in the output plots
#' @param linear, if the "characteristics to compare" is binary(e.g. if T2D) or continuous(e.g. BMI)
#' @return The function return a list: "PCvariance", "PCanfpheno", "object.raw.withinfo", "model", "reg.plot.2d", "model.para"
#' @import Seurat ggplot2 Matrix RColorBrewer gridExtra pscl
#' @export
#' @examples
#' RepACT.obj<-Prepareforpseudoregress.g(OBJ,PCrange=PCrange,phenodic.use=sub_phenotable,pheno=pheno,linear=T)
GetRePACTmodel.cca<-function(ccaWithinfo=cca.L2.info, prefix="CC",pheno="Disease",CCrange=1:10){
require(pscl)
CCnames <- paste("ccaWithinfo$",prefix,"_",CCrange, sep = "")
ccaWithinfo[,pheno]<-as.factor(ccaWithinfo[,pheno])
form <- formula(paste("ccaWithinfo[,pheno]", paste(CCnames, collapse = "+"), sep = "~"))
model <- glm(form, family = "binomial")
# p <- ggplot(ccaWithinfo) + aes_string("CC_1", "CC_2", color = pheno) + geom_point() + scale_color_brewer(palette = "Set2") + geom_abline(slope = model$coefficients[3]/model$coefficients[2])
# model.para <- summary(model)
# lrt.pvalue <- 1 - pchisq(summary(model)$null.deviance - summary(model)$deviance, summary(model)$df.null - summary(model)$df.residual)
return(model)
}
