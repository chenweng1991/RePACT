#' GetRePACTmodel.cca
#'
#' This function is to run logistic regression based on the number of LSIs, and characteristics of samples.
#' @param ccaWithinfo, phenotable
#' @param PCrange, Specified PCs for regression
#' @param prefix
#' @param pheno, the column name of the "characteristics to compare" in the phenodic.use dataframe
#' @param CCrange
#' @return The function return logistic model
#' @import Seurat pscl
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
