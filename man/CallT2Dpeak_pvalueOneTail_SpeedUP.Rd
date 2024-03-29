% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snATAC.RePACT.r
\name{CallT2Dpeak_pvalueOneTail_SpeedUP}
\alias{CallT2Dpeak_pvalueOneTail_SpeedUP}
\title{CallT2Dpeak_pvalueOneTail_SpeedUP}
\usage{
CallT2Dpeak_pvalueOneTail_SpeedUP(
  cellvsPeak.m.aggr,
  depths,
  index,
  doscale = T,
  GlobalSlopes
)
}
\description{
CallT2Dpeak_pvalueOneTail_SpeedUP
}
\examples{
PCAInfo.20bin.donor.PCA <- CallT2Dpeak_pvalueOneTail_SpeedUP(PCAInfo.20bin.donor.ob$cellvsPeak.m.aggr, PCAInfo.20bin.donor.ob$depths, PCAInfo.20bin.donor.ob$index,doscale=T, GlobalSlopes=PCAInfo.20bin.ob.PCA$pseudoregress.all[,1,drop=F])
}
