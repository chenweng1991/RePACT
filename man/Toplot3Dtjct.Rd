% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RePACT.R
\name{Toplot3Dtjct}
\alias{Toplot3Dtjct}
\title{Toplot3Dtjct}
\usage{
Toplot3Dtjct(trajectory.ob = tjct.ob, PCrange = 1:3, pheno,
  linear = T, fam = "binomial", enlag = 10, decided = F,
  theta = 60, phi = 180, singeplotname = NULL, multiplotname,
  titlename = "PC1-3 BMI trajectory regression\\n")
}
\arguments{
\item{trajectory.ob}{This is the object generated from Prepareforpseudoregress.g}

\item{PCrange}{The PC dimensions used for 3D plotting,  usually pick top 3 significant pCs}

\item{pheno}{the column name for regression}

\item{linear}{If true, perform a linear regression, else, perform a logistic regression}

\item{fam}{For logistic regression, "binomial" is used for glm}

\item{enlag}{this is to adjust the length of regression line, default is 10}

\item{decided}{If TRUE, then only plot one based on angles by theta, and phi, otherwise, plot a series of figures in one big PDF.}

\item{theta}{A number indicating angle, I will work if decided}

\item{phi}{A another number indicating angle, I will work if decided}

\item{singeplotname}{a pdf file name if I have a decided single pdf}

\item{multiplotname}{a pdf file name if generating a series of figures}

\item{titlename}{This is the name in figure title}

\item{theta}{A number indicating angle, I will work if decided}
}
\value{
The output include: PCvariance this summrize the percentage of variance chosen PC can explain; PCanfpheno this is a data frame including PC information, phenotype information as well as pseudo.index and residues after regression; $object.raw.withinfo: this is a datafrom of raw data with pseudo.index and residues; $model: this is the model. reg.plot.2d this is the 2d regression plot. model.para: this is the regression parameter.
}
\description{
This function is to make 3D example plot for regression analysis
}
\examples{
Toplot3Dtjct(T2D.tjct.ob,PCrange=c(1,3,4),pheno="disease",linear=F,multiplotname="test.pdf",titlename="tittle")
}
