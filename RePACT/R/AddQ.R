#' AddQ
#'
#' This function is to label the numberic vector by quantile
#' @param x, a numeric vector
#' @param quantiles
#' @return The function return a vector oif labels
#' @import Seurat plyr dplyr Seurat ggrepel ggplot2
#' @export
#' @examples

AddQ<-function(x,quantiles){
        if(x<=quantiles[2]){
                return("Q1")
        }else if(x<=quantiles[3]){
                return("Q2")
        }else if(x<=quantiles[4]){
                return("Q3")
        }else{
                return("Q4")
        }
}

