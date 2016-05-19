library(data.table)
library(Matrix)
library(plyr)
library(reshape2)
library(rpudplus)

#' ClusterCount
#' 
#' Counts the number of clusters at a specfied height on a dendrogram 
#' generated from an ECFP distance matrix
#' 
#' @param compounds data table of compounds from a specific plate
#' @param cpd.identifier character string identifying compounds: Broad ID, InChIKey, etc...
#' @param height numeric value to cut dendrogram defaults to 0.8
#' @export
#' @examples
#' hashes <- round(rnorm(50, 100000000, 1000), 1)
#' 
#' compounds <- data.table('ECFP_4[1]' = sample(hashes, 100, replace=TRUE),
#'                         'ECFP_4[2]' = sample(hashes, 100, replace=TRUE),
#'                         'ECFP_4[3]' = sample(hashes, 100, replace=TRUE),
#'                         'ECFP_4[4]' = sample(hashes, 100, replace=TRUE),
#'                         'ECFP_4[5]' = sample(hashes, 100, replace=TRUE))
#'                         
#' compounds$InChIKey <- sapply(1:dim(compounds)[1], function(x)
#'                              paste(sample(LETTERS, 10, replace=TRUE), collapse=''))
#'                              
#' ClusterCount(compounds, 'InChIKey')


ClusterCount <- function(compounds, cpd.identifier, height=0.8) {

  sparse.fp <- sfp(compounds, cpd.identifier)
  dist.matrix <- rpuDist(sparse.fp, method='binary') / dim(sparse.fp)[2]
  candidate.tree <- rpuHclust(dist.matrix)
  length(unique(cutree(candidate.tree, h=height)))
  
}