library(data.table)
library(Matrix)
library(plyr)
library(reshape2)
library(rpudplus)

#' CrossDiverse
#'
#' Find candidate compounds unique to a reference set
#' 
#' @param candidate data frame or data table of candidate compounds with ECFP finger prints and compound identifier columns
#' @param reference data frame or data table of candidate compounds with ECFP finger prints and compound identifier columns
#' @param height numeric value to cut dendrogram defaults to 0.5
#' @export
#' @examples
#' hashes <- round(rnorm(25, 100000000, 1000), 1)
#' 
#' plate.a <- data.table('ECFP_4[1]' = sample(hashes, 100, replace=TRUE),
#'                       'ECFP_4[2]' = sample(hashes, 100, replace=TRUE),
#'                       'ECFP_4[3]' = sample(hashes, 100, replace=TRUE),
#'                       'ECFP_4[4]' = sample(hashes, 100, replace=TRUE),
#'                       'ECFP_4[5]' = sample(hashes, 100, replace=TRUE))
#'                       
#' plate.a$InChIKey <- sapply(1:dim(plate.a)[1], function(x)
#'                            paste(sample(LETTERS, 10, replace=TRUE), collapse=''))
#'                              
#' plate.b <- data.table('ECFP_4[1]' = sample(hashes, 100, replace=TRUE),
#'                       'ECFP_4[2]' = sample(hashes, 100, replace=TRUE),
#'                       'ECFP_4[3]' = sample(hashes, 100, replace=TRUE),
#'                       'ECFP_4[4]' = sample(hashes, 100, replace=TRUE),
#'                       'ECFP_4[5]' = sample(hashes, 100, replace=TRUE))
#'
#' plate.b$InChIKey <- sapply(1:dim(plate.b)[1], function(x)
#'                            paste(sample(LETTERS, 10, replace=TRUE), collapse=''))
#'
#' CrossDiverse(plate.a, plate.b, cpd.identifier='InChIKey')

CrossDiverse <- function (candidate, selection, cpd.identifier, height=0.5){
  
  #row bind the selection and candidate data
  combine.data <- rbind(selection, candidate)
  sparse.fp <- sfp(combine.data, cpd.identifier)
  dist.matrix <- as.matrix(rpuDist(sparse.fp, method='binary') / dim(sparse.fp)[2])

  id.idx <- grep(cpd.identifier, names(combine.data))
  selection.distance <- dist.matrix[which(row.names(dist.matrix) %in% candidate[[id.idx]]),
                                    -(which(row.names(dist.matrix) %in% candidate[[id.idx]]))]
  non.unique <- row.names(selection.distance)[apply(selection.distance, 1, function(x) any(x < height))]

  #Once a compound is not unique to the selection it will never be unique to the
  #selection, there for we can remove the compound from future consideration
  candidate <- candidate[!(candidate[[id.idx]] %in% non.unique), ]
  
  return(candidate)
}