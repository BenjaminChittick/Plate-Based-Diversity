
library(data.table)
library(Matrix)
library(plyr)
library(reshape2)

#' sfp
#' 
#' Converts ECFP hash codes into a sparse matrix by arbitrarily assigning them as indices.
#' 
#' @param data data.table containing compound identifier and fingerprint hash
#' @param cpd.identifier  character string identifying compounds: Broad ID, InChIKey, etc...
#' @export
#' @examples
#' hashes <- round(rnorm(120, 100000000, 1000), 1)
#' compounds <- data.table('ECFP_4[1]' = sample(hashes, 100, replace=TRUE),
#'                         'ECFP_4[2]' = sample(hashes, 100, replace=TRUE),
#'                         'ECFP_4[3]' = sample(hashes, 100, replace=TRUE),
#'                         'ECFP_4[4]' = sample(hashes, 100, replace=TRUE),
#'                         'ECFP_4[5]' = sample(hashes, 100, replace=TRUE))
#' compounds$InChIKey <- sapply(1:dim(compounds)[1], function(x)
#'                              paste(sample(LETTERS, 10, replace=TRUE), collapse=''))
#' sparseFP <- sfp(compounds, cpd.identifier='InChIKey')
#' sparseFP[1:10, 1:15]

sfp <- function(data, cpd.identifier){

  fp.idx <- grep('CFP_', names(data))
  id.idx <- grep(cpd.identifier, names(data))
 
  data <- data[!duplicated(data[, id.idx, with=FALSE]), ]
  
  data <- melt(data[, c(id.idx, fp.idx), with=FALSE], 
               id.var=cpd.identifier, na.rm=TRUE)
  
  data <- setorderv(data, cpd.identifier)
  setnames(data, 'value', 'col')
  data$col <- as.factor(data$col)
  data$col <- mapvalues(data$col, from=unique(data$col),
                        to=1:length(unique(data$col)))
  
  data$row <- mapvalues(data[[1]], 
                        from=unique(data[[1]]),
                        to=1:length(unique(data[[1]])))
  
  sparsePrint <- sparseMatrix(as.numeric(as.character(data$row)), 
                              as.numeric(as.character(data$col)))
  row.names(sparsePrint) <- unique(data[[1]])
  
  return(sparsePrint)
}