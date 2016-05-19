library(data.table)
library(Matrix)
library(optparse)
library(plyr)
library(reshape2)
library(rpudplus)
library(snow)

library(divplate)

#' PlateWiseDiversity
#' 
#' Main commandline program for plate based diversity selection. The algorithim
#' selectes the most diverse candidate plate relative to a reference plate. Compounds
#' on the candidate plate that are not unique to the reference plate are removed from 
#' future consideration. This approach yields a local maximum diversity for a set of
#' plates. 

#main
option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA,
              type='character', help="pointer to input library file"),
  make_option(c("-c", "--id"), action="store", default='InChIKey',
              type='character', help="compound identifier string"),
  make_option(c("-n", "--number"), action="store", default=NA, 
              type='numeric', help="number of plates to select"),
  make_option(c("-a", "--intra.height"), action="store", default=0.8, 
              type='numeric', help="intra diversity cut height"), 
  make_option(c("-e", "--inter.height"), action="store", default=0.5, 
              type='numeric', help="inter diversity cut height"),  
  make_option(c("-g", "--gpus"), action="store", default=5, 
              type='numeric', help="number of parallel jobs") 
)

opt = parse_args(OptionParser(option_list=option_list))

file.loc <- opt$input
cpd.identifier <- opt$id
intra.height <- opt$intra.height
inter.height <- opt$inter.height
n <- opt$number
n.cores <- opt$gpus

#set up snow
clus <- makeCluster(n.cores)
#export variables
clusterExport(clus, c('cpd.identifier',
                      'intra.height', 
                      'inter.height'), 
              envir=environment())
#export functions
clusterExport(clus, 'sfp')
clusterExport(clus, 'ClusterCount')
clusterExport(clus, 'CrossDiverse')
clusterExport(clus, 'melt')
clusterExport(clus, 'mapvalues')
clusterExport(clus, 'sparseMatrix')
clusterExport(clus, 'rpuDist')
clusterExport(clus, 'rpuHclust')
clusterCall(clus, function() library(data.table))

#read in compound library
plate.set <- fread(file.loc)

#make a full copy of the compound library
plate.repository <- plate.set

#find most internally diverse plate
plates <- split(plate.set, plate.set$PLATE_MAP_NAME)
n.initial <- parSapply(clus, plates, 
                       function(plate) ClusterCount(plate, 
                                                    cpd.identifier, 
                                                    intra.height))

selection <- plates[which(n.initial == max(n.initial))[1]]

#create a tracker for plate diversity
tracker <- data.frame(plate=names(plates))
tracker$plate <- as.character(tracker$plate)
plate.size <- sapply(plates, function(plate) dim(plate)[1])
tracker$plate.size <- plate.size
tracker$n.initial <- n.initial

#iterate through plates
candidates <- plates[names(plates) != names(selection)]
i <- 4
while (1) {

  #compare candidates to the selected plate
  clusterExport(clus, c('selection'), 
                envir=environment())
  candidates <- parLapply(clus, candidates,
                          function(candidate) 
                          CrossDiverse(candidate, 
                                       selection[[1]], 
                                       cpd.identifier, 
                                       height=inter.height))
  
  #find internal diversity of remaining candidates. Uses inter.height
  diverse.count <- parSapply(clus, candidates,
                             function(plate)
                             ClusterCount(plate,
                                          cpd.identifier,
                                          inter.height))

  #append relative diversity count to tracker
  tracker$latest <- sapply(tracker$plate, 
                          function(plate) 
                            diverse.count[which(plate == names(diverse.count))])
  
  #name the selecte plate
  names(tracker)[i] <- names(selection)
  
  #select the new plate
  selection <- candidates[which(diverse.count == max(diverse.count))[1]]
  candidates <- candidates[names(candidates) != names(selection)]
  
  #increase the counter by 1
  i <- i + 1
  
  #break loop when the desired number of plates has been selected
  if (i >= (n + 4)) {break}
}

#tidy  up tracker
plate.size <- array(0, dim(tracker)[2] - 3)
internal.diversity <- array(0, dim(tracker)[2] - 3)
n.unique <- array(0, dim(tracker)[2] - 3)
for (i in 4:dim(tracker)[2]) {
  plate.index <- which(names(tracker)[i] == tracker$plate)
  plate.size[i - 3] <- tracker[plate.index, 2]
  internal.diversity[i - 3] <- tracker[plate.index, 3]
  n.unique[i - 3] <- tracker[plate.index, (i - 1)]
}
selected <- data.frame(plate=names(tracker)[4:dim(tracker)[2]])
selected$plate.size <- plate.size
selected$internal.diversity <- internal.diversity
selected$n.unique <- n.unique

#write selection
selected <- apply(selected, 2, as.character)
write.csv(selected, 'plate.wise.diversity.selection.csv')