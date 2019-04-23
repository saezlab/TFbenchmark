rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


data = read.csv('data/TF_target_sources/curated_databases/cnio/TF-TG_compiled_including_TM-set_190118.csv')
dim(data)
data = subset(data, X.FNL..Confidence == 'High')
data = subset(data, TF %in% load_TFs_census() )
dim(data)
head(data)
network = unique(data[, 2:3])
names(network) = c('TF', 'target')
network = network[ order(network$target), ]
network = network[ order(network$TF), ]
summarize_network(network)
network_similarity(network)
write.table(network, file = 'data/TF_target_sources/curated_databases/cnio/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)
