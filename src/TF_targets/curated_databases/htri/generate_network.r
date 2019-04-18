rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


data = read.delim('data/TF_target_sources/curated_databases/HTRIdb/tf-target_network_052016_literaturecurated.txt', stringsAsFactors = F)
head(data)
nrow(data)
table(data$TECHNIQUE)


network = unique(data[, c(3, 5, 7)])
names(network) = c('TF', 'target', 'pubmed')
network = network[ order(network$target), ]
network = network[ order(network$TF), ]
network = subset(network, TF %in% load_TFs_census() )
summarize_network(network)
network_similarity(network)
head(network)
write.table(unique(network[, 1:2]), file = 'data/TF_target_sources/curated_databases/HTRIdb/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(unique(network), file = 'data/TF_target_sources/curated_databases/HTRIdb/network_pubmed.txt', sep = '\t', col.names = T, row.names = F, quote = F)
