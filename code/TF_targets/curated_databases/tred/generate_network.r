rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


data = read.csv('data/TF_target_sources/curated_databases/tred_via_RegNetwork/export_Fri_Sep_22_11_00_16_UTC_2017.csv', stringsAsFactors = F)
head(data)
nrow(data)
table(data$confidence)
table(data$evidence)
table(data$database)


network = unique(data[, c(1,3)])
names(network) = c('TF', 'target')
network = network[ order(network$target), ]
network = network[ order(network$TF), ]
network = subset(network, TF %in% load_TFs_census() )
summarize_network(network)
network_similarity(network)
write.table(network, file = 'data/TF_target_sources/curated_databases/tred_via_RegNetwork/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)
