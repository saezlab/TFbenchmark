rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


data = read.delim('data/TF_target_sources/curated_databases/trrust/trrust_rawdata.human.v2.tsv', stringsAsFactors = F, header = F)
data$TFtarget = paste(data$V1, data$V2)
head(data)
nrow(data)
# table(data$V3)
data$evidence_count = sapply(strsplit(data$V4, split = ';'), function(x) length(unlist(x)) )
data = data[ order(data$V3), ] # order according the effect sign activation > repression > unknown
data = data[ order(data$evidence_count, decreasing = T), ] # order according the evidence_count
data = rbind(subset(data, V3 != 'Unknown'), subset(data, V3 == 'Unknown')) # prioritize not unknown
data = data[ order(data$V2), ] # order according the target
data = data[ order(data$V1), ] # order according the TF
data = data[ ! duplicated(data[, 1:2]) , ]



network = unique(data[, c(1:3)])
names(network) = c('TF', 'target', 'effect')
table(network$effect)
network$effect = c('1', '0', '-1')[ match(network$effect, c('Activation', 'Unknown', 'Repression') ) ]
network = network[ order(network$target), ]
network = network[ order(network$TF), ]
network = subset(network, TF %in% load_TFs_census() )
write.table(network[, 1:2], file = 'data/TF_target_sources/curated_databases/trrust/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(network, file = 'data/TF_target_sources/curated_databases/trrust_signed/network_signed.sif', sep = '\t', col.names = T, row.names = F, quote = F)
summarize_network(network)
network_similarity(network)
