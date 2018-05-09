rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


data = read.delim('data/TF_target_sources/curated_databases/tfact/Catalogues.txt', stringsAsFactors = F)
head(data)
nrow(data)
table(data$SPECIES)
data = data[ unique(c(grep('homo', data$SPECIES, ignore.case = T), grep('human', data$SPECIES, ignore.case = T))) , ]
nrow(data)
table(data$REF)





tfact = data[ grep('Pubmed', data$REF, ignore.case = T) , ]
nrow(tfact)
network = unique(tfact[, c(2:1, 3)])
names(network) = c('TF', 'target', 'effect')
network$effect = c('1', '-1')[ match(network$effect, c('UP', 'DOWN') ) ]
network = network[ order(network$target), ]
network = network[ order(network$TF), ]
network = subset(network, TF %in% load_TFs_census() )
summarize_network(network)
network_similarity(network)
write.table(network, file = 'data/TF_target_sources/curated_databases/tfact/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)



TRRD = data[ grep('TRRD', data$REF, ignore.case = T) , ]
nrow(TRRD)
network = unique(TRRD[, c(2:1, 3)])
names(network) = c('TF', 'target', 'effect')
network$effect = c('1', '-1')[ match(network$effect, c('UP', 'DOWN') ) ]
network = network[ order(network$target), ]
network = network[ order(network$TF), ]
network = subset(network, TF %in% load_TFs_census() )
summarize_network(network)
network_similarity(network)
write.table(network, file = 'data/TF_target_sources/curated_databases/trrd_via_tfact/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)
