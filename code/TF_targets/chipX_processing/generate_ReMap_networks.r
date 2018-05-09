rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')



TFcensus = load_TFs_census()
remap_network = read.delim('data/TF_target_sources/chipSeq/ReMap/gene_tf_pairs_genesymbol.txt', stringsAsFactors = F, header = F)
names(remap_network) = c('TF', 'target', 'ensembl_target', 'score')
remap_network = subset(remap_network, TF %in% TFcensus )
remap_network = subset(remap_network, score > 100 )
tfs = unique(remap_network$TF) 
length(tfs)
summarize_network(remap_network)

for (i in c(100, 200, 500, 1000) ){
  message(i)
  topEdges = lapply(tfs, function(tf){
    score_at_threshold = min(head(subset(remap_network, TF == tf ), i)$score)
    tfnetwork = subset(remap_network, score > score_at_threshold & TF == tf) [, c('TF', 'target') ]
    return(tfnetwork)
  })
  network = melt(topEdges, id.vars = names(topEdges[[1]]))[, -3]
  summarize_network(network)
  filename = paste('data/TF_target_sources/chipSeq/ReMap/top_', i,'_network.sif', sep = '')
  write.table(network, file = filename, sep = '\t', col.names = T, row.names = F, quote = F)
}




