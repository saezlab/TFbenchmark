rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(viper)



TFs = load_TFs_census()


# Load network_files
network_files = list.files('data/TF_target_sources/', recursive = T, pattern = 'network') %>% grep('consensus', ., value = T, invert = T) %>% grep('sif', ., value=T) %>% grep('old', ., value = T, invert = T) %>% grep('RANDOM_NETWORKS', ., invert = T, value = T)
network_files = setdiff(network_files, 'inferredGTEx/network.sif')
network_files = setdiff(network_files, 'curated_databases/cnio/network.sif')
NETWORKS = list()
for (n in network_files){
  message(n)
  # Load network and fromat regulons
  net = read.delim(paste('data/TF_target_sources/', n, sep = ''), stringsAsFactors = F, sep = '\t')
  net = unique(subset(net, TF %in% TFs))
  NETWORKS[[n]] = net
}


# Define background targets
targets = unlist(sapply(NETWORKS, function(net) net$target))
u_targets = sort(unique(targets))
f_targets = table(targets)
sizes = 1100:3
n_randoms = 1:100

# Build random network
for (si in sizes){
  message('size ', si)
  randoms = list()
  for (ri in n_randoms ){
    randoms[[paste('r', ri, '_', si, sep = '')]] = sample(u_targets, replace = F, prob = f_targets, size = si)
  }
  networks_random = melt(randoms)[, 2:1]
  names(networks_random) = c("TF", "target")
  # write.table(networks_random, file= paste('data/TF_target_sources/RANDOM_NETWORKS_networkPrior/', si, '_network.sif', sep = ''), col.names = T, row.names = F, quote = F, sep = '\t')
  viper_regulon = list()
  for (tf in unique(networks_random$TF) ){
    tf_targets = subset(networks_random, TF == tf)$target
    # Define TF mode
    tf_targets_sign = rep(1, length(tf_targets))
    # Define regulation weight
    tf_targets_weight = rep(1, length(tf_targets))
    # Build VIPER regulon
    viper_regulon[[tf]]$tfmode = tf_targets_sign # NOTE: ttmode has to be the first element of the list. Otherwise msviper will fail.
    names(viper_regulon[[tf]]$tfmode) = tf_targets
    viper_regulon[[tf]]$likelihood = tf_targets_weight
  }
  save(viper_regulon, file = paste('data/TF_target_sources/RANDOM_NETWORKS_networkPrior/', si, '_network_viperRegulon.rdata', sep = '') )
}



