rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/data.r')


TF_census = load_TFs_census()
  
# Load networks
aracne_path = 'data/TF_target_sources/inferred/GTEx/tissue_specific_noCOMBAT/'
networks = list()
for (tissue in list.files(path = aracne_path, full.names = T, recursive = T, pattern = 'viperRegulon')){
  tissue_name = unlist(strsplit(tissue, split = '/'))[6]
  message(tissue_name)
  regulon = get(load(tissue))
  names(regulon) = sapply(strsplit(names(regulon), split = ' '), head, 1)
  regulon_sign = lapply(regulon, function(x) data.frame(target = names(x$tfmode), x$tfmode, stringsAsFactors = F) )
  names(regulon_sign) = names(regulon)
  networks[[tissue_name]] = regulon_sign
}
length(networks)


# MERGE data and format the effect sign
df = melt(networks, id.vars = names(networks[[1]][[1]]) )
names(df) = c('target', 'effect', 'TF', 'tissue')
df = df[, c('TF', 'target', 'effect', 'tissue') ]
df = subset(df, TF %in% TF_census)
df$effect[ df$effect > 0 ] = 1
df$effect[ df$effect < 0 ] = -1



# extract frequency of the signed interaction and remove singletons
df$id = paste(df$TF, df$target, df$effect)
freq = table(df$id)
df$id_freq = freq[ df$id ]
df = subset(df, id_freq > 1)



# BUILD consensus
network = unique(df[ df$id_freq > 9 , 1:3])
network$TFtarget = paste(network$TF, network$target)
network = subset(network, ! network$TFtarget %in% names(which(table(network$TFtarget) > 1))  )  # to remove duplicated of opposite sign
write.table(network[, 1:3], file = 'data/TF_target_sources/inferred/GTEx/pantissue/network_noCOMBAT_i10.sif', sep = '\t', col.names = T, row.names = F, quote = F)


network = unique(df[ df$id_freq > 4 , 1:3])
network$TFtarget = paste(network$TF, network$target)
network = subset(network, ! network$TFtarget %in% names(which(table(network$TFtarget) > 1))  )  # to remove duplicated of opposite sign
write.table(network[, 1:3], file = 'data/TF_target_sources/inferred/GTEx/pantissue/network_noCOMBAT_i5.sif', sep = '\t', col.names = T, row.names = F, quote = F)


network = unique(df[ df$id_freq > 2 , 1:3])
network$TFtarget = paste(network$TF, network$target)
network = subset(network, ! network$TFtarget %in% names(which(table(network$TFtarget) > 1))  )  # to remove duplicated of opposite sign
write.table(network[, 1:3], file = 'data/TF_target_sources/inferred/GTEx/pantissue/network_noCOMBAT_i3.sif', sep = '\t', col.names = T, row.names = F, quote = F)


network = unique(df[ df$id_freq > 1 , 1:3])
network$TFtarget = paste(network$TF, network$target)
network = subset(network, ! network$TFtarget %in% names(which(table(network$TFtarget) > 1)) ) # to remove duplicated of opposite sign
write.table(network[, 1:3], file = 'data/TF_target_sources/inferred/GTEx/pantissue/network_noCOMBAT_i2.sif', sep = '\t', col.names = T, row.names = F, quote = F)
