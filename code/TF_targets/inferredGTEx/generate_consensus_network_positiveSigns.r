rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')



# Load networks
aracne_path = 'data/TF_target_sources/inferredGTEx/aracne/'
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
write.table(network[, 1:3], file = 'data/TF_target_sources/inferredGTEx/network_10.sif', sep = '\t', col.names = T, row.names = F, quote = F)


network = unique(df[ df$id_freq > 4 , 1:3])
network$TFtarget = paste(network$TF, network$target)
network = subset(network, ! network$TFtarget %in% names(which(table(network$TFtarget) > 1))  )  # to remove duplicated of opposite sign
write.table(network[, 1:3], file = 'data/TF_target_sources/inferredGTEx/network_5.sif', sep = '\t', col.names = T, row.names = F, quote = F)


network = unique(df[ df$id_freq > 2 , 1:3])
network$TFtarget = paste(network$TF, network$target)
network = subset(network, ! network$TFtarget %in% names(which(table(network$TFtarget) > 1))  )  # to remove duplicated of opposite sign
write.table(network[, 1:3], file = 'data/TF_target_sources/inferredGTEx/network_3.sif', sep = '\t', col.names = T, row.names = F, quote = F)


network = unique(df[ df$id_freq > 1 , 1:3])
network$TFtarget = paste(network$TF, network$target)
network = subset(network, ! network$TFtarget %in% names(which(table(network$TFtarget) > 1)) ) # to remove duplicated of opposite sign
write.table(network[, 1:3], file = 'data/TF_target_sources/inferredGTEx/network_2.sif', sep = '\t', col.names = T, row.names = F, quote = F)



# Build dynamic consensus
TFxTissue = read.delim('/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_info/expression/GTEx_TF_per_tissue.txt', stringsAsFactors = F, header = F)
TFxTissue$V2 = gsub(' ', '_', TFxTissue$V2)
TFxTissue = subset( TFxTissue, V2 %in% df$tissue )
TFxTissue_Freq = table(TFxTissue$V1)
summary(as.integer(TFxTissue_Freq))
TFTGs = list()
for(tf in sort(unique(TFxTissue$V1)) ){
  n_tissues_expressed = TFxTissue_Freq[tf]
  n = 1
  if( n_tissues_expressed %in% 3:8 )
    n = 2
  if( n_tissues_expressed %in% 9:13 )
    n = 3
  if( n_tissues_expressed %in% 14:18 )
    n = 4
  if( n_tissues_expressed %in% 19:23 )
    n = 5
  if( n_tissues_expressed %in% 24:27 )
    n = 6
  tfnet = unique(subset(df, TF == tf & id_freq >= n)[, 1:3])
  tfnet$TFtarget = paste(tfnet$TF, tfnet$target)
  # to remove duplicated of opposite sign
  tfnet = subset(tfnet, ! tfnet$TFtarget %in% names(which(table(tfnet$TFtarget) > 1)) )
  TFTGs[[tf]] = tfnet
  message('TF: ', tf, ' has ', nrow(TFTGs[[tf]]), ' TF-TG with n > ', n)
}
TFxTissue_Freq[ names(which(sapply(TFTGs, nrow) < 10)) ]
names(which(sapply(TFTGs, nrow) > 1000))
network = melt(TFTGs, id.vars = names(TFTGs[[1]]) )
write.table(unique(network[, 1:3]), file = 'data/TF_target_sources/inferredGTEx/network_dynamic.sif', sep = '\t', col.names = T, row.names = F, quote = F)

