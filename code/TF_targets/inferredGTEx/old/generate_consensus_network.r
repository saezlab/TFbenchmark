rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)



# Load networks
aracne_path = 'data/TF_target_sources/inferredGTEx/aracne/'
networks = list()
for (tissue in list.files(path = aracne_path, full.names = T, recursive = T, pattern = 'network')){
  tissue_name = unlist(strsplit(tissue, split = '/'))[6]
  message(tissue_name)
  networks[[tissue_name]] = read.delim(file = tissue, stringsAsFactors = F, header = T)
}
networks = Filter(Negate(function(x){ nrow(x) == 0}), networks)
length(networks)




df = melt(networks, id.vars = names(networks[[1]]))
df$value = 1
m = acast(df, Regulator~Target, fill = 0)
network = melt(m)
names(network)[1:2] = c('TF', 'target')
network = network[ order(network$target) , ]
network = network[ order(network$TF) , ]
write.table(network[ network$value > 9 , 1:2], file = 'data/TF_target_sources/inferredGTEx/network_10.sif', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(network[ network$value > 4 , 1:2], file = 'data/TF_target_sources/inferredGTEx/network_5.sif', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(network[ network$value > 1 , 1:2], file = 'data/TF_target_sources/inferredGTEx/network_2.sif', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(network[ network$value > 2 , 1:2], file = 'data/TF_target_sources/inferredGTEx/network_3.sif', sep = '\t', col.names = T, row.names = F, quote = F)
