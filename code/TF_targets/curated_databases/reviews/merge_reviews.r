rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


files = list.files('data/TF_target_sources/curated_databases/reviews/', pattern = 'sif', full.names = T, recursive = T)
files = grep('network', files, invert = T, value = T)

data = lapply(files, read.table, stringsAsFactors = F, sep = '\t') %>% lapply(., function(x) x[, 1:2] )
names(data) = strsplit(files, '/') %>%  lapply(., tail ,2) %>% sapply(., head, 1)  %>% strsplit(., '_') %>%  sapply(., tail ,1)
network = melt(data)
names(network) = c("TF", "target", "pubmed")
write.table(network, file = "data/TF_target_sources/curated_databases/reviews/network_pubmed.txt", sep = "\t", quote = F, row.names = F)

