rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(magrittr)




data = read.delim('data/TF_target_sources/curated_databases/IntAct/homo_sapiens_protein2gene_20171031.txt', stringsAsFactors = F, sep = '\t', header = T)
idx.reorder = data$Type.s..interactor.A == 'psi-mi:MI:0250(gene)'
network1 = unname(data[ ! idx.reorder, c(1:2,9)])
network2 = unname(data[   idx.reorder , c(2:1,9)])
names(network1)  = c('TF', 'target', 'pubmed')
names(network2)  = c('TF', 'target', 'pubmed')
network = rbind(network1, network2)
dim(network)

network$TF = sapply(strsplit(network$TF, split = ':'), tail, 1)
network$TF = sapply(strsplit(network$TF, split = '-'), head, 1)
network$target = sapply(strsplit(network$target, split = ':'), tail, 1)
network$pubmed = lapply(strsplit(network$pubmed, split = '\\|'), grep, pattern = "pubmed", value = T) %>% sapply(., paste, collapse = ",") %>% gsub('pubmed:', '', .)


load('/Volumes/GoogleDrive/My Drive/databases/ensemblgenes_annot.Rdata')
network$TF = ensemblgenes_annot$hgnc_symbol[ match(network$TF, ensemblgenes_annot$uniprot_swissprot) ]
network$target = ensemblgenes_annot$hgnc_symbol[ match(network$target, ensemblgenes_annot$ensembl_gene_id) ]
network = network[ ! (is.na(network$TF) | is.na(network$target)), ]
dim(network)

tfs = load_TFs_census()
network = unique(subset(network, TF %in% tfs))
write.table(network, file = 'data/TF_target_sources/curated_databases/IntAct/network_pubmed.txt', quote = F, row.names = F, sep = '\t')
network = unique(network[, 1:2])
dim(network)

network = network[order(network$target), ]
network = network[order(network$TF), ]

write.table(network, file = 'data/TF_target_sources/curated_databases/IntAct/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)
