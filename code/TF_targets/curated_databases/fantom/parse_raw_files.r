rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


gold_standard = read.delim(file = 'data/TF_target_sources/curated_databases/fantom4/edge.GoldStd_TF.tbl.txt', sep ='\t', skip = 24, stringsAsFactors = F)
xref = read.delim(file = 'data/TF_target_sources/curated_databases/fantom4/feature.entrez_gene.tbl.txt', sep ='\t', skip = 29, stringsAsFactors = F)


gold_standard$TF = xref$primary_name[ match(gold_standard$feature1_id, xref$feature_id) ]
gold_standard$target = xref$primary_name[ match(gold_standard$feature2_id, xref$feature_id) ]
gold_standard$pubmed = gold_standard$sym_value


TFs = load_TFs_census()
table(gold_standard$TF %in% TFs)
setdiff(gold_standard$TF, TFs)
table(gold_standard$target %in% TFs)


gold_standard = gold_standard[ gold_standard$TF %in% TFs, ]
gold_standard = gold_standard[ order(gold_standard$target),  ]
gold_standard = gold_standard[ order(gold_standard$TF),  ]
write.table(gold_standard, file = 'data/TF_target_sources/curated_databases/fantom4/edge.GoldStd_TF.tbl.clean.txt', row.names = F, quote = F, sep = '\t')
write.table(unique(gold_standard[, c('TF', 'target') ]), file = 'data/TF_target_sources/curated_databases/fantom4/network.sif', row.names = F, quote = F, sep = '\t')
write.table(unique(gold_standard[, c('TF', 'target', 'pubmed') ]), file = 'data/TF_target_sources/curated_databases/fantom4/network_pubmed.txt', row.names = F, quote = F, sep = '\t')
summarize_network(unique(gold_standard[, c('TF', 'target') ]))
