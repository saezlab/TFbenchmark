rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')




data = read.delim(pipe('grep sapiens data/TF_target_sources/curated_databases/oreganno/ORegAnno_Combined_2016.01.19.tsv | grep TRANSCRIPTION | grep hg38'), stringsAsFactors = F, sep = '\t', header = F)
nrow(data)
head(data)
table(data$V13) # Dataset
table(data$V10) # Dataset
data = subset(data, V10 != "miRBase")
data = subset(data, V8 %in% load_TFs_census() )
data = subset(data, V5 != 'N/A' )
nrow(data)
table(data$V13)


table(data$V13)
data$TF = data$V8
data$target = data$V5
data$effect = 0
data$effect[ data$V3 == 'POSITIVE OUTCOME' ] = 1
data$effect[ data$V3 == 'NEGATIVE OUTCOME' ] = - 1
data = data[ order(data$TF) , ]
data$pubmed = data$V12


oreganno = unique(subset(data, V13 %in% c("N/A", 'STAT1 literature-derived sites'))[, c('TF', 'target', 'effect', 'pubmed') ])
NFIRegulomeDB = unique(subset(data, V13 == "NFIRegulomeDB")[, c('TF', 'target', 'effect') ])
pazar = unique(subset(data, V13 == "PAZAR")[, c('TF', 'target', 'effect') ])
summarize_network(oreganno)
summarize_network(NFIRegulomeDB)
summarize_network(pazar)
write.table(unique(oreganno[, c('TF', 'target', 'effect') ]), file = 'data/TF_target_sources/curated_databases/oreganno/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(oreganno, file = 'data/TF_target_sources/curated_databases/oreganno/network_pubmed.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(NFIRegulomeDB, file = 'data/TF_target_sources/curated_databases/NFIRegulomeDB/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(pazar, file = 'data/TF_target_sources/curated_databases/PAZAR/network.sif', sep = '\t', col.names = T, row.names = F, quote = F)
