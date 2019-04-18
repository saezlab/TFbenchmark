rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


load('data/TF_info/TFrole_genesets.rdata')
TFs = sort(load_TFs_census())
table_S1 = data.frame(TF = TFs, stringsAsFactors = F)
rownames(table_S1) = TFs
dim(table_S1)

for(prop in  names(TFrole_genesets) ){
  table_S1[,prop] = '-'
  for( suprop in names(TFrole_genesets[[prop]]) )
    table_S1[ rownames(table_S1) %in% TFrole_genesets[[prop]][[suprop]], prop ] = suprop
}
head(table_S1)
apply(table_S1, 2, table)
dim(table_S1)
write.csv(table_S1, row.names = F, quote = F, file = 'paper/GarciaAlonso_tableS1.csv')