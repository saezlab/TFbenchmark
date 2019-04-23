rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')



omnipath_TFTG_database = read.csv(file = 'data/TF_target_sources/omnipath_scores/database_PANCANCER.csv', stringsAsFactors = F)
omnipath_TFTG_database = subset(omnipath_TFTG_database, ! is.na(omnipath_TFTG_database$score))

# Network per score
for (sc in unique(omnipath_TFTG_database$score) ){
  message(sc)
  df = subset(omnipath_TFTG_database, score %in% LETTERS[1:which(LETTERS == sc)]  )[, 1:3]
  # df = subset(omnipath_TFTG_database, score %in% sc  )[, 1:3]
  write.table(df, file = paste('data/TF_target_sources/omnipath_scores/', sc, '/network_PANCANCER.sif', sep = ''), col.names = T, row.names = F, sep = '\t', quote = F)
}

# Network from top scored (n>9)
TFs = names(which(table(omnipath_TFTG_database$TF) > 3))
top_interactions = list()
for (tf in TFs){
  tf_int = subset(omnipath_TFTG_database, TF == tf)
  tf_int = subset(tf_int, target != tf)
  counts = sapply(LETTERS[1:5], function(le) nrow(subset(tf_int, score %in%  LETTERS[1:which(LETTERS == le)] )) )
  sc = names(which(counts >= 10)[1])
  if( is.na(sc) )
    sc = names(which.max(counts))[1]
  df = subset(tf_int, score %in% LETTERS[1:which(LETTERS == sc)]  )[, 1:3]
  message(tf, ' ', sc, ' ', nrow(df) )
  top_interactions[[paste(tf, '_', sc, sep = '')]] = df
}
regulons = melt(top_interactions, id.vars = names(top_interactions[[1]]) )
regulons$TF = regulons$L1
write.table(regulons[, 1:3], file = 'data/TF_target_sources/omnipath_scores/TOP/network_PANCANCER.sif', sep = '\t', quote = F, row.names = F)
