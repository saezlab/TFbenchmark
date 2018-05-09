rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')




# Load networks but randoms, cnio and consensus
networks = list.files('data/TF_target_sources/', recursive = T, pattern = 'network') %>% grep('sif', ., value=T) %>% grep('consensus', ., value = T, invert = T)  %>% 
  grep('omnipath_scores', ., value = T, invert = T)  %>% 
  grep('RANDOM_NETWORKS', ., invert = T, value = T) %>% grep('cnio', ., invert = T, value = T)
# also remove unsigned networks (when they have the signed version)
networks = grep('trrust/', networks, invert = T, value = T)  %>% grep('trrd_via_tfact/', ., invert = T, value = T)  %>% grep('TFe/', ., invert = T, value = T)  %>% grep('tfact/', ., invert = T, value = T)  %>% grep('oreganno/', ., invert = T, value = T)  
# Select network/parameters
# inferredGTEx = 3 tissues
networks = grep('inferredGTEx/network_10', networks, invert = T, value = T) %>% grep('inferredGTEx/network_2', ., invert = T, value = T) %>% grep('inferredGTEx/network_5', ., invert = T, value = T)  %>% grep('inferredGTEx/network_dynamic', ., invert = T, value = T)
# TFBS_scanning size = 500
networks = grep('_con.sif', networks, invert = T, value = T) %>% grep('_opRegion.sif', ., invert = T, value = T)  %>% grep('_n100.sif', ., invert = T, value = T)  %>% grep('_n1000.sif', ., invert = T, value = T)  %>% grep('_n200.sif', ., invert = T, value = T)
# chipSeq size = 500
networks = grep('top_100_network.sif', networks, invert = T, value = T) %>% grep('top_200_network.sif', ., invert = T, value = T)  %>% grep('top_1000_network.sif', ., invert = T, value = T)


regulon = list()
for (n in networks){
  message(n)
  # Load network and fromat regulons
  net = read.delim(paste('data/TF_target_sources/', n, sep = ''), stringsAsFactors = F, sep = '\t')
  evidence = unlist(strsplit(n, '/'))[1]
  database = unlist(strsplit(n, '/'))[2]
  if( database == "scanning_output")
    database = unlist(strsplit(n, '/'))[3]
  if( is.null(regulon[[evidence]]) )
    regulon[[evidence]] = list()
  if( is.null(net$effect) )
    net$effect = 0
  for( tf in unique(net$TF) ){
    regulon[[evidence]][[database]][[tf]] = unique(rbind(regulon[[evidence]][[database]][[tf]],
                                                         subset(net, TF == tf)[, c('target', 'effect') ] ))
  }
}
# build merged data.frame ~ database
df = unique(melt(regulon, id.vars = names(regulon[[1]][[1]][[1]]) ))
df = df[, c(3,1,2,5,4)]
names(df) = c('TF', 'target', 'effect', 'evidence', 'database')
write.table(df, file = 'data/TF_target_sources/merged_interactions.txt', quote = F, row.names = F, col.names = T, sep = '\t')
nrow(unique(df[, 1:2]))


# score interactions
df$value = 1
n_review = acast(subset(df, database %in% c('reviews', 'TFe_signed')), TF~target)
n_curated = acast(subset(df, evidence == 'curated_databases'), TF~target)
n_curated_signed = acast(subset(df, evidence == 'curated_databases' & effect != 0 ), TF~target)
n_chipSeq = acast(subset(df, evidence == 'chipSeq'), TF~target)
n_TFBS_scanning = acast(subset(df, evidence == 'TFBS_scanning'), TF~target)
n_inferredGTEx = acast(subset(df, evidence == 'inferredGTEx'), TF~target)
n_noncurated = acast(subset(df, evidence != 'curated_databases'), TF~target)


INTERACTIONS = list()
## A
# > 2 curated databases
INTERACTIONS$A = rbind(INTERACTIONS$A ,
                       subset(melt(n_curated), value > 1)[, 1:2])
# in curated reviews
INTERACTIONS$A = rbind(INTERACTIONS$A,
                       subset(melt(n_review > 0) , value)[, 1:2] )
# # 2 curated databases & any other
# shared_TFs = Reduce(intersect, sapply(list(n_curated, n_noncurated), rownames) )
# shared_targets = Reduce(intersect, sapply(list(n_curated, n_noncurated), colnames) )
# INTERACTIONS$A = rbind(INTERACTIONS$A , 
#                        subset(melt(n_curated[ shared_TFs, shared_targets ] > 1 & n_noncurated[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
# # in curated signed
# INTERACTIONS$A = rbind(INTERACTIONS$A ,
#                        subset(melt(n_curated_signed > 0) , value)[, 1:2] )
# 1 curated_signed & any other
shared_TFs = Reduce(intersect, sapply(list(n_curated_signed, n_noncurated), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_curated_signed, n_noncurated), colnames) )
INTERACTIONS$A = rbind(INTERACTIONS$A ,
                       subset(melt(n_curated_signed[ shared_TFs, shared_targets ] > 0 & n_noncurated[ shared_TFs, shared_targets ] > 0) , value)[, 1:2] )
# 4 evidences
shared_TFs = Reduce(intersect, sapply(list(n_chipSeq, n_curated, n_inferredGTEx, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_chipSeq, n_curated, n_inferredGTEx, n_TFBS_scanning), colnames) )
INTERACTIONS$A = rbind(INTERACTIONS$A ,
                       subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_chipSeq[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0 & n_inferredGTEx[ shared_TFs, shared_targets ] > 0) , value)[, 1:2] )
nrow(unique(INTERACTIONS$A))
length(unique(INTERACTIONS$A$Var1))



## B
# curated & chipSeq
shared_TFs = Reduce(intersect, sapply(list(n_chipSeq, n_curated), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_chipSeq, n_curated), colnames) )
INTERACTIONS$B = subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_chipSeq[ shared_TFs, shared_targets ] > 0) , value)[, 1:2]
# curated & predictions
shared_TFs = Reduce(intersect, sapply(list(n_curated, n_inferredGTEx, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_curated, n_inferredGTEx, n_TFBS_scanning), colnames) )
INTERACTIONS$B = rbind(INTERACTIONS$B,
                       subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0 & n_inferredGTEx[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
# chipSeq & predictions
shared_TFs = Reduce(intersect, sapply(list(n_chipSeq, n_inferredGTEx, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_chipSeq, n_inferredGTEx, n_TFBS_scanning), colnames) )
INTERACTIONS$B = rbind(INTERACTIONS$B,
                       subset(melt(n_chipSeq[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0 & n_inferredGTEx[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
## C
# curated & TFBS
shared_TFs = Reduce(intersect, sapply(list(n_curated, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_curated, n_TFBS_scanning), colnames) )
INTERACTIONS$C = subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0) , value)[, 1:2]
# curated & inferred
# shared_TFs = Reduce(intersect, sapply(list(n_curated, n_inferredGTEx), rownames) )
# shared_targets = Reduce(intersect, sapply(list(n_curated, n_inferredGTEx), colnames) )
# INTERACTIONS$C = rbind(INTERACTIONS$C,
                       # subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_inferredGTEx[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
# chipSeq & any prediction
shared_TFs = Reduce(intersect, sapply(list(n_chipSeq, n_inferredGTEx), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_chipSeq, n_inferredGTEx), colnames) )
INTERACTIONS$C = rbind(INTERACTIONS$C,
                       subset(melt(n_chipSeq[ shared_TFs, shared_targets ] > 0 & n_inferredGTEx[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
shared_TFs = Reduce(intersect, sapply(list(n_chipSeq, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_chipSeq, n_TFBS_scanning), colnames) )
INTERACTIONS$C = rbind(INTERACTIONS$C,
                       subset(melt(n_chipSeq[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
## D
# curated > 0
INTERACTIONS$D = subset(melt(n_curated), value > 0)[, 1:2]
# chipSeq > 0
INTERACTIONS$D = rbind(INTERACTIONS$D,
                       subset(melt(n_chipSeq), value > 0)[, 1:2])
## E
# TFBS > 0
INTERACTIONS$E = subset(melt(n_TFBS_scanning), value > 0)[, 1:2]
# inferredGTEx > 0
INTERACTIONS$E = rbind(INTERACTIONS$E,
                       subset(melt(n_inferredGTEx), value > 0)[, 1:2])


INTERACTIONS$A = unique(INTERACTIONS$A)
INTERACTIONS$B = unique(INTERACTIONS$B)
INTERACTIONS$C = unique(INTERACTIONS$C)
INTERACTIONS$D = unique(INTERACTIONS$D)
SCORES = melt(INTERACTIONS)
names(SCORES) = c('TF', 'target', 'score')
SCORES$TFtarget = paste(SCORES$TF, SCORES$target)
table(SCORES$score)
SCORES = subset(SCORES, ! duplicated(TFtarget) ) 
table(SCORES$score)
SCORES = SCORES[ order(SCORES$TFtarget), ]
write.csv(SCORES, file = 'data/TF_target_sources/omnipath_scores/TF_target_scores.csv', row.names = F)
