rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/functions.r')





# Define networks files to merge
# Load curated networks but cnio
networks2merge = c()
networks = list.files('data/TF_target_sources/curated_databases', pattern = 'network', recursive = T, full.names = T) %>% 
  grep('sif', ., value=T) %>% 
  grep('cnio', ., invert = T, value = T)
# remove unsigned networks (when they have the signed version)
networks = grep('trrust/', networks, invert = T, value = T)  %>% 
  grep('trrd_via_tfact/', ., invert = T, value = T)  %>% 
  grep('TFe/', ., invert = T, value = T)  %>% 
  grep('tfact/', ., invert = T, value = T)  %>% 
  grep('oreganno/', ., invert = T, value = T)  
networks2merge = append(networks2merge, values = networks)
# Select network/parameters
# inferred = 3 tissues
networks = list.files('data/TF_target_sources/inferred/GTEx/pantissue',  pattern = 'network_i3.sif', recursive = T, full.names = T)
networks2merge = append(networks2merge, values = networks)
# TFBS_scanning size = 500
networks = list.files('data/TF_target_sources/TFBS_scanning', pattern = 'n500.sif', recursive = T, full.names = T)
networks2merge = append(networks2merge, values = networks)
# ChIP_Seq size = 500
networks = list.files('data/TF_target_sources/ChIP_Seq', pattern = 'n500.sif', recursive = T, full.names = T)
networks2merge = append(networks2merge, values = networks)
networks2merge




# Merge networks into a data frame
df = load_and_merge_regulons_asdf(networks2merge)
write.table(df, file = 'data/TF_target_sources/merged_interactions.txt', quote = F, row.names = F, col.names = T, sep = '\t')
nrow(unique(df[, 1:2]))




# score interactions
df$value = 1
n_review = acast(subset(df, database %in% c('reviews', 'TFe_signed')), TF~target)
n_curated = acast(subset(df, evidence == 'curated_databases'), TF~target)
n_curated_signed = acast(subset(df, evidence == 'curated_databases' & effect != 0 ), TF~target)
n_ChIP_Seq = acast(subset(df, evidence == 'ChIP_Seq'), TF~target)
n_TFBS_scanning = acast(subset(df, evidence == 'TFBS_scanning'), TF~target)
n_inferred = acast(subset(df, evidence == 'inferred'), TF~target)
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
shared_TFs = Reduce(intersect, sapply(list(n_ChIP_Seq, n_curated, n_inferred, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_ChIP_Seq, n_curated, n_inferred, n_TFBS_scanning), colnames) )
INTERACTIONS$A = rbind(INTERACTIONS$A ,
                       subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_ChIP_Seq[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0 & n_inferred[ shared_TFs, shared_targets ] > 0) , value)[, 1:2] )
nrow(unique(INTERACTIONS$A))
length(unique(INTERACTIONS$A$Var1))



## B
# curated & ChIP_Seq
shared_TFs = Reduce(intersect, sapply(list(n_ChIP_Seq, n_curated), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_ChIP_Seq, n_curated), colnames) )
INTERACTIONS$B = subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_ChIP_Seq[ shared_TFs, shared_targets ] > 0) , value)[, 1:2]
# curated & predictions
shared_TFs = Reduce(intersect, sapply(list(n_curated, n_inferred, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_curated, n_inferred, n_TFBS_scanning), colnames) )
INTERACTIONS$B = rbind(INTERACTIONS$B,
                       subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0 & n_inferred[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
# ChIP_Seq & predictions
shared_TFs = Reduce(intersect, sapply(list(n_ChIP_Seq, n_inferred, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_ChIP_Seq, n_inferred, n_TFBS_scanning), colnames) )
INTERACTIONS$B = rbind(INTERACTIONS$B,
                       subset(melt(n_ChIP_Seq[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0 & n_inferred[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
## C
# curated & TFBS
shared_TFs = Reduce(intersect, sapply(list(n_curated, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_curated, n_TFBS_scanning), colnames) )
INTERACTIONS$C = subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0) , value)[, 1:2]
# curated & inferred
# shared_TFs = Reduce(intersect, sapply(list(n_curated, n_inferred), rownames) )
# shared_targets = Reduce(intersect, sapply(list(n_curated, n_inferred), colnames) )
# INTERACTIONS$C = rbind(INTERACTIONS$C,
                       # subset(melt(n_curated[ shared_TFs, shared_targets ] > 0 & n_inferred[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
# ChIP_Seq & any prediction
shared_TFs = Reduce(intersect, sapply(list(n_ChIP_Seq, n_inferred), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_ChIP_Seq, n_inferred), colnames) )
INTERACTIONS$C = rbind(INTERACTIONS$C,
                       subset(melt(n_ChIP_Seq[ shared_TFs, shared_targets ] > 0 & n_inferred[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
shared_TFs = Reduce(intersect, sapply(list(n_ChIP_Seq, n_TFBS_scanning), rownames) )
shared_targets = Reduce(intersect, sapply(list(n_ChIP_Seq, n_TFBS_scanning), colnames) )
INTERACTIONS$C = rbind(INTERACTIONS$C,
                       subset(melt(n_ChIP_Seq[ shared_TFs, shared_targets ] > 0 & n_TFBS_scanning[ shared_TFs, shared_targets ] > 0) , value)[, 1:2])
## D
# curated > 0
INTERACTIONS$D = subset(melt(n_curated), value > 0)[, 1:2]
# ChIP_Seq > 0
INTERACTIONS$D = rbind(INTERACTIONS$D,
                       subset(melt(n_ChIP_Seq), value > 0)[, 1:2])
## E
# TFBS > 0
INTERACTIONS$E = subset(melt(n_TFBS_scanning), value > 0)[, 1:2]
# inferred > 0
INTERACTIONS$E = rbind(INTERACTIONS$E,
                       subset(melt(n_inferred), value > 0)[, 1:2])


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
