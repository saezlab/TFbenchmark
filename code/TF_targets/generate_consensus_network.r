rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')



fix_sign_conflict = function (ddf){
  # Identify TF-target interactions whith effect sign conflict (i.e. 1 and -1)
  x = unique(subset(ddf, effect != 0 & evidence == 'curated_databases')[, c('TF', 'target', 'effect', 'TFtarget')])
  TFtargte_with_effect_confict = unique(x$TFtarget[ duplicated(x$TFtarget) ])
  message('NOTE: ', length(TFtargte_with_effect_confict), ' whith effect sign conflict.\nPrioritizing those with matching inferredGTEx effect sign')
  print(TFtargte_with_effect_confict)
  # If in inferredGTEx, then use this effect
  with_coexpression = subset(ddf, evidence == 'inferredGTEx' & TFtarget %in% intersect(TFtargte_with_effect_confict, ddf$TFtarget[ ddf$evidence == 'inferredGTEx' ]) )
  # If NOT in inferredGTEx, then prioritize according database
  without_coexpression = subset(ddf, evidence == 'curated_databases' & TFtarget %in% setdiff(TFtargte_with_effect_confict, ddf$TFtarget[ ddf$evidence == 'inferredGTEx' ]) )
  without_coexpression = rbind(subset(without_coexpression, database == 'TFe_signed'),
                               subset(without_coexpression, database == 'tfact_signed'),
                               subset(without_coexpression, database == 'trrust_signed'),
                               subset(without_coexpression, database == 'trrd_via_tfact_signed'),
                               subset(without_coexpression, database == 'NFIRegulomeDB'),
                               subset(without_coexpression, database == 'oreganno_signed')  )
  # generate hash-maps
  with_coexpression = unique(with_coexpression[, c('TF', 'target', 'effect', 'TFtarget')])
  without_coexpression = unique(without_coexpression[ ! duplicated(without_coexpression$TFtarget) , c('TF', 'target', 'effect', 'TFtarget')])
  # fix effect sign for curated & signed
  idx = ddf$TFtarget %in% with_coexpression$TFtarget & ddf$evidence == 'curated_databases' & ddf$effect != 0
  ddf$effect[ idx ] = with_coexpression$effect[ match(ddf$TFtarget[ idx ], with_coexpression$TFtarget) ]
  idx = ddf$TFtarget %in% without_coexpression$TFtarget & ddf$evidence == 'curated_databases' & ddf$effect != 0
  ddf$effect[ idx ] = without_coexpression$effect[ match(ddf$TFtarget[ idx ], without_coexpression$TFtarget) ]
  return(ddf)
}


which_unsigned_duplicated = function(ddf){
  ddf$TFtargetEvidence = paste(ddf$TF, ddf$target, ddf$evidence)
  unsigned = ddf$effect == 0
  dup = ddf$TFtargetEvidence %in% ddf$TFtargetEvidence[duplicated(ddf$TFtargetEvidence)]
  which(unsigned & dup)
}


intersect_TFTG = function(df, evidences, match_sign = T, any_of = NULL){
  # filter data
  ddf = subset(df, evidence %in% evidences )
  ddf$TFtarget = paste(ddf$TF, ddf$target)
  # select overlapping tf-target interactions
  int_counts = table( unlist(sapply(evidences, function(e) unique(subset(ddf, evidence == e)$TFtarget) )) )
  if( ! is.null(any_of) ){
    interactions_of_interest = names(which(int_counts >= any_of))
  }else{
    interactions_of_interest = names(which(int_counts == length(evidences) ))
  }
  message(length(interactions_of_interest), ' potential interactions of interest')
  ddf = subset(ddf, TFtarget %in% interactions_of_interest)
  # order: prioritize evidence curated (curated > inferred) and signed effect (1/-1 > 0)
  ddf = ddf[ order(ddf$evidence), ]
  ddf = ddf[ order(abs(ddf$effect), decreasing = T), ]
  ddf = ddf[ order(ddf$TFtarget), ]
  if( ! match_sign ){
    new_ddf = unique(subset(ddf, TFtarget %in% interactions_of_interest)[, 1:2])
  }else{
    # remove evidences with unknown effects
    is_evidence_unknowneffect = sapply(evidences, function(e) all(subset(ddf, evidence == e)$effect == 0) )
    message('Is unknown effect:')
    print(is_evidence_unknowneffect)
    # if all unknown effect, then ignore sign
    if(all(is_evidence_unknowneffect))
      new_ddf = unique(subset(ddf, TFtarget %in% interactions_of_interest)[, 1:2])
    # if only 1 with knwon effect, then consider the sign of this evidence
    evidence_with_unknown_effect = names(which(is_evidence_unknowneffect))
    if( length(is_evidence_unknowneffect) - length(evidence_with_unknown_effect) == 1){
      new_ddf = unique(subset(ddf, TFtarget %in% interactions_of_interest & ! evidence %in% evidence_with_unknown_effect )[, 1:3])
      return(new_ddf)
    }
    # if more than 1 evidence with sign, select those evidences and find discrepant/opposite signed TFtarget interactions
    ddf = subset(ddf, ! evidence %in% evidence_with_unknown_effect )
    discrepant = which(table(  unique(subset(ddf, effect %in% c(1,-1) )[, -4 ])$TFtarget)    == 2)
    non_discrepant = setdiff(interactions_of_interest, names(discrepant) )
    ddf = subset(ddf, TFtarget %in% non_discrepant)
    new_ddf = ddf[ ! duplicated(ddf$TFtarget) , 1:3]
  }
  return(new_ddf)
}




# Load networks but randoms, cnio and consensus
networks = list.files('data/TF_target_sources/', recursive = T, pattern = 'network') %>% grep('sif', ., value=T) %>% grep('consensus', ., value = T, invert = T) %>% grep('omnipath_scores', ., value = T, invert = T) %>% grep('RANDOM_NETWORKS', ., invert = T, value = T) %>% grep('cnio', ., invert = T, value = T)
# also remove unsigned networks (when they have the signed version)
networks = grep('trrust/', networks, invert = T, value = T)  %>% grep('trrd_via_tfact/', ., invert = T, value = T)  %>% grep('TFe/', ., invert = T, value = T)  %>% grep('tfact/', ., invert = T, value = T)  %>% grep('oreganno/', ., invert = T, value = T)  
# Select network/parameters
# inferredGTEx = 3 tissues
networks = grep('inferredGTEx/network_10', networks, invert = T, value = T) %>% grep('inferredGTEx/network_2', ., invert = T, value = T) %>% grep('inferredGTEx/network_5', ., invert = T, value = T) %>% grep('inferredGTEx/network_dynamic', ., invert = T, value = T)
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
  if( is.null(regulon[[evidence]]) )
    regulon[[evidence]] = list()
  if( is.null(net$effect) )
    net$effect = 0
  for( tf in unique(net$TF) ){
    regulon[[evidence]][[tf]] = unique(rbind(regulon[[evidence]][[tf]],
                                             subset(net, TF == tf)[, c('target', 'effect') ] ))
  }
}
# build merged data.frame ~ database
df = unique(melt(regulon, id.vars = names(regulon$chipSeq[[1]]) ))
df = df[, c(3,1,2,4)]
names(df) = c('TF', 'target', 'effect', 'evidence')
write.table(df, file = 'data/TF_target_sources/merged_interactions.txt', quote = F, row.names = F, col.names = T, sep = '\t')




# Fix opposite/conflict effect sign in curated databases
df$TFtarget = paste(df$TF, df$target)
df$TFtarget_effect = paste(df$TF, df$target, df$effect)
df = fix_sign_conflict(df)



# remove unsigned & duplicated TF-targets from the same evidence
int2remove = which_unsigned_duplicated(df)
df = df[ - int2remove, ]



# BUILD CONSENSUS NETWORKS
# specific evidence pairs
unique(df$evidence)
cons2 = intersect_TFTG(df, evidences = c("chipSeq", "TFBS_scanning"), match_sign = F)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_chipSeq_TFBSscanning_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')



cons2 = intersect_TFTG(df, evidences = c("chipSeq", "inferredGTEx"), match_sign = T)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_chipSeq_inferredGTEx_signed_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("chipSeq", "inferredGTEx"), match_sign = F)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_chipSeq_inferredGTEx_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')



cons2 = intersect_TFTG(df, evidences = c("chipSeq", "curated_databases"), match_sign = T)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_chipSeq_curateddatabases_signed_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("chipSeq", "curated_databases"), match_sign = F)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_chipSeq_curateddatabases_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')




cons2 = intersect_TFTG(df, evidences = c("TFBS_scanning", "inferredGTEx"), match_sign = T)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_TFBSscanning_inferredGTEx_signed_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("TFBS_scanning", "inferredGTEx"), match_sign = F)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_TFBSscanning_inferredGTEx_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')



cons2 = intersect_TFTG(df, evidences = c("TFBS_scanning", "curated_databases"), match_sign = T)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_TFBSscanning_curateddatabases_signed_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("TFBS_scanning", "curated_databases"), match_sign = F)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_TFBSscanning_curateddatabases_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')



cons2 = intersect_TFTG(df, evidences = c("inferredGTEx", "curated_databases"), match_sign = T)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_curateddatabases_inferredGTEx_signed_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("inferredGTEx", "curated_databases"), match_sign = F)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_curateddatabases_inferredGTEx_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')






# any 2 evidences
cons2 = intersect_TFTG(df, evidences = c("inferredGTEx", "curated_databases", "chipSeq", "TFBS_scanning"), match_sign = F, any_of = 2)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("inferredGTEx", "curated_databases", "chipSeq", "TFBS_scanning"), match_sign = T, any_of = 2)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence2_signed_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')


# any 3 evidences
cons3 = intersect_TFTG(df, evidences = c("inferredGTEx", "curated_databases", "chipSeq", "TFBS_scanning"), match_sign = F, any_of = 3)
summarize_network(cons3)
write.table(cons3, file= 'data/TF_target_sources/consensus/evidence3_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons3 = intersect_TFTG(df, evidences = c("inferredGTEx", "curated_databases", "chipSeq", "TFBS_scanning"), match_sign = T, any_of = 3)
summarize_network(cons3)
write.table(cons3, file= 'data/TF_target_sources/consensus/evidence3_signed_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')



# 4 evidences
cons4 = intersect_TFTG(df, evidences = c("inferredGTEx", "curated_databases", "chipSeq", "TFBS_scanning"), match_sign = F)
summarize_network(cons4)
write.table(cons4, file= 'data/TF_target_sources/consensus/evidence4_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
# 4 evidences
cons4 = intersect_TFTG(df, evidences = c("inferredGTEx", "curated_databases", "chipSeq", "TFBS_scanning"), match_sign = T)
summarize_network(cons4)
write.table(cons4, file= 'data/TF_target_sources/consensus/evidence4_signed_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
