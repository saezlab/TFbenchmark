rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')



fix_sign_conflict = function (ddf){
  # Identify TF-target interactions whith effect sign conflict (i.e. 1 and -1)
  x = unique(subset(ddf, effect != 0 & evidence == 'curated_databases')[, c('TF', 'target', 'effect', 'TFtarget')])
  TFtargte_with_effect_confict = unique(x$TFtarget[ duplicated(x$TFtarget) ])
  message('NOTE: ', length(TFtargte_with_effect_confict), ' whith effect sign conflict.\nPrioritizing those with matching inferredGTEx effect sign')
  print(TFtargte_with_effect_confict)
  # If in inferredGTEx, then use this effect
  with_coexpression = subset(ddf, evidence == 'inferred' & TFtarget %in% intersect(TFtargte_with_effect_confict, ddf$TFtarget[ ddf$evidence == 'inferred' ]) )
  # If NOT in inferredGTEx, then prioritize according database
  without_coexpression = subset(ddf, evidence == 'curated_databases' & TFtarget %in% setdiff(TFtargte_with_effect_confict, ddf$TFtarget[ ddf$evidence == 'inferred' ]) )
  # without_coexpression = rbind(subset(without_coexpression, database == 'TFe_signed'),
  #                              subset(without_coexpression, database == 'tfact_signed'),
  #                              subset(without_coexpression, database == 'trrust_signed'),
  #                              subset(without_coexpression, database == 'trrd_via_tfact_signed'),
  #                              subset(without_coexpression, database == 'NFIRegulomeDB'),
  #                              subset(without_coexpression, database == 'oreganno_signed')  )
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
# inferredGTEx = 3 tissues
networks = list.files('data/TF_target_sources/inferred/GTEx/pantissue',  pattern = 'network_i3.sif', recursive = T, full.names = T)
networks2merge = append(networks2merge, values = networks)
# TFBS_scanning size = 500
networks = list.files('data/TF_target_sources/TFBS_scanning', pattern = 'n500.sif', recursive = T, full.names = T)
networks2merge = append(networks2merge, values = networks)
# ChIPSeq size = 500
networks = list.files('data/TF_target_sources/ChIP_Seq', pattern = 'n500.sif', recursive = T, full.names = T)
networks2merge = append(networks2merge, values = networks)
networks2merge




# Load and merge
regulon = list()
for (n in networks2merge){
  message(n)
  # Load network and fromat regulons
  net = read.delim(n, stringsAsFactors = F, sep = '\t')
  evidence = unlist(strsplit(n, '/'))[3]
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
df = unique(melt(regulon, id.vars = names(regulon$ChIP_Seq[[1]]) ))
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
cons2 = intersect_TFTG(df, evidences = c("ChIP_Seq", "TFBS_scanning"), match_sign = F)
write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_ChIPSeq_TFBSscanning.sif', col.names = T, row.names = F, quote = F, sep = '\t')



# cons2 = intersect_TFTG(df, evidences = c("ChIP_Seq", "inferred"), match_sign = T)
# write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_ChIPSeq_inferredGTEx_signed.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("ChIP_Seq", "inferred"), match_sign = F)
write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_ChIPSeq_inferredGTEx.sif', col.names = T, row.names = F, quote = F, sep = '\t')



# cons2 = intersect_TFTG(df, evidences = c("ChIP_Seq", "curated_databases"), match_sign = T)
# write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_ChIPSeq_curateddatabases_signed.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("ChIP_Seq", "curated_databases"), match_sign = F)
write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_ChIPSeq_curateddatabases.sif', col.names = T, row.names = F, quote = F, sep = '\t')




# cons2 = intersect_TFTG(df, evidences = c("TFBS_scanning", "inferred"), match_sign = T)
# write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_TFBSscanning_inferredGTEx_signed.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("TFBS_scanning", "inferred"), match_sign = F)
write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_TFBSscanning_inferredGTEx.sif', col.names = T, row.names = F, quote = F, sep = '\t')



# cons2 = intersect_TFTG(df, evidences = c("TFBS_scanning", "curated_databases"), match_sign = T)
# write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_TFBSscanning_curateddatabases_signed.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("TFBS_scanning", "curated_databases"), match_sign = F)
write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_TFBSscanning_curateddatabases.sif', col.names = T, row.names = F, quote = F, sep = '\t')



# cons2 = intersect_TFTG(df, evidences = c("inferred", "curated_databases"), match_sign = T)
# write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_curateddatabases_inferredGTEx_signed.sif', col.names = T, row.names = F, quote = F, sep = '\t')
cons2 = intersect_TFTG(df, evidences = c("inferred", "curated_databases"), match_sign = F)
write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_curateddatabases_inferredGTEx.sif', col.names = T, row.names = F, quote = F, sep = '\t')






# any 2 evidences
cons2 = intersect_TFTG(df, evidences = c("inferred", "curated_databases", "ChIP_Seq", "TFBS_scanning"), match_sign = F, any_of = 2)
write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2.sif', col.names = T, row.names = F, quote = F, sep = '\t')
# cons2 = intersect_TFTG(df, evidences = c("inferred", "curated_databases", "ChIP_Seq", "TFBS_scanning"), match_sign = T, any_of = 2)
# write.table(cons2, file= 'data/TF_target_sources/consensus/network_e2_signed.sif', col.names = T, row.names = F, quote = F, sep = '\t')


# any 3 evidences
cons3 = intersect_TFTG(df, evidences = c("inferred", "curated_databases", "ChIP_Seq", "TFBS_scanning"), match_sign = F, any_of = 3)
write.table(cons3, file= 'data/TF_target_sources/consensus/network_e3.sif', col.names = T, row.names = F, quote = F, sep = '\t')
# cons3 = intersect_TFTG(df, evidences = c("inferred", "curated_databases", "ChIP_Seq", "TFBS_scanning"), match_sign = T, any_of = 3)
# write.table(cons3, file= 'data/TF_target_sources/consensus/network_e3_signed.sif', col.names = T, row.names = F, quote = F, sep = '\t')



# 4 evidences
cons4 = intersect_TFTG(df, evidences = c("inferred", "curated_databases", "ChIP_Seq", "TFBS_scanning"), match_sign = F)
write.table(cons4, file= 'data/TF_target_sources/consensus/network_e4.sif', col.names = T, row.names = F, quote = F, sep = '\t')
# # 4 evidences
# cons4 = intersect_TFTG(df, evidences = c("inferred", "curated_databases", "ChIP_Seq", "TFBS_scanning"), match_sign = T)
# write.table(cons4, file= 'data/TF_target_sources/consensus/network_e4_signed.sif', col.names = T, row.names = F, quote = F, sep = '\t')
