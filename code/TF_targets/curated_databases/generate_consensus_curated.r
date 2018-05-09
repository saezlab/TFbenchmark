rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


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
    # ddf = subset(ddf, ! evidence %in% evidence_with_unknown_effect )
    discrepant = which(table(  unique(subset(ddf, effect %in% c(1,-1) )[, -4 ])$TFtarget)    == 2)
    non_discrepant = setdiff(interactions_of_interest, names(discrepant) )
    ddf = subset(ddf, TFtarget %in% non_discrepant)
    new_ddf = ddf[ ! duplicated(ddf$TFtarget) , 1:3]
  }
  return(new_ddf)
}




# Load networks
networks = list.files('data/TF_target_sources/', recursive = T, pattern = 'network') %>% grep('sif', ., value=T) %>% grep('curated', ., value=T) %>% grep('consensus', ., value = T, invert = T) %>%
  grep('oreganno/', ., value = T, invert = T) %>% grep('tfact/', ., value = T, invert = T) %>% grep('TFe/', ., value = T, invert = T) %>% grep('trrust/', ., value = T, invert = T) %>% grep('trrd_via_tfact/', ., value = T, invert = T) %>% grep('cnio/', ., value = T, invert = T)

regulon = list()
for (n in networks){
  message(n)
  # Load network and fromat regulons
  net = read.delim(paste('data/TF_target_sources/', n, sep = ''), stringsAsFactors = F, sep = '\t')
  evidence = unlist(strsplit(n, '/'))[2]
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
df = unique(melt(regulon, id.vars = names(regulon[[1]][[1]]) ))
df = df[, c(3,1,2,4)]
names(df) = c('TF', 'target', 'effect', 'evidence')
head(df)


# remove unsigned & duplicated TF-targets from the same evidence
# int2remove = which_unsigned_duplicated(df)
# df = df[ - int2remove, ]



# any 1 evidences
cons1 = intersect_TFTG(df, evidences = names(regulon), match_sign = T, any_of = 1)
summarize_network(cons1)
write.table(cons1, file= 'data/TF_target_sources/consensus/evidence_1curateddatabases_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')


# any 2 evidences
cons2 = intersect_TFTG(df, evidences = names(regulon), match_sign = T, any_of = 2)
summarize_network(cons2)
write.table(cons2, file= 'data/TF_target_sources/consensus/evidence_2curateddatabases_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')



# any 3 evidences
cons3 = intersect_TFTG(df, evidences = names(regulon), match_sign = T, any_of = 3)
summarize_network(cons3)
write.table(cons3, file= 'data/TF_target_sources/consensus/evidence_3curateddatabases_network.sif', col.names = T, row.names = F, quote = F, sep = '\t')
