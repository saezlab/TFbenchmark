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
  # fix effect sign
  idx = ddf$TFtarget %in% with_coexpression$TFtarget & ddf$evidence == 'curated_databases' & ddf$effect != 0
  ddf$effect[ idx ] = with_coexpression$effect[ match(ddf$TFtarget[ idx ], with_coexpression$TFtarget) ]
  idx = ddf$TFtarget %in% without_coexpression$TFtarget & ddf$evidence == 'curated_databases' & ddf$effect != 0
  ddf$effect[ idx ] = without_coexpression$effect[ match(ddf$TFtarget[ idx ], without_coexpression$TFtarget) ]
  return(ddf)
}





assign_TFTG_sign = function(ddf){
  ddf$TFtarget = paste(ddf$TF, ddf$target)
  # order: prioritize evidence curated (curated > inferred) and signed effect (1/-1 > 0)
  ddf = ddf[ order(ddf$database, decreasing = T), ]
  ddf = ddf[ order(ddf$evidence), ]
  ddf = ddf[ order(abs(ddf$effect), decreasing = T), ]
  ddf = ddf[ order(ddf$TFtarget), ]
  new_ddf = ddf[ ! duplicated(ddf$TFtarget) , 1:3]
  return(new_ddf)
}




# Load networks but randoms, cnio and consensus
networks = list.files('data/TF_target_sources/', recursive = T, pattern = 'network') %>% grep('sif', ., value=T) %>% grep('consensus', ., value = T, invert = T) %>% grep('omnipath_scores', ., value = T, invert = T) %>% grep('omnipath_scores', ., value = T, invert = T)  %>% grep('RANDOM_NETWORKS', ., invert = T, value = T) %>% grep('cnio', ., invert = T, value = T)
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

# Load TF census
TFcensus = load_TFs_census()


# build merged data.frame ~ database
TFTG_information = unique(melt(regulon, id.vars = names(regulon[[1]][[1]][[1]]) ))
TFTG_information = TFTG_information[, c(3,1,2,5,4)]
names(TFTG_information) = c('TF', 'target', 'effect', 'evidence', 'database')
TFTG_information = subset(TFTG_information, TF %in% TFcensus )
TFTG_information$TFtarget = paste(TFTG_information$TF, TFTG_information$target)
TFTG_information$TFtarget_effect = paste(TFTG_information$TF, TFTG_information$target, TFTG_information$effect)
nrow(unique(TFTG_information[, 1:2]))
length(unique(TFTG_information$TF))
length(unique(TFTG_information$target))


# identify unique interactions and add the SIGN
TFTG_information = fix_sign_conflict(TFTG_information)
TFTG_information$effect[ TFTG_information$evidence == 'inferredGTEx'  ] = 0 # Ignore sign from coexpression
omnipath_TFTG_database = assign_TFTG_sign(TFTG_information)
omnipath_TFTG_database$TFtarget = paste(omnipath_TFTG_database$TF, omnipath_TFTG_database$target)
omnipath_TFTG_database$TFtarget_effect = paste(omnipath_TFTG_database$TF, omnipath_TFTG_database$target, omnipath_TFTG_database$effect)
head(omnipath_TFTG_database)


# add SCORE
TFTG_scores = read.csv(file = 'data/TF_target_sources/omnipath_scores/TF_target_scores.csv', stringsAsFactors = F)
table(TFTG_scores$score)
omnipath_TFTG_database$score = TFTG_scores$score[ match(omnipath_TFTG_database$TFtarget, TFTG_scores$TFtarget) ]
head(omnipath_TFTG_database)



# add evidences
omnipath_TFTG_database$is_evidence_curateddatabase = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, evidence == 'curated_databases' )$TFtarget
omnipath_TFTG_database$is_evidence_chipSeq = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, evidence == 'chipSeq' )$TFtarget
omnipath_TFTG_database$is_evidence_TFbindingMotif = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, evidence == 'TFBS_scanning' )$TFtarget
omnipath_TFTG_database$is_evidence_coexpression = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, evidence == 'inferredGTEx' )$TFtarget
head(omnipath_TFTG_database)
write.csv(omnipath_TFTG_database[, -c(4:5) ], file = 'data/TF_target_sources/omnipath_scores/database.csv', row.names = F)



# add sources
#
omnipath_TFTG_database$which_curateddatabase = 'none'
idx = which(omnipath_TFTG_database$is_evidence_curateddatabase )
curated_TFTG_information = subset(TFTG_information, evidence == 'curated_databases')
omnipath_TFTG_database$which_curateddatabase[ idx ] = sapply(omnipath_TFTG_database$TFtarget[ idx ], function(TFTG)
  paste(subset(curated_TFTG_information, TFtarget == TFTG)$database, collapse = ',')
  )
#
omnipath_TFTG_database$which_chipSeq = 'none'
idx = which(omnipath_TFTG_database$is_evidence_chipSeq )
omnipath_TFTG_database$which_chipSeq[ idx ] = 'ReMap'
#
omnipath_TFTG_database$which_TFbindingMotif = 'none'
idx.1 = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, database  == 'network_jaspar2018_n500.sif')$TFtarget
omnipath_TFTG_database$which_TFbindingMotif[ idx.1  ] = 'jaspar_v2018'
idx.2 = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, database  == 'network_hocomocov11_n500.sif')$TFtarget
omnipath_TFTG_database$which_TFbindingMotif[ idx.2  ] = 'hocomoco_v11'
omnipath_TFTG_database$which_TFbindingMotif[ idx.1 & idx.2  ] = 'hocomoco_v11,jaspar_v2018'
#
omnipath_TFTG_database$which_coexpression = 'none'
idx = which(omnipath_TFTG_database$is_evidence_coexpression )
omnipath_TFTG_database$which_coexpression[ idx ] = 'ARACNe-GTEx'
write.csv(omnipath_TFTG_database[, -c(4:5) ], file = 'data/TF_target_sources/omnipath_scores/database.csv', row.names = F)






