rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/functions.r')
source('code/lib/data.r')




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




# Load TF census
TFcensus = load_TFs_census()


# build merged data.frame ~ database
TFTG_information = load_and_merge_regulons_asdf(networks2merge)
TFTG_information = subset(TFTG_information, TF %in% TFcensus )
TFTG_information$TFtarget = paste(TFTG_information$TF, TFTG_information$target)
TFTG_information$TFtarget_effect = paste(TFTG_information$TF, TFTG_information$target, TFTG_information$effect)
# nrow(unique(TFTG_information[, 1:2]))
# length(unique(TFTG_information$TF))
# length(unique(TFTG_information$target))
# table(TFTG_information$evidence)


# identify unique interactions and add the SIGN
TFTG_information = fix_sign_conflict(TFTG_information)
TFTG_information$effect[ TFTG_information$evidence == 'inferred'  ] = 0 # Ignore sign from coexpression
omnipath_TFTG_database = assign_TFTG_sign(TFTG_information)
omnipath_TFTG_database$TFtarget = paste(omnipath_TFTG_database$TF, omnipath_TFTG_database$target)
omnipath_TFTG_database$TFtarget_effect = paste(omnipath_TFTG_database$TF, omnipath_TFTG_database$target, omnipath_TFTG_database$effect)
# head(omnipath_TFTG_database)



# add SCORE
TFTG_scores = read.csv(file = 'data/TF_target_sources/omnipath_scores/TF_target_scores.csv', stringsAsFactors = F)
omnipath_TFTG_database$score = TFTG_scores$score[ match(omnipath_TFTG_database$TFtarget, TFTG_scores$TFtarget) ]
# table(TFTG_scores$score)
# head(omnipath_TFTG_database)



# add evidences
omnipath_TFTG_database$is_evidence_curateddatabase = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, evidence == 'curated_databases' )$TFtarget
omnipath_TFTG_database$is_evidence_ChIPSeq = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, evidence == 'ChIP_Seq' )$TFtarget
omnipath_TFTG_database$is_evidence_TFbindingMotif = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, evidence == 'TFBS_scanning' )$TFtarget
omnipath_TFTG_database$is_evidence_coexpression = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, evidence == 'inferred' )$TFtarget
# head(omnipath_TFTG_database)
write.csv(omnipath_TFTG_database[, -c(4:5) ], file = 'data/TF_target_sources/omnipath_scores/database.csv', row.names = F)



# add sources database
# curateddatabase
omnipath_TFTG_database$which_curateddatabase = 'none'
idx = which(omnipath_TFTG_database$is_evidence_curateddatabase )
curated_TFTG_information = subset(TFTG_information, evidence == 'curated_databases')
omnipath_TFTG_database$which_curateddatabase[ idx ] = sapply(omnipath_TFTG_database$TFtarget[ idx ], function(TFTG)
  paste(subset(curated_TFTG_information, TFtarget == TFTG)$database, collapse = ',')
  )
omnipath_TFTG_database$which_curateddatabase = gsub('_signed', '', omnipath_TFTG_database$which_curateddatabase)
# ChIP_Seq
omnipath_TFTG_database$which_ChIPSeq = 'none'
idx = which(omnipath_TFTG_database$is_evidence_ChIPSeq )
omnipath_TFTG_database$which_ChIPSeq[ idx ] = 'ReMap'
# TFbindingMotif
omnipath_TFTG_database$which_TFbindingMotif = 'none'
idx.1 = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, database  == 'network_jaspar2018_n500.sif')$TFtarget
omnipath_TFTG_database$which_TFbindingMotif[ idx.1  ] = 'jaspar_v2018'
idx.2 = omnipath_TFTG_database$TFtarget %in% subset(TFTG_information, database  == 'network_hocomocov11_n500.sif')$TFtarget
omnipath_TFTG_database$which_TFbindingMotif[ idx.2  ] = 'hocomoco_v11'
omnipath_TFTG_database$which_TFbindingMotif[ idx.1 & idx.2  ] = 'hocomoco_v11,jaspar_v2018'
# coexpression
omnipath_TFTG_database$which_coexpression = 'none'
idx = which(omnipath_TFTG_database$is_evidence_coexpression )
omnipath_TFTG_database$which_coexpression[ idx ] = 'ARACNe-GTEx'
write.csv(omnipath_TFTG_database[, -c(4:5) ], file = 'data/TF_target_sources/omnipath_scores/database.csv', row.names = F)






## add pubmed IDS from curated data sources
data = list.files('data/TF_target_sources/curated_databases/', pattern = 'pubmed', full.names = T, recursive = T) %>%
  lapply(read.table, stringsAsFactors = F, header = T) %>% lapply(., function(x) x[, c("TF", "target", "pubmed")] )
mdata = unique(rbind(data[[1]], data[[2]], data[[3]], data[[4]], data[[5]], data[[6]], data[[7]]))
mdata$id = paste(mdata$TF, mdata$target)
xref = data.frame(id=unique(mdata$id), stringsAsFactors = F)
xref$pubmed = sapply(xref$id, function(id) paste(mdata$pubmed[ which(mdata$id == id) ], collapse = ',') )
xref$pubmed_count =  lapply(xref$pubmed, strsplit, ',') %>% lapply(., unlist) %>% lapply(., grep, pattern ='[0-9]', value = T) %>% lapply(., unique) %>% sapply(., length)
xref$pubmed =  lapply(xref$pubmed, strsplit, ',') %>% lapply(., unlist) %>% lapply(., grep, pattern ='[0-9]', value = T) %>% lapply(., unique) %>% sapply(., paste, collapse = ",")
xref = xref[ order(xref$pubmed_count, decreasing = T) ,]
table(xref$pubmed_count)
xref = subset(xref, pubmed_count > 0)
head(xref)

omnipath_TFTG_database$pubmedID_from_curated_resources = '-'
omnipath_TFTG_database$pubmedID_from_curated_resources[ omnipath_TFTG_database$TFtarget %in% xref$id ] = xref$pubmed[ match(omnipath_TFTG_database$TFtarget[ omnipath_TFTG_database$TFtarget %in% xref$id ],
                                                                                                              xref$id) ]
table(omnipath_TFTG_database$which_curateddatabase[ omnipath_TFTG_database$is_evidence_curateddatabase & omnipath_TFTG_database$pubmedID_from_curated_resources == 'none' ])
length(which(omnipath_TFTG_database$is_evidence_curateddatabase & omnipath_TFTG_database$pubmedID_from_curated_resources == 'none' ))
length(which(omnipath_TFTG_database$is_evidence_curateddatabase & omnipath_TFTG_database$pubmedID_from_curated_resources != 'none' ))
write.csv(omnipath_TFTG_database[, -c(4:5) ], file = 'data/TF_target_sources/omnipath_scores/database.csv', row.names = F)



## add kegg
kegg = read.delim('data/TF_target_sources/curated_databases/kegg/GErel_hsa.txt', stringsAsFactors = F)
kegg$id = paste(kegg$Entry1Name, kegg$Entry2Name)
kegg$pathway = paste(kegg$L1, gsub(' - Homo sapiens \\(human)', '', kegg$path_name), sep = '|')
xref = data.frame(id=unique(kegg$id), stringsAsFactors = F)
xref$pathway = sapply(xref$id, function(id) paste(kegg$pathway[ which(kegg$id == id) ], collapse = ',') )
xref$pathway_count =  lapply(xref$pathway, strsplit, ',') %>% lapply(., unlist) %>% lapply(., grep, pattern ='[0-9]', value = T) %>% lapply(., unique) %>% sapply(., length)
xref$pathway =  lapply(xref$pathway, strsplit, ',') %>% lapply(., unlist) %>% lapply(., grep, pattern ='[0-9]', value = T) %>% lapply(., unique) %>% sapply(., paste, collapse = ",")
table(xref$pathway_count)
xref = xref[ order(xref$pathway_count, decreasing = T) ,]
head(xref)

omnipath_TFTG_database$kegg_pathway = '-'
omnipath_TFTG_database$kegg_pathway[ omnipath_TFTG_database$TFtarget %in% xref$id ] = xref$pathway[ match(omnipath_TFTG_database$TFtarget[ omnipath_TFTG_database$TFtarget %in% xref$id ],
                                                                                                                            xref$id) ]
table(omnipath_TFTG_database$kegg_pathway)
write.csv(omnipath_TFTG_database[, -c(4:5) ], file = 'data/TF_target_sources/omnipath_scores/database.csv', row.names = F)

