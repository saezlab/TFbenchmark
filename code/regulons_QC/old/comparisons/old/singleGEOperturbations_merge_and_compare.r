rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(plyr)



get_percentile_in_randoms = function(ob_nes, ra_nes){
  per = ecdf(ra_nes)(ob_nes)
}



load_and_filter_activities = function(f){
  gsea_results = get(load(f))
  gsea_results = subset(gsea_results, Size >= Nmin)
  gsea_results = subset(gsea_results, Size <= Nmax)
  gsea_results$Regulon = as.character(gsea_results$Regulon)
  gsea_results$regulon = sapply(strsplit(gsea_results$Regulon, split = ' - '), tail, 1)
  gsea_results$regulon[ grep('^r[1-9]', gsea_results$regulon) ] = 'random'
  return(gsea_results)
}
  

Nmax = 1000
Nmin = 4
activities_list = list()



activity_files = list.files(path = 'data/regulons_QC/pertubations/GEO/activities/real', pattern = 'viper', full.names = T)
for (f in activity_files){
  # Load data
  message(f)
  observed_results = load_and_filter_activities(f)
  observed_results = observed_results[ grep('aracne', observed_results$Regulon, invert = T), names(observed_results)  != 'Size_adaptative' ]
  random_results = load_and_filter_activities(gsub('real/', 'random_coexpressionPriors/', f))
  random_results$Regulon = 'random - random'
  random_results$perturbation_tissue = observed_results$perturbation_tissue[1]
  random_results$regulon = paste('random/random', random_results$Size, sep = '_')
  random_results = random_results[, names(observed_results)]
  observed_results = rbind(observed_results, subset(random_results, Size %in% c(100, 200, 500, 1000)) )
  observed_results = subset(observed_results, Size %in% random_results$Size)
  observed_results = observed_results[ order(observed_results$Size),  ]
  # Compute percentiles and ranks
  percentile = unlist(sapply(unique(observed_results$Size),
                    function(size)  get_percentile_in_randoms(subset(observed_results, Size == size)$NES, subset(random_results, Size == size)$NES))
  )
  observed_results$pvalue_1tailed = 1-percentile
  observed_results$pvalue_2tailed = 1-abs(percentile - 0.5)*2
  # observed_results$TF_rank_1tailed = rank(observed_results$pvalue_1tailed, ties.method = "first") / max(rank(observed_results$pvalue_1tailed, ties.method = "first"))
  # observed_results$TF_rank_2tailed = rank(observed_results$pvalue_2tailed, ties.method = "first") / max(rank(observed_results$pvalue_2tailed, ties.method = "first"))

  # Format output
  # 1) TF regulon name and tissue specificity
  observed_results$TF = unlist(sapply(strsplit(as.character(observed_results$Regulon), split = ' - '), head, 1))
  observed_results$regulon_is_tissuespecific = sapply(strsplit(observed_results$TF, split = ' '), length) > 1
  observed_results$regulon_tissue = 'none'
  observed_results$regulon_tissue[ observed_results$regulon_is_tissuespecific ] = sapply(strsplit(observed_results$TF[ observed_results$regulon_is_tissuespecific ], split = ' '), tail, 1)
  observed_results$TF[ observed_results$regulon_is_tissuespecific ] = sapply(strsplit(observed_results$TF[ observed_results$regulon_is_tissuespecific ], split = ' '), head, 1)

  # 2) Experiment parameters
  observed_results$experiment = gsub('data/regulons_QC/pertubations/GEO/activities/real/viper_', '', f) %>% gsub('.rdata', '', .)
  observed_results$perturbed_TF = sapply(strsplit(observed_results$experiment, split = '\\.'), head, 1)

  # 3) Regulon dataset names
  observed_results$regulon_type = sapply(strsplit(observed_results$regulon, split = '\\/'), head, 1)
  observed_results$regulon_shortname = sapply(strsplit(as.character(observed_results$regulon), split = '\\/'), function(x) x[2])
  observed_results$regulon_shortname[ observed_results$regulon_is_tissuespecific ] = 'inferredGTEx_tissue_specific'
  observed_results$regulon_shortname[ observed_results$regulon_shortname == 'scanning_output' ] = sapply(strsplit(observed_results$regulon[ observed_results$regulon_shortname == 'scanning_output' ], split = "\\/"), tail, 1)
  observed_results$regulon_shortname[ observed_results$regulon_shortname == 'ReMap' ] = sapply(strsplit(observed_results$regulon[ observed_results$regulon_shortname == 'ReMap' ], split = "\\/"), tail, 1)
  observed_results$regulon_shortname = gsub('network_', '', observed_results$regulon_shortname) %>% gsub('.rdata', '', .) %>% gsub('_network', '', .) %>% gsub('network', '', .)
  observed_results$regulon_shortname = gsub('viperRegulon_', '', observed_results$regulon_shortname) %>% gsub('_viperRegulon', '', .) %>% gsub('viperRegulon', '', .)
  observed_results$regulon_shortname[ observed_results$regulon_type == 'reverse_engenieered' ] = paste('>', observed_results$regulon_shortname[ observed_results$regulon_type == 'reverse_engenieered' ], 'tissues')

  activities_list[[f]] = observed_results
}
save(activities_list, file = 'data/regulons_QC/pertubations/GEO/results/viper_isolated_results_coexpressionPriors.rdata')
# df = melt(activities_list, id.vars = names(activities_list[[1]]))
# df$TF = unlist(df$TF)
# df$is_perturbed = df$TF == df$perturbed_TF
# df$is_tissue_match = T
# df$is_tissue_match[ df$regulon_is_tissuespecific ] = df$perturbation_tissue[ df$regulon_is_tissuespecific ] == df$regulon_tissue[ df$regulon_is_tissuespecific ] 
# save(df, file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors.rdata')
# df$TF_rank_1tailed = NA
# df$TF_rank_2tailed = NA
# for( re in sort(unique(df$Regulon[df$is_perturbed])) ){
#   message(re)
#   idx = df$Regulon == re 
#   df$TF_rank_1tailed[ idx ] = rank(subset(df, idx)$pvalue_1tailed, ties.method = "first") / max(rank(subset(df, idx)$pvalue_1tailed), ties.method = "first"))
#   df$TF_rank_2tailed[ idx ] = rank(subset(df, idx)$pvalue_2tailed, ties.method = "first") / max(rank(subset(df, idx)$pvalue_2tailed, ties.method = "first"))
# }
# save(df, file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors.rdata')
# 
# 

activities_list = list()
activity_files = list.files(path = 'data/regulons_QC/pertubations/GEO/activities/real_tissuespecific', pattern = 'viper', full.names = T)
for (f in activity_files){
  # Load data
  message(f)
  observed_results = load_and_filter_activities(f)
  observed_results = observed_results[ grep('aracne', observed_results$Regulon, invert = T), names(observed_results)  != 'Size_adaptative' ]
  random_results = load_and_filter_activities(gsub('real_tissuespecific/', 'random_tissuespecific_coexpressionPriors/', f))
  random_results$perturbation_tissue = observed_results$perturbation_tissue[1]
  random_results$regulon = paste('random/random', random_results$Size, sep = '_')
  observed_results = subset(observed_results, Size %in% random_results$Size)
  observed_results = observed_results[ order(observed_results$Size),  ]
  # Compute percentiles and ranks
  percentile = unlist(sapply(unique(observed_results$Regulon), 
                             function(r_name)  get_percentile_in_randoms(subset(observed_results, Regulon == r_name)$NES, random_results[grep(r_name, random_results$Regulon),]$NES)
  ))
  observed_results$pvalue_1tailed = 1-percentile
  observed_results$pvalue_2tailed = 1-abs(percentile - 0.5)*2
  observed_results$TF_rank_1tailed = rank(observed_results$pvalue_1tailed, ties.method = "first") / max(rank(observed_results$pvalue_1tailed, ties.method = "first"))
  observed_results$TF_rank_2tailed = rank(observed_results$pvalue_2tailed, ties.method = "first") / max(rank(observed_results$pvalue_2tailed, ties.method = "first"))
  
  # Format output
  # 1) TF regulon name and tissue specificity
  observed_results$TF = unlist(sapply(strsplit(as.character(observed_results$Regulon), split = ' - '), head, 1))
  observed_results$regulon_is_tissuespecific = sapply(strsplit(observed_results$TF, split = ' '), length) > 1
  observed_results$regulon_tissue = 'none'
  observed_results$regulon_tissue[ observed_results$regulon_is_tissuespecific ] = sapply(strsplit(observed_results$TF[ observed_results$regulon_is_tissuespecific ], split = ' '), tail, 1)
  observed_results$TF[ observed_results$regulon_is_tissuespecific ] = sapply(strsplit(observed_results$TF[ observed_results$regulon_is_tissuespecific ], split = ' '), head, 1)
  
  # 2) Experiment parameters 
  observed_results$experiment = gsub('data/regulons_QC/pertubations/GEO/activities/real/viper_', '', f) %>% gsub('.rdata', '', .)
  observed_results$perturbed_TF = sapply(strsplit(observed_results$experiment, split = '\\.'), head, 1)
  
  # 3) Regulon dataset names
  observed_results$regulon_type = sapply(strsplit(observed_results$regulon, split = '\\/'), head, 1)
  observed_results$regulon_shortname = sapply(strsplit(as.character(observed_results$regulon), split = '\\/'), function(x) x[2])
  observed_results$regulon_shortname[ observed_results$regulon_is_tissuespecific ] = 'inferredGTEx_tissue_specific'
  observed_results$regulon_shortname[ observed_results$regulon_shortname == 'scanning_output' ] = sapply(strsplit(observed_results$regulon[ observed_results$regulon_shortname == 'scanning_output' ], split = "\\/"), tail, 1)
  observed_results$regulon_shortname[ observed_results$regulon_shortname == 'ReMap' ] = sapply(strsplit(observed_results$regulon[ observed_results$regulon_shortname == 'ReMap' ], split = "\\/"), tail, 1)
  observed_results$regulon_shortname = gsub('network_', '', observed_results$regulon_shortname) %>% gsub('.rdata', '', .) %>% gsub('_network', '', .) %>% gsub('network', '', .)
  observed_results$regulon_shortname = gsub('viperRegulon_', '', observed_results$regulon_shortname) %>% gsub('_viperRegulon', '', .) %>% gsub('viperRegulon', '', .)
  observed_results$regulon_shortname[ observed_results$regulon_type == 'reverse_engenieered' ] = paste('>', observed_results$regulon_shortname[ observed_results$regulon_type == 'reverse_engenieered' ], 'tissues')
  observed_results = observed_results[order(observed_results$pvalue_1tailed), ]
  activities_list[[f]] = observed_results
}
df = melt(activities_list, id.vars = names(activities_list[[1]]))
df$TF = unlist(df$TF)
df$is_perturbed = df$TF == df$perturbed_TF
df$is_tissue_match = T
df$is_tissue_match[ df$regulon_is_tissuespecific ] = df$perturbation_tissue[ df$regulon_is_tissuespecific ] == df$regulon_tissue[ df$regulon_is_tissuespecific ] 
for( re in sort(unique(df$Regulon[df$is_perturbed])) ){
  idx = df$Regulon == re 
  df$TF_rank_1tailed[ idx ] = rank(subset(df, idx)$pvalue_1tailed, ties.method = "first") / max(rank(subset(df, idx)$pvalue_1tailed, ties.method = "first"))
  df$TF_rank_2tailed[ idx ] = rank(subset(df, idx)$pvalue_2tailed, ties.method = "first") / max(rank(subset(df, idx)$pvalue_2tailed, ties.method = "first"))
}
df$TF_rank_1tailed = NA
df$TF_rank_2tailed = NA
save(df, file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors_tissuespecific.rdata')

