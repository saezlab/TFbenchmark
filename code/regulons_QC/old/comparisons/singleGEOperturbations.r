rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(plyr)



load_and_filter_activities = function(f){
  gsea_results = get(load(f))
  gsea_results = subset(gsea_results, Size >= Nmin & Size <= Nmax)
  return(gsea_results)
}

Nmax = 1000
Nmin = 5



# Load UNSIGNED randoms
random_files = list.files(path = 'data/regulons_QC/pertubations/GEO/activities/random_coexpressionPriors', full.names = T)
randoms = vector("list", length = Nmax)
for (f in random_files){
  message(f)
  random_results = load_and_filter_activities(f)
  for( si in unique(random_results$Size) ){
    randoms[[si]] = c(randoms[[si]], subset(random_results, Size == si)$NES)
  }
}
names(randoms) = 1:Nmax



# Load UNSIGNED activities
activities_list = list()
activity_files = list.files(path = 'data/regulons_QC/pertubations/GEO/activities/real', pattern = 'viper', full.names = T)
for (f in activity_files){
  # Load data
  experiment = strsplit(f, '/') %>% unlist(.) %>%  tail(., 1) %>% gsub('viper_', '', .) %>% gsub('.rdata', '', .)
  message(experiment)
  observed_results = load_and_filter_activities(f)
  observed_results = subset(observed_results, Size %in% Nmin:Nmax)
  observed_results = observed_results[ order(observed_results$Size), ]
  # Compute percentiles and ranks
  percentile = unlist(sapply(unique(observed_results$Size),
                    function(size)
                      get_percentile(subset(observed_results, Size == size)$NES, randoms[[as.character(size)]])
  ))
  observed_results$pvalue_1tailed = 1-percentile
  activities_list[[experiment]] = observed_results
}


# Generate matrixes activities
df = melt(activities_list, id.vars = names(activities_list[[1]]))
df$value = df$NES
activity_nes = acast(df, Regulon ~ L1, fill = NA)
df$value = df$pvalue_1tailed
activity_pvalue1tailed = acast(df, Regulon ~ L1, fill = NA)
df$value = df$Size
regulon_size = acast(df, Regulon ~ L1, fill = NA)


# Generate ranks
activity_rank_pvalue1tailed = t(apply(activity_pvalue1tailed, 1, rank, ties.method = 'min'))
activity_rank_pvalue1tailed = activity_rank_pvalue1tailed / apply(activity_rank_pvalue1tailed, 1, max)
activity_rank_pvalue1tailed[ is.na(activity_nes) ] = NA



# CHECK SIZE bias
# x = melt(regulon_size)
# x$percentile = melt(activity_pvalue1tailed)$value
# xx = subset(x, value %in% c(5, 10, 50, 100, 500, 1000))
# xx$value = factor(xx$value, levels = c(5, 10, 50, 100, 500, 1000))
# ggplot(xx, aes(y = percentile, x = value )) + geom_boxplot() + xlab('regulon size')
# 




# Regulons information - ROWS
regulons = rownames(activity_nes)
rows_regulons_annot = data.frame(regulon_id = regulons,
                                TF = manage_datasets_names(regulons, what = 'regulon_id2TF'),
                                regulon_dataset_full = manage_datasets_names(regulons, what = 'regulon_id2regulon_dataset_full'),
                                stringsAsFactors = F)
rows_regulons_annot$regulon_evidence = manage_datasets_names(rows_regulons_annot$regulon_dataset_full, what = 'regulon_dataset_full2evidence')
rows_regulons_annot$regulon_dataset = manage_datasets_names(rows_regulons_annot$regulon_dataset_full, what = 'regulon_dataset_full2dataset')


# PErturbation information - COLUMNS
perturbation = colnames(activity_nes)
columns_perturbation_annot = data.frame(perturbation_id = perturbation,
                                perturbed_TF = sapply(strsplit(perturbation, '\\.'), head, 1),
                                perturbation_GEOid = strsplit(perturbation, '\\.') %>% sapply(., tail, 1) %>% strsplit(., '-') %>%  sapply(., head, 1),
                                stringsAsFactors = F)

# Save data
save(activity_nes, 
     activity_pvalue1tailed, 
     activity_rank_pvalue1tailed, 
     columns_perturbation_annot, rows_regulons_annot,
     file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors.rdata')





############################################################## 
# Signed 
############################################################## 
# NOTE: we need to do this separately because we haven't tested all the inferretdGTEx regulons in allperturbations
# just those with mathcing tissue
message('##########\nLoading tissue-specific networks ')

# Load SIGNED activities
activities_list = list()
activity_files = list.files(path = 'data/regulons_QC/pertubations/GEO/activities/real_tissuespecific', pattern = 'viper', full.names = T)
for (f in activity_files){
  # Load data
  experiment = strsplit(f, '/') %>% unlist(.) %>%  tail(., 1) %>% gsub('viper_', '', .) %>% gsub('.rdata', '', .)
  message(experiment)
  observed_results = load_and_filter_activities(f)
  observed_results = observed_results[ order(observed_results$Regulon),  ]
  random_results = load_and_filter_activities(gsub('real_tissuespecific', 'random_tissuespecific_coexpressionPriors', f))
  random_results$Regulon = sapply(strsplit(as.character(random_results$Regulon), split = ' - '), head, 1)
  
  # Compute percentiles and ranks
  percentile = unlist(sapply(unique(observed_results$Regulon),
                             function(r_name)  get_percentile(subset(observed_results, Regulon == r_name)$NES, 
                                                                         random_results[grep(r_name, random_results$Regulon),]$NES)
  ))
  observed_results$pvalue_1tailed = 1-percentile
  activities_list[[experiment]] = observed_results
}
# Generate matrixes activities
df = melt(activities_list, id.vars = names(activities_list[[1]]))
df$value = df$NES
activity_nes = acast(df, Regulon ~ L1, fill = NA)
df$value = df$pvalue_1tailed
activity_pvalue1tailed = acast(df, Regulon ~ L1, fill = NA)
# Generate ranks
activity_rank_pvalue1tailed = t(apply(activity_pvalue1tailed, 1, rank, ties.method = 'min'))
activity_rank_pvalue1tailed = activity_rank_pvalue1tailed / apply(activity_rank_pvalue1tailed, 1, max)
activity_rank_pvalue1tailed[ is.na(activity_nes) ] = NA



# Regulons information - ROWS
regulons = rownames(activity_nes)
rows_regulons_annot = data.frame(regulon_id = regulons,
                                 TF = sapply(strsplit(regulons, ' '), head, 1),
                                 regulon_dataset_full = paste('inferredGTEX', sapply(strsplit(regulons, ' '), tail, 1)) ,
                                 stringsAsFactors = F)
rows_regulons_annot$regulon_evidence = 'inferredGTEx'
rows_regulons_annot$regulon_dataset = rows_regulons_annot$regulon_dataset_full


# PErturbation information - COLUMNS
perturbation = colnames(activity_nes)
columns_perturbation_annot = data.frame(perturbation_id = perturbation,
                                        perturbed_TF = sapply(strsplit(perturbation, '\\.'), head, 1),
                                        stringsAsFactors = F)

# Save data
save(activity_nes, 
     activity_pvalue1tailed, 
     activity_rank_pvalue1tailed, 
     columns_perturbation_annot, rows_regulons_annot,
     file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors_tissuespecific.rdata')

