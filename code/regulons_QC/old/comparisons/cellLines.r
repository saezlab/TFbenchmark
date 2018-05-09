rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


# Set regulon size range
N_min = 5
N_max = 1000




# Load unsigned activities
# load randoms
random_activity_files = list.files('data/regulons_QC/cell_lines/activities', pattern = 'RANDOM_NETWORKS_coexpressionPriors', full.names = T) %>% grep('_tissueSpecific', ., value =T, invert = T)
randoms = vector("list", length = N_max)
for (ra in random_activity_files){
  message(ra)
  nes = get(load(ra))
  nes = nes[ nes[,1] <= N_max,  ]
  for( si in unique(nes[,1]) ){
    randoms[[si]] = c(randoms[[si]], nes[ nes[,1] == si , -1])
  }
}
sapply(randoms, length)
names(randoms) = 1:N_max
randoms = randoms[ sapply(randoms, length) > 100 ]

# Normalize activities using randoms
activity_files = list.files('data/regulons_QC/cell_lines/activities', full.names = T) %>% grep(.,  pattern = 'RANDOM', value = T, invert= T) %>% grep(.,  pattern = 'aracne', value = T, invert= T)
activities_list = list()
for (f in activity_files){
  message(f)
  nes = get(load(f))
  nes = nes[ as.character(nes[,'Size']) %in% names(randoms) & nes[,'Size'] >= N_min , ]
  if( ! is.matrix(nes) )
    next
  if( nrow(nes) == 0)
    next
  if( is.matrix(nes) ) {
    p1t = nes[,-1]
    for( tf in rownames(nes) ){
      tf_size = as.character(nes[tf,'Size'])
      p1t[tf,] = get_percentile(p1t[tf,], randoms[[tf_size]])
    }
    df = melt(p1t, value.name = 'pvalue_1tailed')
    df$nes = melt(nes[ , colnames(p1t) ], value.name = 'nes')$nes
    activities_list[[f]] = df
  }
}




# Generate activities matrixes
df = melt(activities_list, id.vars = names(activities_list[[1]]))
df$value = df$nes
activity_nes = acast(df, Var1 ~ Var2, fill = NA)
df$value = df$pvalue_1tailed
activity_pvalue1tailed = acast(df, Var1 ~ Var2, fill = NA)



# Generate ranks
activity_rank_pvalue1tailed = t(apply(activity_pvalue1tailed, 1, rank, ties.method = 'min'))
activity_rank_pvalue1tailed = activity_rank_pvalue1tailed / apply(activity_rank_pvalue1tailed, 1, max)
activity_rank_pvalue1tailed[ is.na(activity_nes) ] = NA




# Regulons information - ROWS
regulons = rownames(activity_nes)
rows_regulons_annot = data.frame(regulon_id = regulons,
                                 TF = manage_datasets_names(regulons, what = 'regulon_id2TF'),
                                 regulon_dataset_full = manage_datasets_names(regulons, what = 'regulon_id2regulon_dataset_full'),
                                 stringsAsFactors = F)
rows_regulons_annot$regulon_evidence = manage_datasets_names(rows_regulons_annot$regulon_dataset_full, what = 'regulon_dataset_full2evidence')
rows_regulons_annot$regulon_dataset = manage_datasets_names(rows_regulons_annot$regulon_dataset_full, what = 'regulon_dataset_full2dataset')
rows_regulons_annot$regulon_group = manage_datasets_names(rows_regulons_annot$regulon_dataset, what = 'dataset2group')



# PErturbation information - COLUMNS
perturbation = colnames(activity_nes)
columns_perturbation_annot = data.frame(perturbation_id = perturbation,
                                        perturbed_TF = 'unknown',
                                        stringsAsFactors = F)

# Save data
save(activity_nes, 
     activity_pvalue1tailed, 
     activity_rank_pvalue1tailed, 
     columns_perturbation_annot, rows_regulons_annot,
     file = 'data/regulons_QC/cell_lines/activity_comparison_results.rdata')
