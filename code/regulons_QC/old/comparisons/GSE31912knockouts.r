rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


N = 5


# Load unsigned activities
# build unsigned randoms
random_activity_files = list.files('data/regulons_QC/pertubations/GSE31912/activities', pattern = 'RANDOM_NETWORKS_coexpressionPriors', full.names = T) %>% grep('tissueSpecific', ., value =T, invert = T)
randoms = vector("list", length = 1000)
for (ra in random_activity_files){
  message(ra)
  nes = get(load(ra))
  nes = nes[ nes[,1] <= 1000,  ]
  for( si in unique(nes[,1]) ){
    randoms[[si]] = c(randoms[[si]], nes[ nes[,1] == si , -1])
  }
}
sapply(randoms, length)
names(randoms) = 1:1000
randoms = randoms[ sapply(randoms, length) > 100 ]
# normalize unsigned activities
activity_files = list.files('data/regulons_QC/pertubations/GSE31912/activities', full.names = T) %>% grep(.,  pattern = 'RANDOM', value = T, invert= T) %>% grep(.,  pattern = 'aracne', value = T, invert= T)
activities_list = list()
for (f in activity_files){
  message(f)
  nes = get(load(f))
  nes = nes[ as.character(nes[,'Size']) %in% names(randoms) & nes[,'Size'] >= N , ]
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




################################################ 
# Load signed activities
tissue = 'Breast'
random_activity_file = list.files('data/regulons_QC/pertubations/GSE31912/activities', pattern = 'RAN', full.names = T) %>% grep(tissue, ., value = T)
nes = get(load(random_activity_file))
nes = nes[ nes[,1] <= 1000,  ]
nes[,1] = sapply(strsplit(rownames(nes), split = ' - '), head, 1)
randoms = list()
for( si in unique(nes[,1]) ){
  randoms[[si]] = as.numeric(nes[ nes[,1] == si, -1])
}
# normalize signed activities
activity_files = list.files('data/regulons_QC/pertubations/GSE31912/activities', full.names = T) %>% grep(.,  pattern = tissue, value = T) %>% grep(.,  pattern = 'aracne', value = T)
for (f in activity_files){
  message(f)
  nes = get(load(f))
  nes = nes[ rownames(nes) %in%  names(randoms) , ]
  if( ! is.matrix(nes) )
    next
  if( nrow(nes) == 0)
    next
  if( is.matrix(nes) ) {
    p1t = nes[,-1]
    for( tf in rownames(nes) ){
      p1t[tf,] = get_percentile(p1t[tf,], randoms[[tf]])
    }
    df = melt(p1t, value.name = 'pvalue_1tailed')
    df$Var1 = as.character(df$Var1)
    df$nes = melt(nes[ , colnames(p1t) ], value.name = 'nes')$nes
    df$perturbed_gene = sapply(strsplit(as.character(df$Var2), split = '_'), tail, 1)
    activities_list[[f]] = df
  }
}





# Generate matrixes activities
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
rows_regulons_annot$regulon_evidence[ rows_regulons_annot$regulon_dataset_full == tissue] = 'inferredGTEx'
rows_regulons_annot$regulon_dataset[ rows_regulons_annot$regulon_dataset_full == tissue] =  paste('inferredGTEx', tissue) 
rows_regulons_annot$regulon_group = manage_datasets_names(rows_regulons_annot$regulon_dataset, what = 'dataset2group')



# PErturbation information - COLUMNS
perturbation = colnames(activity_nes)
columns_perturbation_annot = data.frame(perturbation_id = perturbation,
                                        perturbed_TF = sapply(strsplit(perturbation, '_'), tail, 1),
                                        stringsAsFactors = F)

# Save data
save(activity_nes, 
     activity_pvalue1tailed, 
     activity_rank_pvalue1tailed, 
     columns_perturbation_annot, rows_regulons_annot,
     file = 'data/regulons_QC/pertubations/GSE31912/activity_comparison_results.rdata')
