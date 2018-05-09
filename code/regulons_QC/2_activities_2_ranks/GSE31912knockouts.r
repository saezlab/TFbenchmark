rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


N_min = 5
N_max = 20000


# Load unsigned activities
activity_files = list.files('data/regulons_QC/pertubations/GSE31912/activities', full.names = T) %>% grep(.,  pattern = 'RANDOM', value = T, invert= T) %>% grep(.,  pattern = 'aracne', value = T, invert= T)
activities_list = list()
for (f in activity_files){
  message(f)
  nes = get(load(f))
  nes = nes[ nes[,'Size'] <= N_max & nes[,'Size'] >= N_min , ]
  if( ! is.matrix(nes) )
    next
  if( nrow(nes) == 0)
    next
  if( is.matrix(nes) ) {
    df = melt(nes[,-1])
    names(df)[3] = 'nes'
    df$perturbed_gene = sapply(strsplit(as.character(df$Var2), split = '_'), tail, 1)
    activities_list[[f]] = df
  }
}




################################################ 
# Load signed activities
tissue = 'Breast'
f = list.files('data/regulons_QC/pertubations/GSE31912/activities', full.names = T) %>% grep(.,  pattern = tissue, value = T) %>% grep(.,  pattern = 'aracne', value = T)
message(f)
nes = get(load(f))
if( ! is.matrix(nes) )
  next
if( nrow(nes) == 0)
  next
if( is.matrix(nes) ) {
  df = melt(nes[,-1])
  names(df)[3] = 'nes'
  df$perturbed_gene = sapply(strsplit(as.character(df$Var2), split = '_'), tail, 1)
  activities_list[[f]] = df
}


# Generate matrixes activities
df = melt(activities_list, id.vars = names(activities_list[[1]]))
df$value = df$nes
activity_nes = acast(df, Var1 ~ Var2, fill = NA)

# Generate ranks
activity_ranks = t(apply(activity_nes, 1, rank, ties.method = 'min'))
activity_ranks = activity_ranks / apply(activity_ranks, 1, max)
activity_ranks[ is.na(activity_nes) ] = NA



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
     activity_ranks, 
     columns_perturbation_annot, rows_regulons_annot,
     file = 'data/regulons_QC/pertubations/GSE31912/activity_comparison_results.rdata')
