#' Integration and formating of the TF activities files
#' corresponding to the B2 benchmarck dataset GSE31534 
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#'
#'#' This code loads the individual TF activity files (tables; one file per REGULON DATASET) 
#' and integrates them into:
#'    1. "activity_nes" a TF activities matrix (rows: TF regulons; columns: perturbtion experiments; values NES scores from VIPER)
#'    2. "activity_ranks" a matrix rank-transformed NES (i.e. for each TF regulon, we rank the perturbation experiments according the NES)
#'        the lower the NES (i.e. lower TF activities), the lower the rank.
#'    3. "columns_perturbation_annot" data frame with relevant information from the regulons (TF name, dataset, evidence type, group, etc)
#'    4. "rows_regulons_annot" data frame with relevant information from the perturbations (perturbed TF, GEO id of the experiment, etc).
#' 



rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/functions.r')



# Set regulon size thresolds
N_min = 5
N_max = 20000


#  Activity files - global regulons
activity_files_g = list.files('data/regulons_QC/B2_perturbations/GSE31534/activities', full.names = T, pattern = 'rdata') %>% 
  grep(.,  pattern = 'specific', value = T, invert= T)
# Activity files - tissue specific regulons
tissue = 'Skin'
tcga = 'skcm'
activity_files_ts = c(list.files('data/regulons_QC/B2_perturbations/GSE31534/activities', full.names = T, pattern = 'rdata') %>% 
                        grep(.,  pattern = tissue, value = T),
                      list.files('data/regulons_QC/B2_perturbations/GSE31534/activities', full.names = T, pattern = 'rdata') %>% 
                        grep(.,  pattern = tcga, value = T))
activity_files = c(activity_files_g, activity_files_ts)




# Load  activities 
activities_list = list()
for (f in activity_files){
  df = load_activity_VIPERfile(f, Nmin, Nmax)
  if(is.null(df))
    next
  df$perturbed_gene = sapply(strsplit(as.character(df$Sample), split = '_'), tail, 1)
  activities_list[[f]] = df
}



# Generate activities matrixes integrating all the regulons
activity_nes = aggregate_activities_NESmatrix(activities_list)
# Generate ranks
activity_ranks = activity_nes2ranks(activity_nes)







# Regulons information - ROWS
regulons = rownames(activity_nes)
rows_regulons_annot = annotate_regulons(regulons)
length(regulons)


# PErturbation information - COLUMNS
perturbation = colnames(activity_nes)
columns_perturbation_annot = data.frame(perturbation_id = perturbation,
                                        perturbed_TF = sapply(strsplit(perturbation, '_'), tail, 1),
                                        stringsAsFactors = F)

# Save data
save(activity_nes, 
     activity_ranks, 
     columns_perturbation_annot, rows_regulons_annot,
     file = 'data/regulons_QC/B2_perturbations/GSE31534/activity_comparison_results.rdata')

