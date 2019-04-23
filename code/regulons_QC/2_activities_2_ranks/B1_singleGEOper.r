#' Integration and formating of the TF activities files
#' from individual TF perturbation experiments
#' corresponding to the B1 benchmarck datasets
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code loads the individual TF activity files (tables; one file per PERTURBATION EXPERIMENT) 
#' and integrates them into:
#'    1. "activity_nes" a TF activities matrix (rows: TF regulon; columns: perturbtion experiments; values NES scores from VIPER)
#'    2. "activity_ranks" a matrix rank-trnasformed NES (i.e. for each TF regulon, we rank the perturbation experiments according the NES)
#'        the lower the NES (i.e. negative changes in TF activities), the lower the rank.
#'    3. "columns_perturbation_annot" data frame with relevant information from the regulons (TF name, dataset, evidence type, group, etc)
#'    4. "rows_regulons_annot" data frame with relevant information from the perturbations (perturbed TF, GEO id of the experiment, etc).
#' 
#' NOTE: Global and tissue specific regulons are processed separately because we haven't tested all the inferretdGTEx regulons
#'  in all perturbations, just those with mathcing tissue




rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/functions.r')





# Set regulon size thresolds
Nmax = 25000
Nmin = 5





# List activity files
activity_files = list.files(path = 'data/regulons_QC/B1_perturbations/activities/', pattern = 'rdata', full.names = T, recursive = T)
# Load activities data
activities_list = list()
for (f in activity_files){
  df = load_activity_msVIPERfile(f, Nmin, Nmax)
  experiment = strsplit(f, '/') %>% unlist(.) %>% tail(., 1) %>%
    gsub('.rdata', '', .) %>% strsplit(., '_') %>% unlist(.) %>% tail(., 1)
  df$Sample = experiment
  if( length(grep('noCOMBAT', f)) == 1 ){ #noCOMBAT analysis requested by a referee. Ignore it otherwise
    df$Regulon = paste( df$Regulon, '_noCOMBAT', sep ='') # This adds the "_noCOMBAT" label to the regulon name. 
    df$perturbation_tissue = paste( df$perturbation_tissue, '_noCOMBAT', sep ='')
  }
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


# Perturbation information - COLUMNS
perturbation = colnames(activity_nes)
columns_perturbation_annot = data.frame(perturbation_id = perturbation,
                                perturbed_TF = sapply(strsplit(perturbation, '\\.'), head, 1),
                                perturbation_GEOid = strsplit(perturbation, '\\.') %>% sapply(., tail, 1) %>% strsplit(., '-') %>%  sapply(., head, 1),
                                stringsAsFactors = F)

# Save data
save(activity_nes, 
     activity_ranks,
     columns_perturbation_annot, rows_regulons_annot,
     file = 'data/regulons_QC/B1_perturbations/results/aggregated_activities.rdata')





# ############################################################## 
# # Tissue or cancer specific regulons  
# ############################################################## 
# activities_list = list()
# activity_files = list.files(path = 'data/regulons_QC/B1_perturbations/activities/tissue_specific/', pattern = 'rdata', full.names = T)
# # Load  activities 
# activities_list = list()
# for (f in activity_files){
#   df = load_activity_msVIPERfile(f, Nmin, Nmax)
#   experiment = strsplit(f, '/') %>% unlist(.) %>% tail(., 1) %>% 
#     gsub('.rdata', '', .) %>% strsplit(., '_') %>% unlist(.) %>% tail(., 1)  
#   df$Sample = experiment
#   if( length(grep('noCOMBAT', f)) == 1 ){ #noCOMBAT analysis requested by a referee. Ignore it otherwise
#     df$Regulon = paste( df$Regulon, '_noCOMBAT', sep ='') # This adds the "_noCOMBAT" label to the regulon name. 
#     df$perturbation_tissue = paste( df$perturbation_tissue, '_noCOMBAT', sep ='')
#   }
#   activities_list[[f]] = df
# }
# 
# 
# # Generate activities matrixes integrating all the regulons
# activity_nes = aggregate_activities_NESmatrix(activities_list)
# # Generate ranks
# activity_ranks = activity_nes2ranks(activity_nes)
# 
# 
# 
# 
# # Regulons information - ROWS
# regulons = rownames(activity_nes)
# rows_regulons_annot = data.frame(regulon_id = regulons,
#                                  TF = sapply(strsplit(regulons, ' '), head, 1),
#                                  regulon_dataset_full = paste('inferred', sapply(strsplit(regulons, ' '), tail, 1)) ,
#                                  stringsAsFactors = F)
# rows_regulons_annot$regulon_evidence = 'inferred'
# rows_regulons_annot$regulon_group = manage_datasets_names(rows_regulons_annot$regulon_dataset_full, what = 'dataset2group')
# rows_regulons_annot$regulon_dataset = rows_regulons_annot$regulon_group
# 
# 
# 
# 
# 
# # Perturbation information - COLUMNS
# perturbation = colnames(activity_nes)
# columns_perturbation_annot = data.frame(perturbation_id = perturbation,
#                                         perturbed_TF = sapply(strsplit(perturbation, '\\.'), head, 1),
#                                         stringsAsFactors = F)
# 
# # Save data
# save(activity_nes, 
#      activity_ranks,
#      columns_perturbation_annot, rows_regulons_annot,
#      file = 'data/regulons_QC/B1_perturbations/results/viper_aggregated_results_tissuespecific.rdata')
# 
