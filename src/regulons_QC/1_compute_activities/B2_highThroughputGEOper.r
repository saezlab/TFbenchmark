#' Compute activities from perturbed and basal gene expression in 2 cell lines
#' corresponding to the B2 benchmarck dataset (GEO accession ids: GSE31534 and GSE31912)
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code loads the regulons (in viper format; see "0_format_viper_regulons")
#' and an gene expression matrix containing zscore-transformed measurements (samples as columns; genes as rows)
#' to estimate TF activities using the aREA method from the VIPER R package (www.bioconductor.org/packages/release/bioc/html/viper.html)
#' Please, cite (Alvarez et al. Nature Genetics 2016)
#' 
#' This code calls the viper() function, which output is a TF activities matrix (NES values; with samples as columns nd TFs as rows)
#' NOTE viper() and msviper() functionw return different format and thus, B1 -B2&B3 are processed in different ways.


# Set enviroment
rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/matrix.r')



# Set folders
expression_signature_GSE31534 = 'data/regulons_QC/B2_perturbations/GSE31534/expression_signatures/A375_zscores.rdata'
expression_signature_GSE31912 = 'data/regulons_QC/B2_perturbations/GSE31912/expression_signatures/MCF7_zscores.rdata'
act_folder_GSE31534 = 'data/regulons_QC/B2_perturbations/GSE31534/activities/'
act_folder_GSE31912 = 'data/regulons_QC/B2_perturbations/GSE31912/activities/'
regulons_folder ='data/TF_target_sources/'



# Load gene expression data (zscores) and average duplicated genes
message('\n\nLoading zscore-transformed gene expression data')
E_pert_GSE31534 = get(load(expression_signature_GSE31534))
E_pert_GSE31534 = average_duplicated_genes(E_pert_GSE31534)
E_pert_GSE31912 = get(load(file = expression_signature_GSE31912))
E_pert_GSE31912 = average_duplicated_genes(E_pert_GSE31912)




# Set min gene set size
N = 4


# Compute activities for each regulon dataset
message('Computing activities')
regulons = list.files(regulons_folder, recursive = T, pattern = 'viperRegulon.rdata')
for (reg in regulons ){
  message(' - ', reg)
  # Format regulon-specific activities outfile
  act_file = reg %>% gsub('/', '.', .) %>% gsub('viperRegulon_', '', .) %>%  gsub('viperRegulon', '', .) %>% 
    gsub('sif', 'activities.rdata', .)  %>% gsub('\\.\\.', '.', .)  %>% gsub('_\\.', '.', .)
  # Load regulons
  regulon = load(paste(regulons_folder, reg, sep = '')) %>% get(.)
  # Perturbation 1 activities
  activities = viper(eset = E_pert_GSE31534, regulon = regulon, minsize = N, nes = T, method = 'none', eset.filter = F, pleiotropy = F, verbose = F)
  activities = cbind(Size = sapply(rownames(activities), function(tf) round(sum(regulon[[tf]]$likelihood), digits = 0)  ) , activities)
  save(activities, file = paste(act_folder_GSE31534, act_file, sep = ''))
  # # Perturbation 2 activities
  activities = viper(eset = E_pert_GSE31912, regulon = regulon, minsize = N, nes = T, method = 'none', eset.filter = F, pleiotropy = F, verbose = F)
  activities = cbind(Size = sapply(rownames(activities), function(tf) round(sum(regulon[[tf]]$likelihood), digits = 0)  ) , activities)
  save(activities, file = paste(act_folder_GSE31912, act_file, sep = ''))
}



