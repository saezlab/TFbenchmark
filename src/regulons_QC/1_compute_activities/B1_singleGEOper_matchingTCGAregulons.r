#' Compute activities from differential gene expression signatures 
#' from individual TF perturbation experiments
#' corresponding to the B1 benchmarck datasets 
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code loads the TCGA regulons (in viper format; see "code/TF_targets/inferred/download_aracne_networks_tcga.r")
#' and the differential gene expression tables generated with limma R package (lmFit > contrasts.fit > eBayes > topTable)
#' to estimate differential TF activities using the aREA method from the VIPER R package (www.bioconductor.org/packages/release/bioc/html/viper.html)
#' Please, cite (Alvarez et al. Nature Genetics 2016)
#' 
#' It geneares a data.frame containing TF activities (NES values), together with other experimental information




rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/functions.r')



# Set folders enviroment
contrast_folder = 'data/regulons_QC/B1_perturbations/contrasts'
design_folder = 'data/regulons_QC/B1_perturbations/design'
TCGA_folder = 'data/TF_target_sources/inferred/TCGA/cancer_specific/'



# Set minimum regulon size
Nmin = 4



# Find perturbed TFs
perturbed_TFs = list.dirs(path = design_folder, full.names = F) %>% setdiff(., '')



# List TCGA networks
network_files = list.files(TCGA_folder, recursive = T, pattern = 'viperRegulon.rdata')



# Load perturbations, matching TCGA regulon and compute activities with msVIPER
perturbation_files = list.files(path = contrast_folder, pattern = 'rdata', full.names = T)
for (DE_file in perturbation_files){
  
  perturbation_id = tail(unlist(strsplit(DE_file, split = '/')),1) %>% gsub('.rdata', '', .) %>% strsplit(., split = '\\.') %>% unlist(.)
  design_file = list.files(design_folder, recursive = T, full.names = T) %>% grep(perturbation_id[1], ., value = T) %>% grep(perturbation_id[2], ., value = T)
  
  TCGA_regulon = find_matching_TCGA_regulon(design_file, TCGA_folder, perturbed_TFs)
  if( is.null(TCGA_regulon) )
    next()
  activities = DEs2activities(DE_file, design_folder, TCGA_regulon)
  
  analysis_name = gsub('contrasts/', 'activities/tissue_specific/viper_tcga_', DE_file)
  save(activities, file = analysis_name)
  
}





