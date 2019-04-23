#' Compute activities from differential gene expression signatures 
#' from individual TF perturbation experiments
#' corresponding to the B1 benchmarck datasets
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code loads the GTEx regulons (in viper format; see "code/TF_targets/inferred/" and  '0_format_viper_regulons')
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
GTEx_folder = 'data/TF_target_sources/inferred/GTEx/tissue_specific_noCOMBAT/'



# Set minimum regulon size
Nmin = 4



# Find perturbed TFs
perturbed_TFs = list.dirs(path = design_folder, full.names = F) %>% setdiff(., '')



# List GTEx networks
network_files = list.files(GTEx_folder, recursive = T, pattern = 'viperRegulon.rdata')



# Load perturbations, matching GTEx regulon and compute activities with msVIPER
perturbation_files = list.files(path = contrast_folder, pattern = 'rdata', full.names = T)
for (DE_file in perturbation_files){
  
  perturbation_id = tail(unlist(strsplit(DE_file, split = '/')),1) %>% gsub('.rdata', '', .) %>% strsplit(., split = '\\.') %>% unlist(.)
  design_file = list.files(design_folder, recursive = T, full.names = T) %>% grep(perturbation_id[1], ., value = T) %>% grep(perturbation_id[2], ., value = T)
  
  GTEx_regulon = find_matching_GTEx_regulon(design_file, GTEx_folder, perturbed_TFs)
  if( is.null(GTEx_regulon) )
    next()
  activities = DEs2activities(DE_file, design_folder, GTEx_regulon)
  
  analysis_name = gsub('contrasts/', 'activities/tissue_specific/viper_noCOMBAT_GTEx_', DE_file)
  save(activities, file = analysis_name)
  
}





