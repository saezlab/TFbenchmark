#' Compute activities from differential gene expression signatures 
#' from individual TF perturbation experiments
#' corresponding to the B1 benchmarck datasets
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code loads the regulons (in viper format; see "0_format_viper_regulons")
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




# Set enviroment folders
contrast_folder = 'data/regulons_QC/B1_perturbations/contrasts'
design_folder = 'data/regulons_QC/B1_perturbations/design'
regulons_folder ='data/TF_target_sources'



# Set minimum regulon size
Nmin = 4


# Identify which TFs are perturbed
perturbed_TFs = list.dirs(path = design_folder, full.names = F) %>% setdiff(., '')




message('Merging regulons into one object')
network_files = list.files(regulons_folder, recursive = T, pattern = 'viperRegulon.rdata') %>% 
  grep("specific", ., invert = T, value = T)
networks = load_and_merge_regulons(regulon_files, regulons_folder, filter_TFs = perturbed_TFs)




message('Computing differential activities from DE signatures')
perturbation_files = list.files(path = contrast_folder, pattern = 'rdata', full.names = T)
for (DE_file in perturbation_files){
  activities = DEs2activities(DE_file, design_folder, networks)
  analysis_name = gsub('contrasts/', 'activities/global/viper_', DE_file)
  save(activities, file = analysis_name)
}





