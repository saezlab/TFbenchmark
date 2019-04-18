#' Generate and format TF regulons for the TF resources benchmark 
#' available at: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code loads the individual TF-target networks in .sif format
#' and transforms them into R objects with a format compatible with VIPER regulons (www.bioconductor.org/packages/release/bioc/html/viper.html)
#' Each row indicates an individual TF-target interaction
#' Columns indicate TF (mandatory), target (mandatory), effect (optional) and weight (optional), respectively
#' The supplied expression object has to contain HGNC symbols in rows.
#' 
#' 
#' The processing includes the following steps:
#' 1. Remove duplicated interactions
#' 2. Remove self interactions
#' 3. Remove emply interactions
#' 4. Remove overrepresented targets (i.e. targets regulateb by more than X TFs; here >10% of the TFs in the dataset). 
#'    Set "remove_overrepresented_targets" and "overrepresented_targets_th" variables.
#' 5. Exclude TFs not in the TFcensus
#' 6. Remove TFs with less than Nmin targets
#' 7. Define the TF/interaction mode (i.e. effect of the regulation: activation, repression)
#' 8. Define the TF/interaction weigth 
#' 9. Build VIPER regulon object





### Set enviroment
rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/functions.r')
source('code/lib/data.r')



### Variables
regulons_folder = 'data/TF_target_sources/'
Nmin = 4  # minimum number of targets per TF required
remove_overrepresented_targets = T # logical, remove overrepresented targets
overrepresented_targets_th = 0.1 # % of TFs allowed to regulate a target
repressors = load_repressors() # list of global repressors (i.e. TFs with gene expression repressor role)
TFcensus = load_TFs_census() # list of proteins considered to have a TF role




# Generate global regulons
network_files = list.files(regulons_folder, recursive = T, pattern = 'sif') %>% grep('network', ., value=T) 
for (n in network_files){
  viper_regulon = sif2viper(n, regulons_folder, TFcensus, repressors)
  save(viper_regulon, file = paste(regulons_folder, gsub('.sif', '_viperRegulon.rdata', n) , sep = '') )
}



