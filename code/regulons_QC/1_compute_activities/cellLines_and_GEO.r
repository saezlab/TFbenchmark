rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/project_data_loads.r')
library(viper)
source('code/utils.r')


# Load data
# E_cell_lines = t(scale(t(load_expression_gdsc(RNASEQ = T, norm_method = 'none'))))
E_pert_GSE31534 = get(load(file = 'data/regulons_QC/pertubations/GSE31534/expression_signatures/A375_zscores.rdata'))
E_pert_GSE31534 = average_duplicated_genes(E_pert_GSE31534)
E_pert_GSE31912 = get(load(file = 'data/regulons_QC/pertubations/GSE31912/expression_signatures/MCF7_zscores.rdata'))
E_pert_GSE31912 = average_duplicated_genes(E_pert_GSE31912)




# Set min gene set size
N = 4


# Load networks
networks = list.files('data/TF_target_sources/', recursive = T, pattern = 'viperRegulon.rdata') %>% 
  grep('2', ., value = T, invert = F) %>% 
  grep('tissueSpecific', ., value = T, invert = T) %>% 
   grep('RAND', . , value = T, invert = F) #%>% 
  # grep('omni', ., value = T, invert = F)
message('\n\nComputing activities')
for (n in networks ){
  message(n)
  # Load network and fromat regulons
  regulon = get(load(paste('data/TF_target_sources/', n, sep = '')))
  outfile = n %>% gsub('/', '.', .) %>% gsub('viperRegulon_', '', .) %>%  gsub('viperRegulon', '', .) %>% gsub('sif', 'activities.rdata', .)  %>% gsub('\\.\\.', '.', .)  %>% gsub('_\\.', '.', .)
  
  # Perturbation 1 activities
  activities = viper(eset = E_pert_GSE31534, regulon = regulon, minsize = N, nes = T, method = 'none', eset.filter = F, pleiotropy = F, verbose = T)
  activities = cbind(Size = sapply(rownames(activities), function(tf) round(sum(regulon[[tf]]$likelihood), digits = 0)  ) , activities)
  save(activities, file = paste('data/regulons_QC/pertubations/GSE31534/activities/', outfile, sep = ''))
  # Perturbation 2 activities
  activities = viper(eset = E_pert_GSE31912, regulon = regulon, minsize = N, nes = T, method = 'none', eset.filter = F, pleiotropy = F, verbose = T)
  activities = cbind(Size = sapply(rownames(activities), function(tf) round(sum(regulon[[tf]]$likelihood), digits = 0)  ) , activities)
  save(activities, file = paste('data/regulons_QC/pertubations/GSE31912/activities/', outfile, sep = ''))
  # # Cell lines activities
  # activities = viper(eset = E_cell_lines, regulon = regulon, minsize = N, nes = T, method = 'none', eset.filter = F, pleiotropy = F, verbose = T)
  # activities = cbind(Size = sapply(rownames(activities), function(tf) round(sum(regulon[[tf]]$likelihood), digits = 0)  ) , activities)
  # save(activities, file = paste('data/regulons_QC/cell_lines/activities/', outfile, sep = ''))
  # message('\n\n\n')
}



