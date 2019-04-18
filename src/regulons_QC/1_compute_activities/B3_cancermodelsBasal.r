#' Compute activities from basal gene expression in cancer cell lines
#' corresponding to the B3 benchmarck dataset (Garcia-Alonso et al. Cancer Research 2018)
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code loads the regulons (in viper format; see "0_format_viper_regulons")
#' and an gene expression matrix containing zscore-transformed basal measurements (samples as columns; genes as rows)
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



# Set folders
expression_signature_cell_lines = 'data/regulons_QC/B3_cell_lines/expression/voom_batchcor_duplicates_merged.RData'
act_folder = 'data/regulons_QC/B3_cell_lines/activities/'
regulons_folder = 'data/TF_target_sources/'




# Load expression data and compute zscores
message('\n\nLoading zscore-transformed gene expression data')
E_cell_lines = load(expression_signature_cell_lines) %>% get(.) %>% t(.) %>% scale(.) %>% t(.)
# Explude samples without COSMIC id since we do not have phenotypic data for these
load('data/regulons_QC/B3_cell_lines/cellines2tissues_mapping.rdata')
E_cell_lines = E_cell_lines[, colnames(E_cell_lines) %in% celllines_annotation$COSMIC_ID]




# Set min gene set size
N = 4



# Compute activities for each regulon dataset
message('Computing activities')
networks = list.files(regulons_folder, recursive = T, pattern = 'viperRegulon.rdata')
for (reg in networks ){
  message(' - ', reg)
  # Format regulon-specific activities outfile
  act_file = reg %>% gsub('/', '.', .) %>% gsub('viperRegulon_', '', .) %>%  gsub('viperRegulon', '', .) %>% 
    gsub('sif', 'activities.rdata', .)  %>% gsub('\\.\\.', '.', .)  %>% gsub('_\\.', '.', .)
  # Load regulons
  regulon = load(paste(regulons_folder, reg, sep = '')) %>% get(.)
  ## Cell lines activities
  activities = viper(eset = E_cell_lines, regulon = regulon, minsize = N, nes = T, method = 'none', eset.filter = F, pleiotropy = F, verbose = F)
  activities = cbind(Size = sapply(rownames(activities), function(tf) round(sum(regulon[[tf]]$likelihood), digits = 0)  ) , activities)
  save(activities, file = paste(act_folder, act_file, sep = ''))
}



