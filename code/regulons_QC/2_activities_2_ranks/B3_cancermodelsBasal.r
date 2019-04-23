#' Integration and formating of the TF activities files
#' corresponding to the B3 benchmarck dataset GSE31912
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




# Set regulon size range
N_min = 5
N_max = 25000



#  Activity files - global regulons
activity_files_g = list.files('data/regulons_QC/B3_cell_lines/activities', full.names = T, pattern = "rdata") %>% 
  grep(.,  pattern = 'specific', value = T, invert= T)
# Activity files - tissue specific regulons --- NOTE only those mathcing our samples
load('data/regulons_QC/B3_cell_lines/cellines2tissues_mapping.rdata')
GTEx_labels = setdiff(celllines_annotation$GTEx, 'none')
tcga_labels = setdiff(celllines_annotation$Study.Abbreviation, '')
activity_files_ts = list.files('data/regulons_QC/B3_cell_lines/activities', full.names = T, pattern = "rdata") %>% 
  grep(.,  pattern = 'specific', value = T)
activity_files_ts = sapply(c(GTEx_labels, tcga_labels), function(lab) grep(paste('\\.', lab,  sep = ''), activity_files_ts, ignore.case = T, value = T)) %>% 
  unlist(.) %>% unique(.)
activity_files = c(activity_files_g, activity_files_ts)



# Load  activities 
activities_list = list()
for (f in activity_files){
  df = load_activity_VIPERfile(f, Nmin, Nmax)
  if(is.null(df))
    next
  tissue = intersect(c(GTEx_labels, tolower(tcga_labels)), unlist(strsplit(f,'\\.')))
  if( length(tissue) == 1 ){ #  if tissue-specific regulon; then match regulon to cell-type
    which_samples_in_tissue = which(celllines_annotation$GTEx %in% tissue | 
                                    celllines_annotation$Study.Abbreviation %in% toupper(tissue))
    cell_lines_in_tissue =  celllines_annotation$COSMIC_ID[ which_samples_in_tissue ]
    cell_lines_in_tissue =  setdiff(cell_lines_in_tissue, NA) # Remove samples without cosmic id since we do not have phenotypic data for these
    df = subset(df, Sample %in% cell_lines_in_tissue  ) # remove samples not matching the regulon tissue-type
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


# Samples information - COLUMNS
# Note that 'perturbed' information will be extracted from CNA and shRNA phenotypic data
perturbation = colnames(activity_nes)
columns_perturbation_annot = data.frame(perturbation_id = perturbation,
                                        perturbed_TF = 'unknown',
                                        stringsAsFactors = F)

# Save data
save(activity_nes, 
     activity_ranks, 
     columns_perturbation_annot, rows_regulons_annot,
     file = 'data/regulons_QC/B3_cell_lines/aggregated_activities.rdata')

