#' This code runs a differential expression analysis
#' using limma pipeline
#' for each GEO perturbation experiment




# Set enviroment
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/contrast.r')
source('code/lib/GEO.r')



# Define folders
DE_folder = 'data/regulons_QC/B1_perturbations/contrasts/'
norm_folder = 'data/regulons_QC/B1_perturbations/normalized/'
design_folder = 'data/regulons_QC/B1_perturbations/design'




# Read experimental design
desig_files = list.files(design_folder, recursive = T, pattern = 'G', full.names = T)
designs = lapply(desig_files, read_desigfile)
names(designs) = sapply(designs, function(x) x$id)
length(designs)



# List normalized expression files
norm_files = list.files(norm_folder, full.names = T)
# Run differential expression for each experiment
problematic_experiments = list()
for (de in designs){
  results = NULL
  print(de$id)
  norm_exp = try(get(load(grep(paste(de$accession, '.rdata', sep = ''), norm_files, value = T))), silent = T)
  if( class(norm_exp) == 'try-error')
    next
  try(results <- run_differential_expression(norm_exp, de, path_plots = DE_folder))
  if ( class(results) == "data.frame" )
    save(results, file = paste(DE_folder, de$id, '.rdata', sep = ''))
  if ( class(results) != "data.frame" ){
    cat('\n\n', rep('-', 20), '\nCheck experiment: ', de$id, '\n', rep('-', 20), '\n\n')
    problematic_experiments[[de$id]] = de
  }
}
length(problematic_experiments)
