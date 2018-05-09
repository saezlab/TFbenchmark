home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
  
source('code/TFperturbations/contrast_lib.r')



# Read normalized expression files
norm_files = list.files('data/regulons_QC/pertubations/GEO/normalized', full.names = T)


# Read designs
desig_files = list.files('data/regulons_QC/pertubations/GEO/design', recursive = T, pattern = 'G', full.names = T)
designs = lapply(desig_files, read_desigfile)
names(designs) = sapply(designs, function(x) x$id)
length(designs)

outdir = 'data/regulons_QC/pertubations/GEO/contrasts/'

problematic_experiments = list()
for (de in designs){
  results = NULL
  print(de$id)
  norm_exp = try(get(load(grep(paste(de$accession, '.rdata', sep = ''), norm_files, value = T))), silent = T)
  if( class(norm_exp) == 'try-error')
    next
  try(results <- run_differential_expression(norm_exp, de, path_plots = outdir))
  if ( class(results) == "data.frame" )
    save(results, file = paste(outdir, de$id, '.rdata', sep = ''))
  if ( class(results) != "data.frame" ){
    cat('\n\n', rep('-', 20), '\nCheck experiment: ', de$id, '\n', rep('-', 20), '\n\n')
    problematic_experiments[[de$id]] = de
  }
}
length(problematic_experiments)
