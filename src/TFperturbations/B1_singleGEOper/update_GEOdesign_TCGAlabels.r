#' This code modifies the design files to add a new field 
#' with the corresponding matching TCGA cancer type
#' if samples are cancer




# Set enviroment
rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/contrast.r')



# Define folders
design_folder = 'data/regulons_QC/B1_perturbations/design/'



# Read experimental design files
design_files = list.files(design_folder, recursive = T, pattern = 'G', full.names = T)
designs = lapply(design_files, read_desigfile)
names(designs) = sapply(designs, function(x) x$id)
sort(table(sapply(designs, function(x) x$GTEx_tissue)))


# Check each experiment
for (nd in names(designs )){
  d = designs[[nd]]
  print(d)
  print(new_designs[[nd]])
  n = readline(prompt="TCGA label: ")
  n = as.character(n)
  designs[[nd]]$TCGA_label = n
  save(designs, file = '~/tmp/designs.rdata')
}


sort(table(sapply(designs, function(x) x$GTEx_tissue)))
sort(table(sapply(designs, function(x) x$TCGA_label)))
any(is.na(sapply(designs, function(x) x$TCGA_label)))

for (nd in names(designs )){
  d = designs[[nd]]
  d$negative_samples = paste(d$negative_samples, collapse = ',')
  d$positive_samples = paste(d$positive_samples, collapse = ',')
  df = melt(d)[,2:1]
  filename_s = unlist(strsplit(d$id, '\\.'))[2]
  filename = paste(design_folder, d$pathway, '/', filename_s, '.txt', sep = '')
  print(filename)
  write.table(df, file = filename, col.names = F, row.names = F, sep = '\t', quote = F)
}
