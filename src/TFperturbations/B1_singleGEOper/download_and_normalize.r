#' This code retrieves the GEO perturbation experiments 
#' curated in data/regulons_QC/B1_perturbations/attribute_list_entries_GEOsignatures_TF_perturbation.csv
#' and normalizes the expression datasets
#' 
#' normalized data is saved in the normalized_folder folder




# Set enviroment
rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/GEO.r')
source('code/lib/functions.r')


# Define folders
experiments2download = "data/regulons_QC/B1_perturbations/attribute_list_entries_GEOsignatures_TF_perturbation.csv"
normalized_folder = 'data/regulons_QC/B1_perturbations/normalized/'
experiments2download_success = "data/regulons_QC/B1_perturbations/attribute_list_entries_GEOsignatures_TF_perturbation_final.csv"
tmp_directory = '/Users/luzgaral/tmp/GEO/'


# Read experiments from file
GEOdatasets = read.csv(file = experiments2download, stringsAsFactors = F)
GEOdatasets = subset(GEOdatasets, IS_TF)



# Load and normalize each experiment
for (my_GEOid in sort(unique(GEOdatasets$GEOaccession), decreasing = T) ){
  data = try(download_and_normalize_dataset(my_GEOid, tmp_directory), silent = T)
  if( ! is.null(data) & class(data) != "" & class(data) != 'try-error'){
    outfile = paste(normalized_folder, my_GEOid, '.rdata', sep = '')
    save(data, file = outfile)
  }
  if ( class(data) == 'try-error' | is.null(data)  )
    cat('\n\nCheck experiment ', my_GEOid, '\n\n')
}



# Check experiments successfully downloaded
df = subset(GEOdatasets, GEOaccession %in% gsub('.rdata', '', list.files()) )
write.csv(df, file = "data/regulons_QC/B1_perturbations/attribute_list_entries_GEOsignatures_TF_perturbation_final.csv", row.names = F)








