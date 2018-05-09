rm(list = ls())
home = '~/Google Drive/projects/TFbenchmark/'
setwd(home)

source('code/download_and_process_TF_perturbations/download_and_DGE_lib_experimental_designs.r')



# Load GEO accessions from harmonizome
GEOdatasets = read.csv(file = "data/TF_target_sources/perturbations/Harmonizome/attribute_list_entries_GEOsignatures_TF_perturbation.csv", stringsAsFactors = F)



ANALYSIS = list()
pdf(file = 'data/TF_target_sources/perturbations/GEO/nocels_boxplots.pdf')
for (study_id in rev(GEOdatasets$study) ){
  message(study_id)
  my_GEOid = GEOdatasets$GEOaccession[ GEOdatasets$study == study_id ]
  ANALYSIS[[study_id]] = try(load_eset_and_phAnnot(my_GEOid), silent = T)
  if( ! is.null(ANALYSIS[[study_id]]) ){
    ANALYSIS[[study_id]]$GEOid = my_GEOid
    ANALYSIS[[study_id]]$name = GEOdatasets$study[ GEOdatasets$study == study_id ]
    ANALYSIS[[study_id]]$TF = GEOdatasets$GeneSym[ GEOdatasets$study == study_id ]
    # Run DEA
    dea = run_differential_expression(eset = ANALYSIS[[study_id]]$eset,
                                      sAnnot = ANALYSIS[[study_id]]$sAnnot,
                                      study_id = study_id)
    ANALYSIS[[study_id]]$dea =  dea$dea
    ANALYSIS[[study_id]]$contrast = dea$contrast
    # Save
    outfile = paste(ANALYSIS[[study_id]]$study_files_directory, 'DEx.rdata', sep = '')
    my_analysis = ANALYSIS[[study_id]]
    save(my_analysis, file = outfile)
    cat('Done!\n', outfile,'\n', rep('*', 50),'\n\n')
  }
}
dev.off()
save(ANALYSIS, file='data/TF_target_sources/perturbations/geo_overview.rdata')


