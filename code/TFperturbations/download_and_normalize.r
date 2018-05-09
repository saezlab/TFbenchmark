rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)


download_and_normalize_dataset = function(my_GEOid){
  cat('\n\nDownloading data for\n', my_GEOid, '\n')
  # Create study_files_directory
  study_files_directory = paste('/Users/luzgaral/tmp/GEO/', my_GEOid, '/', sep = '')
  dir.create(study_files_directory, showWarnings = F)
  # Download corresponding files (soft and cel) from GEO in a directory
  GEOobj = try(getGEO(my_GEOid, GSEMatrix = F, destdir = study_files_directory), silent = T)
  try(getGEOSuppFiles(my_GEOid, baseDir = study_files_directory), silent = T)
  metaData = Meta(GEOobj)
  platform_id = grep('^G', unname(metaData[ grep('platform', names(metaData)) ]), value = T)
  cel_tarfile = list.files(study_files_directory, pattern = '_RAW.tar', full.names = T, ignore.case = T, recursive = T)
  eset = NULL
  has_cells =  length(cel_tarfile) > 0
  # Load cel files and normalize
  if( ! class(cel_tarfile) == 'try-error' & has_cells  ){
    untar(cel_tarfile[1], exdir = study_files_directory, extras = '-xvz')
    cel_files = list.files(study_files_directory, pattern = 'cel', full.names = T, ignore.case = T)
    sampleNames = sapply(list.files(study_files_directory, pattern = 'cel', ignore.case = T), function(x) unlist(strsplit(unlist(strsplit(x, split = '_')), split = '\\.'))[1] )
    if( length(cel_files) > 1 ){
      cat('Loading CEL files\n')
      my_cels = try(ReadAffy(filenames = cel_files, 
                             sampleNames = sampleNames,
                             verbose = F, 
                             compress = T), silent = T)
      if( class(my_cels) != 'try-error' )
        # Normalize cel files
        eset = try(affy::rma(my_cels), silent = T)
    }
  }
  if( is.null(eset) | class(eset) == 'try-error' ){
    cat('Not CELs found. Downloading GSEMatrix\n')
    eset = try(ExpressionSet(exprs(getGEO(my_GEOid, destdir = study_files_directory)[[1]])), silent = T)
    if( class(eset) == 'try-error' ){
      GEOdata_matrix = try(getGEO(my_GEOid, GSEMatrix = T, destdir = study_files_directory), silent = T)
      eset = try(GDS2eSet(GEOdata_matrix, do.log2=T, getGPL = F), silent = T)
    }
    if( class(eset) == 'try-error' )
      eset = try(ExpressionSet(exprs(GEOdata_matrix[[1]])), silent = T)
    if( class(eset) == 'try-error' )
      return(NULL)
    cat('GSEMatrix loaded. Normalizing data using normalizeQuantiles() function\n')
    boxplot(exprs(eset), main = my_GEOid)
    eset = exprs(eset) %>%  normalizeQuantiles(.) %>% ExpressionSet(.)
    boxplot(exprs(eset), main = my_GEOid)
  }
  # Fix gene symbols in fData
  gene_annot = find.geneannot(platform_id)
  if( ! is.null(gene_annot) ){
    ID = featureNames(eset)
    Symbol = try(getSYMBOL(ID, data = gene_annot), silent = T)
    if( class(Symbol) != 'try-error' )
      fData(eset) = data.frame(ID=ID,Symbol=Symbol)
  }
  return(list(eset=eset, GEOobj=GEOobj, study_files_directory=study_files_directory, has_cells=has_cells))
}





source('code/TFperturbations/old/download_and_DGE_lib_experimental_designs.r')
GEOdatasets = read.csv(file = "data/regulons_QC/pertubations/attribute_list_entries_GEOsignatures_TF_perturbation.csv", stringsAsFactors = F)
GEOdatasets = subset(GEOdatasets, IS_TF)
nrow(GEOdatasets)


for (my_GEOid in sort(unique(GEOdatasets$GEOaccession), decreasing = T) ){
  data = try(download_and_normalize_dataset(my_GEOid), silent = T)
  if( ! is.null(data) & class(data) != "" & class(data) != 'try-error'){
    outfile = paste('data/regulons_QC/pertubations/GEO/normalized/', my_GEOid, '.rdata', sep = '')
    save(data, file = outfile)
  }
  if ( class(data) == 'try-error' | is.null(data)  )
    cat('\n\nCheck experiment ', my_GEOid, '\n\n')
}
df = subset(GEOdatasets, GEOaccession %in% gsub('.rdata', '', list.files('data/regulons_QC/pertubations/GEO/normalized/')) )
nrow(df)
write.csv(df, file = "data/regulons_QC/pertubations/attribute_list_entries_GEOsignatures_TF_perturbation_final.csv", row.names = F)








