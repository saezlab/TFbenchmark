rm(list = ls())
home = '~/Google Drive/projects/TFbenchmark/'
setwd(home)

library(reshape2)
library(limma)
library(GEOquery)
library(GEOsearch)
library(biomaRt)
library(GEOmetadb)
library(affy)
library(hgu133plus2.db) # GPL570
library(hgu133a.db)  # GPL96
library(hgu133a2.db) # GPL571
library(illuminaHumanv4.db) # GPL10558
library(pd.hugene.1.0.st.v1) # GPL6244
library(annotate)




find.geneannot = function(GDSDataobj){
  platform = Meta(GDSDataobj)$platform
  gene_annot = NULL
  if( platform == 'GPL570')
    gene_annot = 'hgu133plus2.db'
  if( platform == 'GPL96')
    gene_annot = 'hgu133a.db'
  if( platform == 'GPL571')
    gene_annot = 'hgu133a2.db'
  if( platform == 'GPL10558')
    gene_annot = 'illuminaHumanv4.db'
  if( platform == 'GPL6244')
    gene_annot = 'pd.hugene.1.0.st.v1'
  return(gene_annot)
}


design.pairs = function(levels) {
  n <- length(levels)
  design <- matrix(0,n,choose(n,2))
  rownames(design) <- levels
  colnames(design) <- 1:choose(n,2)
  k <- 0
  for (i in 1:(n-1))
    for (j in (i+1):n) {
      k <- k+1
      design[i,k] <- 1
      design[j,k] <- -1
      colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
    }
  design
}


run_differential_expression = function(my_GEOid){
  cat('\n\n', rep('*', 50),'\nDownloading data for\n', my_GEOid, '\n')
  # Create study_files_directory
  study_files_directory = paste('data/TF_target_sources/perturbations/GEO/', my_GEOid, '/', sep = '')
  dir.create(study_files_directory, showWarnings = F)
  # Download corresponding files (soft and cel) from GEO in a directory
  GDSDataobj = try(getGEO(my_GEOid, GSEMatrix = TRUE, destdir = study_files_directory), silent = T)
  if( class(GDSDataobj) == 'try-error' )
    return(NULL)
  GSEid = Meta(GDSDataobj)$reference_series
  cel_tarfile = rownames(getGEOSuppFiles(GEO = GSEid, baseDir = study_files_directory))[1]
  if( is.null(cel_tarfile) )
    return(NULL)
  untar(cel_tarfile, exdir = study_files_directory)
  # Design phenotipic contrasts
  samples_annot = Columns(GDSDataobj)
  names(samples_annot)[2] = 'phenotype'
  samples_annot$phenotype = gsub(' ', '_', samples_annot$phenotype)
  samples_annot$value = 1
  design = acast(samples_annot, sample~phenotype, fill = 0)
  # Load cel files and normalize
  cel_files = list.files(study_files_directory, pattern = 'cel', full.names = T, ignore.case = T)
  sampleNames = sapply(list.files(study_files_directory, pattern = 'cel', ignore.case = T), function(x) unlist(strsplit(unlist(strsplit(x, split = '_')), split = '\\.'))[1] )
  if( length(cel_files) == 0 )
    return(NULL)
  my_cels = try(ReadAffy(filenames = cel_files[ sampleNames %in% rownames(design) ], verbose = F, 
                     sampleNames = sampleNames[ sampleNames %in% rownames(design) ], # clean names
                     compress = T), silent = T)
  if( class(my_cels) == 'try-error' )
    return(NULL)
  # Normalize cel files
  eset = affy::rma(my_cels)
  # Fix gene symbols in fData
  ID = featureNames(eset)
  gene_annot = find.geneannot(GDSDataobj)
  if(is.null(gene_annot))
    return(NULL)
  Symbol = getSYMBOL(ID, data = gene_annot)
  fData(eset) = data.frame(ID=ID,Symbol=Symbol)
  # Fit model with contrasts
  fit = lmFit(eset, design[colnames(exprs(eset)), ])
  RESULTS = list()
  pairwise_contrasts = design.pairs(colnames(design))
  for ( id in colnames(pairwise_contrasts) ){
    fit2 = contrasts.fit(fit, contrasts=pairwise_contrasts[,id])
    ebayes = eBayes(fit2, trend=TRUE)
    RESULTS[[id]] = topTable(ebayes, adjust="fdr", number = nrow(eset))
  }
  save(RESULTS, file = paste(study_files_directory, 'DEx.rdata'))
  cat('Done!\n', rep('*', 50),'\n\n')
  return(colnames(pairwise_contrasts))
}

# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# annot = getBM(attributes=c("hgnc_symbol", "entrezgene", 'ensembl_gene_id', 'gene_biotype', 'ensembl_transcript_id', 'uniprotswissprot'), mart = ensembl)
# 

# Load GEO accessions from harmonizome
GEO_ids = read.delim('data/TF_target_sources/perturbations/Harmonizome/attribute_list_entries_GEOsignatures_TF_perturbation.txt', stringsAsFactors = F)
GEO_ids$GeneSym = toupper(GEO_ids$GeneSym)
GEO_ids = GEO_ids[ grep('human', GEO_ids[,1]), ]
GEO_ids = GEO_ids[ order(GEO_ids$GeneSym), ]
GEO_ids$GEOaccession = sapply(GEO_ids[,1], function(x) toupper(tail(unlist(strsplit(x, split = '_')), 1)) )
GEO_ids$GEOplatform = sapply(GEO_ids[,1], function(x) toupper(tail(unlist(strsplit(x, split = '_')), 2))[1] )
results = list()
for (my_GEOid in GEO_ids$GEOaccession ){
  results[[my_GEOid]] = run_differential_expression(my_GEOid)
}




