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




find.geneannot = function(platform){
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
  idx = grep('control', levels) 
  levels = c(levels[ ! 1:length(idx) %in% idx ], levels[ 1:length(idx) %in% idx ]) 
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
  cat('\n\n\n\n\n', rep('*', 50),'\nDownloading data for\n', my_GEOid, '\n')
  # Create study_files_directory
  study_files_directory = paste('/Users/luzgaral/tmp/GEO/', my_GEOid, '/', sep = '')
  dir.create(study_files_directory, showWarnings = F)
  samples_annot = NULL
  # Download corresponding files (soft and cel) from GEO in a directory
  if( length(grep('GSE', my_GEOid)) == 1 ){
    GSEDataobj = try(getGEO(my_GEOid, GSEMatrix = F, destdir = study_files_directory), silent = T)
    GEOdata_matrix = try(getGEO(my_GEOid, GSEMatrix = T, destdir = study_files_directory), silent = T)
    GSEid = Meta(GSEDataobj)$geo_accession
    GSMid = Meta(GSEDataobj)$sample_id
    platform_id = Meta(GSEDataobj)$platform_id
    if( ! class(GEOdata_matrix) == 'try-error' ){
      samples_annot = pData(phenoData(GEOdata_matrix[[1]]))[, 2:1]
    }else{
      samples_annot = t(sapply(Meta(GSEDataobj)$sample_id, function(x)  c(x, Meta( getGEO(x, destdir = study_files_directory))$title) ))
      samples_annot = as.data.frame(samples_annot)
    }
  }
  if( length(grep('GDS', my_GEOid)) == 1 ){
    GDSDataobj = try(getGEO(my_GEOid, GSEMatrix = F, destdir = study_files_directory), silent = T)
    GEOdata_matrix = try(getGEO(my_GEOid, GSEMatrix = T, destdir = study_files_directory), silent = T)
    GSEid = Meta(GDSDataobj)$reference_series
    platform_id = Meta(GDSDataobj)$platform
    samples_annot = Columns(GDSDataobj)
  }
  # Define samples 
  names(samples_annot)[1] = 'sample'
  names(samples_annot)[2] = 'phenotype'
  samples_annot$phenotype = gsub('rep',  '', gsub('replicate',  '', samples_annot$phenotype, ignore.case = T), ignore.case = T)
  samples_annot$phenotype = gsub('*[0-9]$', '', as.character(samples_annot$phenotype))
  samples_annot$phenotype = gsub(',', '_', as.character(samples_annot$phenotype))
  samples_annot$phenotype = gsub('-', '_', as.character(samples_annot$phenotype))
  samples_annot$phenotype = gsub(' ', '_', samples_annot$phenotype)
  samples_annot$phenotype = gsub('=', '_', samples_annot$phenotype)
  samples_annot$phenotype = gsub('__', '_', samples_annot$phenotype)
  samples_annot$phenotype = toupper(gsub('_$', '', samples_annot$phenotype))
  samples_annot$value = 1
  samples = as.character(samples_annot$sample[ samples_annot$phenotype %in% names(which(table(samples_annot$phenotype)>1)) ])
  rownames(samples_annot) = samples_annot$sample
  try(getGEOSuppFiles(GEO = GSEid, baseDir = study_files_directory), silent = T)
  cel_tarfile = list.files(study_files_directory, pattern = '_RAW.tar', full.names = T, ignore.case = T, recursive = T)
  eset = NULL
  # Load cel files and normalize
  if( length(cel_tarfile) > 0 ){
    untar(cel_tarfile, exdir = study_files_directory, extras = '-vz')
    cel_files = list.files(study_files_directory, pattern = 'cel', full.names = T, ignore.case = T)
    sampleNames = sapply(list.files(study_files_directory, pattern = 'cel', ignore.case = T), function(x) unlist(strsplit(unlist(strsplit(x, split = '_')), split = '\\.'))[1] )
    if( length(cel_files) != 0 ){
      my_cels = try(ReadAffy(filenames = cel_files[ sampleNames %in% samples ], #[ sampleNames %in% rownames(design) ], 
                           sampleNames = sampleNames[ sampleNames %in% samples ],#[ sampleNames %in% rownames(design) ], # clean names
                           verbose = F, 
                           compress = T), silent = T)
      if( class(my_cels) != 'try-error' )
        # Normalize cel files
        eset = try(affy::rma(my_cels), silent = T)
    }
  }
  if( is.null(eset) | class(eset) == 'try-error' ){
      no_cels <<- my_GEOid
      eset = try(GDS2eSet(GEOdata_matrix,do.log2=T), silent = T)
      if( class(eset) == 'try-error' )
        return(NULL)
      boxplot(exprs(eset), main = my_GEOid)
  }
  # Fix gene symbols in fData
  gene_annot = find.geneannot(platform_id)
  if( ! is.null(gene_annot) ){
    ID = featureNames(eset)
    Symbol = getSYMBOL(ID, data = gene_annot)
    fData(eset) = data.frame(ID=ID,Symbol=Symbol)
  }
  # Design phenotipic contrasts
  design = acast(samples_annot[ colnames(eset) , ], sample~phenotype, fill = 0)
  # Fit model with contrasts
  fit = lmFit(eset, design)
  RESULTS = list()
  pairwise_contrasts = design.pairs(colnames(design))
  for ( id in colnames(pairwise_contrasts) ){
    fit2 = contrasts.fit(fit, contrasts=pairwise_contrasts[ ,id])
    ebayes = eBayes(fit2, trend=TRUE)
    RESULTS[[id]] = topTable(ebayes, adjust="fdr", number = nrow(eset))
  }
  outfile = paste(study_files_directory, 'DEx.rdata', sep = '')
  save(RESULTS, file = outfile)
  cat('Done!\n', outfile,'\n', rep('*', 50),'\n\n')
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
write.csv(GEO_ids, file = "data/TF_target_sources/perturbations/Harmonizome/attribute_list_entries_GEOsignatures_TF_perturbation.csv", row.names = F)


no_cels = c()
manually_chechek_GEOs = c('GDS4280')
results = list()
pdf(file = 'data/TF_target_sources/perturbations/GEO/nocels_boxplots.pdf')
for (my_GEOid in rev(GEO_ids$GEOaccession ))
  results[[my_GEOid]] = run_differential_expression(my_GEOid)
dev.off()


