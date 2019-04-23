

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
library(dplyr)




# --------------  
# -------------- Download and normalize an GEO experiment according to their GEO id
# --------------
# input: my_GEOid
# returns: 
#     eset = expression matrix (eset format)
#     GEOobj = corresponfing R GEO object
#     study_files_directory = files location
#     has_cells = logical

download_and_normalize_dataset = function(my_GEOid, tmp_directory='~/tmp/'){
  cat('\n\nDownloading data for\n', my_GEOid, '\n')
  # Create study_files_directory
  study_files_directory = paste(tmp_directory, my_GEOid, '/', sep = '')
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



format_phenoname = function(pheno_name){
  pheno_name = gsub(',', '_', as.character(pheno_name))
  pheno_name = gsub('-', '_', as.character(pheno_name))
  pheno_name = gsub(' ', '_', pheno_name)
  pheno_name = gsub('=', '_', pheno_name)
  pheno_name = gsub('rep',  '', gsub('replicate',  '', pheno_name, ignore.case = T), ignore.case = T)
  pheno_name = gsub('_[0-9]$', '', as.character(pheno_name))
  pheno_name = gsub('_[A-G]$', '', as.character(pheno_name))
  pheno_name = gsub('__', '_', pheno_name)
  pheno_name = toupper(gsub('_$', '', pheno_name))
  pheno_name = gsub('_[1-9]_', '', as.character(pheno_name))
  pheno_name = gsub('_[A-S,U-Z]_', '', as.character(pheno_name))
}




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
  if( platform == 'GPL10558')
    gene_annot = 'illuminaHumanv4.db'
  if( platform == 'GPL8300')
    gene_annot = 'hgu95av2.db'
  return(gene_annot)
}

design.pairs = function(levels) {
  idx = grep('pos_sample', levels) 
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


load_eset_and_phAnnot = function(my_GEOid, tmp_directory='~/tmp'){
  cat('\n\n\n\n\n', rep('*', 50),'\nDownloading data for\n', my_GEOid, '\n')
  # Create study_files_directory
  study_files_directory = paste(tmp_directory, my_GEOid, '/', sep = '')
  dir.create(study_files_directory, showWarnings = F)
  sAnnot = NULL
  # Download corresponding files (soft and cel) from GEO in a directory
  if( length(grep('GSE', my_GEOid)) == 1 ){
    GSEDataobj = try(getGEO(my_GEOid, GSEMatrix = F, destdir = study_files_directory), silent = T)
    GEOdata_matrix = try(getGEO(my_GEOid, GSEMatrix = T, destdir = study_files_directory), silent = T)
    GSEid = Meta(GSEDataobj)$geo_accession
    GSMid = Meta(GSEDataobj)$sample_id
    platform_id = Meta(GSEDataobj)$platform_id
    if( ! class(GEOdata_matrix) == 'try-error' ){
      sAnnot = pData(phenoData(GEOdata_matrix[[1]]))[, 2:1]
    }else{
      sAnnot = t(sapply(Meta(GSEDataobj)$sample_id, function(x)  c(x, Meta( getGEO(x, destdir = study_files_directory))$title) ))
      sAnnot = as.data.frame(sAnnot)
    }
  }
  if( length(grep('GDS', my_GEOid)) == 1 ){
    GDSDataobj = try(getGEO(my_GEOid, GSEMatrix = F, destdir = study_files_directory), silent = T)
    GEOdata_matrix = try(getGEO(my_GEOid, GSEMatrix = T, destdir = study_files_directory), silent = T)
    GSEid = Meta(GDSDataobj)$reference_series
    platform_id = Meta(GDSDataobj)$platform
    sAnnot = Columns(GDSDataobj)
  }
  # Define samples 
  names(sAnnot)[1] = 'sample'
  names(sAnnot)[2] = 'phenotype'
  sAnnot$phenotype = format_phenoname(pheno_name = sAnnot$phenotype)
  samples = as.character(sAnnot$sample[ sAnnot$phenotype %in% names(which(table(sAnnot$phenotype)>1)) ])
  rownames(sAnnot) = sAnnot$sample
  try(getGEOSuppFiles(GEO = GSEid, baseDir = study_files_directory), silent = T)
  cel_tarfile = list.files(study_files_directory, pattern = '_RAW.tar', full.names = T, ignore.case = T, recursive = T)
  eset = NULL
  has_cells = T
  # Load cel files and normalize
  if( ! class(cel_tarfile) == 'try-error' & length(cel_tarfile) > 0 ){
    untar(cel_tarfile[1], exdir = study_files_directory, extras = '-vz')
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
    has_cells = F
    eset = try(GDS2eSet(GEOdata_matrix,do.log2=T), silent = T)
    if( class(eset) == 'try-error' )
      return(NULL)
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
  return(list(eset=eset, sAnnot=sAnnot, study_files_directory=study_files_directory, has_cells=has_cells))
}


get_experiment_design = function(eset, sAnnot, study_id){
  class_column = EXPERIMENTAL_DESIGNS[[study_id]]$class_column
  pos_sample = EXPERIMENTAL_DESIGNS[[study_id]]$pos_sample
  neg_sample = EXPERIMENTAL_DESIGNS[[study_id]]$neg_sample
  covariates = EXPERIMENTAL_DESIGNS[[study_id]]$covariates
  
  sAnnot = sAnnot[ sAnnot[,class_column] %in% c(pos_sample, neg_sample), ]
  sAnnot[,class_column] = (sAnnot[,class_column] %in% pos_sample) + 0
  sAnnot = sAnnot[ order(sAnnot[,class_column]) , ]
  
  phenotype = sAnnot[,class_column]
  design = model.matrix(~phenotype)
  for( i in covariates ){
    cov_ =  sAnnot[,i]
    design = cbind(design, model.matrix(~cov_))
  }
  design =  design[ , ! duplicated(colnames(design)) ]
  rownames(design) = sAnnot$sample
  return(list(design=design, my_contrasts=list(pos_sample=pos_sample, neg_sample=neg_sample)))
}



run_differential_expression = function(eset, sAnnot, study_id){
  if( is.null(EXPERIMENTAL_DESIGNS[[study_id]]) )
    return(NULL)
  # Design phenotipic contrasts
  my_design = get_experiment_design(eset, sAnnot, study_id)
  design = my_design$design
  design = design[ rownames(design) %in% sampleNames(eset), ]
  # Fit model with contrasts
  fit = lmFit(eset[, sampleNames(eset) %in% rownames(design) ], design)
  fit = eBayes(fit)
  dea = topTable(fit, coef="phenotype", adjust="fdr", number = nrow(eset))
  return(list(dea=dea, contrast=my_design$my_contrasts))
}

EXPERIMENTAL_DESIGNS = list("myc_17159920_cancer_cell_lines_lof_human_gpl570_gds2526" = list(class_column = 'phenotype', 
                                                                                             pos_sample = 'CONTROL',
                                                                                             neg_sample = 'C_MYC_KNOCKDOWN',
                                                                                             covariates = c('other','cell.line')),
                            
                            "ar_21330406_lncap_lof_human_gpl570_gds4113" = list(class_column = 'phenotype', 
                                                                                pos_sample = 'CONTROL_SHRNA_+_VEHICLE',
                                                                                neg_sample = 'ANDROGEN_RECEPTOR_SHRNA_+_VEHICLE'),
                            
                            "elk1_23426362_lncap_lof_human_gpl570_gse34589" = list(class_column = 'phenotype', 
                                                                                   pos_sample = 'LNCAP_+_VEHICLE',
                                                                                   neg_sample = 'LNCAP_+_ELK1_SHRNA'),
                            
                            "znf217_22593193_mda_mb_231_gof_human_gpl570_gse35511" = list(class_column = 'phenotype', 
                                                                                          pos_sample = 'MDA_MB_231_PCDNA6_CELL_CULTURE', 
                                                                                          neg_sample = 'MDA_MB_231_ZNF217_CELL_CULTURE'),
                            
                            'bcor_22012066_cn_aml_lof_human_gpl570_gds4280' = list(class_column = 'phenotype', 
                                                                                   pos_sample = 'WILD_TYPE', 
                                                                                   neg_sample = 'BCOR_MUTATION'),
                            
                            'ccnd1_18413728_imr_neuroblastoma_lof_human_gpl570_gse8866' = list(class_column = 'phenotype', 
                                                                                               pos_sample = 'IMR_32_GFP_SIRNA_T_48', 
                                                                                               neg_sample = c('IMR_32_CCND1_SIRNA_T_48', 'IMR_32_CCND1SIRNA_T_48')),
                            
                            'cebpa_22442349_pbmc_lof_human_gpl570_gds4407' = list(class_column = 'phenotype', 
                                                                                  pos_sample = 'CONTROL', 
                                                                                  neg_sample = 'BIALLELIC_AML_MLD_POSITIVE'),
                            
                            'creb1_18801183_k562_lof_human_gpl570_gds3487' = list(class_column = 'phenotype', 
                                                                                  pos_sample = 'CONTROL', 
                                                                                  neg_sample = 'CREB_DEPLETION'),
                            
                            'ctnnb1_19652203_myeloma_lof_human_gpl570_gds3578' = list(class_column = 'phenotype', 
                                                                                      pos_sample = 'CONTROL',
                                                                                      neg_sample = 'BETA_CATENIN_DEPLETION'),
                            
                            'ctnnb1_25082960_aspc1_lof_human_gpl570_gds4386' = list(class_column = 'phenotype', 
                                                                                    pos_sample = 'CONTROL',
                                                                                    neg_sample = 'BETA_CATENIN_DEPLETION'),
                            
                            'dlx4_00000000_mcfdash7_gof_human_gpl201_gse21657' = list(class_column = 'phenotype', 
                                                                                      pos_sample = 'MCF_TRANSFECTED_WITH_PCDNA3.2_BP1_.',
                                                                                      neg_sample = 'MCF_TRANSFECTED_WITH_PCDNA3.2_.'),
                            
                            'elk1_23426362_lncap_lof_human_gpl570_gse34589' = list(class_column = 'phenotype', 
                                                                                   pos_sample = 'LNCAP_+_VEHICLE',
                                                                                   neg_sample = 'LNCAP_+_ELK1_SHRNA'),
                            
                            'erg_19359602_huvec_lof_human_gpl570_gds3557' = list(class_column = 'phenotype', 
                                                                                 pos_sample = 'CONTROL_VECTOR', 
                                                                                 neg_sample = 'ERG_SIRNA'),
                            
                            'esr1_21713035_mcf7_lof_human_gpl570_gds4061' = list(class_column = 'phenotype', 
                                                                                 pos_sample = 'CONTROL',
                                                                                 neg_sample = 'ESTROGEN_RECEPTOR_KNOCKDOWN'),
                            
                            'gata3_21892208_mda_mb_231_gof_human_gpl570_gds4080' = list(class_column = 'phenotype', 
                                                                                        pos_sample = 'GATA',
                                                                                        neg_sample = 'CONTROL'),
                            
                            'hif1a_16565084_mcf7_lof_human_gpl2507_gds2761' = list(class_column = 'phenotype', 
                                                                                   pos_sample = 'CONTROL',
                                                                                   neg_sample = 'HIF_1ALPHA_DEPLETION'),
                            
                            'hif2a_16565084_mcf7_lof_human_gpl2507_gds2761' = list(class_column = 'phenotype', 
                                                                                   pos_sample = 'CONTROL',
                                                                                   neg_sample = 'HIF_2ALPHA_DEPLETION'),
                            
                            
                            'hivep2_22294689_hek293_gof_human_gpl6244_gds5213' = list(class_column = 'phenotype', 
                                                                                      pos_sample = 'MIBP1_OVEREXPRESSION',
                                                                                      neg_sample = 'CONTROL'),
                            
                            
                            'hnf4a_21852396_wao9_lof_human_gpl570_gds3926' = list(class_column = 'phenotype', 
                                                                                  pos_sample = 'CONTROL',
                                                                                  neg_sample = 'HNF4Α_KNOCKDOWN'),
                            
                            'hsf1_17216044_hela_lof_human_gpl571_gds1733' = list(class_column = 'phenotype', 
                                                                                 pos_sample = 'CONTROL',
                                                                                 neg_sample = 'HSF1_KNOCKDOWN',
                                                                                 covariates = c('time')),
                            
                            'klf4_17017123_rko_gof_human_gpl96_gds1942' = list(class_column = 'phenotype',
                                                                               pos_sample = 'INDUCED',
                                                                               neg_sample = 'CONTROL',
                                                                               covariates = c('time')),
                            
                            'lmo4_17452977_mcf7_lof_human_gpl570_gds2789' = list(class_column = 'phenotype', 
                                                                                 pos_sample = 'BASELINE',
                                                                                 neg_sample = 'DOMINANT_NEGATIVE_CLIM_INDUCTION'),
                            
                            'ncoa_22072566_ly2_lof_human_gpl570_gds4095' = list(class_column = 'phenotype', 
                                                                                pos_sample = 'CONTROL',
                                                                                neg_sample = 'SRC1_KNOCKDOWN',
                                                                                covariates = c('agent')),
                            
                            'nod2_21335489_hek293_lof_human_gpl570_gds4416' = list(class_column = 'phenotype', 
                                                                                   pos_sample = 'CONTROL',
                                                                                   neg_sample = 'NOD2_WILD_TYPE',
                                                                                   covariates = c('agent', 'time')),
                            
                            'bmi1_17452456_medulloblastoma_lof_human_gpl570_gds2724' = list(class_column = 'phenotype', 
                                                                                            pos_sample = 'CONTROL',
                                                                                            neg_sample = 'BMI_KNOCKDOWN'),
                            
                            'esr1_21299862_mcf7_lof_human__gds4065' = list(class_column = 'phenotype', 
                                                                           pos_sample = 'CONTROL',
                                                                           neg_sample = 'ESTROGEN_RECEPTOR_KNOCKDOWN',
                                                                           covariates = c('agent')),
                            
                            'rara_21299862_mcf7_lof_human__gds4065' = list(class_column = 'phenotype', 
                                                                           pos_sample = 'CONTROL',
                                                                           neg_sample = 'RETINOIC_ACID_RECEPTOR_KNOCKDOWN',
                                                                           covariates = c('agent')),
                            
                            'rent1_15448691_hela_lof_human_gpl8300_gds705' = list(class_column = 'phenotype', 
                                                                                  pos_sample = 'CONTROL',
                                                                                  neg_sample = 'RENT1_KNOCKDOWN'),
                            
                            'sin3a_22783022_mcf7_lof_human_gpl570_gds4388' = list(class_column = 'phenotype', 
                                                                                  pos_sample = 'CONTROL',
                                                                                  neg_sample = 'SIN3A_KNOCKDOWN',
                                                                                  covariates = c('agent')),
                            
                            'rnf2_20805357_u2os_osteosarcoma_lof_human_gpl570_gse23035' = list(class_column = 'phenotype', 
                                                                                               pos_sample = 'SHCONTROL_CELLS',
                                                                                               neg_sample = 'SHBAP1_#2_CELLS'),
                            
                            'sox4_16636670_acc3_lof_human_gpl96_gds2193' = list(class_column = 'phenotype', 
                                                                                pos_sample = 'CONTROL',
                                                                                neg_sample = 'SOX4_RNAI'),
                            
                            'sox7_18682240_esc_gof_human_gpl570_gds3300' = list(class_column = 'protocol', 
                                                                                pos_sample = 'SOX7 overexpression',
                                                                                neg_sample = 'control'),
                            
                            'sox17_18682240_esc_gof_human_gpl570_gds3300' = list(class_column = 'protocol', 
                                                                                 pos_sample = 'SOX17 overexpression',
                                                                                 neg_sample = 'control'),
                            
                            
                            'stat3_00000000_a549_lof_human_gpl571_gse42979' = list(class_column = 'phenotype', 
                                                                                   pos_sample = 'A549_EGFPLUS',
                                                                                   neg_sample = 'A549_SIRNA_STAT3'),
                            
                            
                            'tardbp_19910924_hek293e_lof_human_gpl570_gds3730' = list(class_column = 'phenotype', 
                                                                                      pos_sample = 'CONTROL',
                                                                                      neg_sample = 'TDP_43_KNOCKDOWN'),
                            
                            'tp63_16885358_squamous_epithelia_lof_human_gpl570_gds2088' = list(class_column = 'protocol', 
                                                                                               pos_sample = 'pos_sample',
                                                                                               neg_sample = 'p63 knockdown',
                                                                                               covariates = c('phenotype', 'cell.type')),
                            
                            'wt1_17088532_huvec_lof_human_gpl570_gds2010' = list(class_column = 'phenotype', 
                                                                                 pos_sample = 'CONTROL',
                                                                                 neg_sample = 'WTAP_KNOCKDOWN'),
                            
                            'yap1_18413746_mcf10a_gof_human_gpl570_gds3220' = list(class_column = 'phenotype', 
                                                                                   pos_sample = 'YAP_OVEREXPRESSION',
                                                                                   neg_sample = 'CONTROL'),
                            
                            'yy1_20215434_hela_lof_human_gpl570_gds3788' = list(class_column = 'phenotype', 
                                                                                pos_sample = 'CONTROL',
                                                                                neg_sample = 'SHYY1_KNOCK_DOWN'),
                            
                            'znf148_21828133_erythroblast_lof_human_gpl571_gse31092' = list(class_column = 'phenotype', 
                                                                                            pos_sample = 'EMPTY_VECTOR_BIOLOGICAL',
                                                                                            neg_sample = 'SHZBP_89_BIOLOGICAL'),
                            
                            'znf217_22593193_mda_mb_231_gof_human_gpl570_gse35511' = list(class_column = 'phenotype', 
                                                                                          pos_sample = 'MDA_MB_231_ZNF217_CELL_CULTURE',
                                                                                          neg_sample = 'MDA_MB_231_PCDNA6_CELL_CULTURE'),
                            
                            'ctnnb1_21914722_ls174t_lof_human_gpl570_gds4386' = list(class_column = 'protocol', 
                                                                                     pos_sample = 'β-catenin shRNA, uninduced',
                                                                                     neg_sample = 'β-catenin shRNA, induced'),
                            
                            'esr1_22536322_mcf7_lof_human_gpl570_gse31912' = list(class_column = 'phenotype', 
                                                                                  pos_sample = 'MCF7_SIRNA_SICONTROL48H',
                                                                                  neg_sample = 'MCF7_SIRNA_ESR148H'),
                            
                            'nfkb1_22536322_mcf7_lof_human_gpl570_gse31912' = list(class_column = 'phenotype', 
                                                                                   pos_sample = 'MCF7_SIRNA_SICONTROL48H',
                                                                                   neg_sample = 'MCF7_SIRNA_NFKB148H'),
                            
                            'pax3_22536322_mcf7_lof_human_gpl570_gse31912' = list(class_column = 'phenotype', 
                                                                                  pos_sample = 'MCF7_SIRNA_SICONTROL48H',
                                                                                  neg_sample = 'MCF7_SIRNA_PAX348H')
                            
)
