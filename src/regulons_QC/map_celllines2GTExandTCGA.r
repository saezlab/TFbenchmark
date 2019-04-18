#' Annotate GDSC cancer cell lines to their corresponging TCGA-GTEx tissues
#' by manual curation



rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)



celllines_annotation = get(load('data/regulons_QC/B3_cell_lines/cosmic_cancercell_lines_information_LGA2018CancerResearch.rdata'))
head(celllines_annotation)
table(celllines_annotation$Study.Abbreviation)
tissues_map = list('ALL'='Blood',
                  'CLL'='Blood',
                  'DLBC'='Blood',
                  'LAML'='Blood',
                  'LCML'='Blood',
                  'BLCA'='Blader',
                  'LGG'='Brain',
                  'GBM'='Brain',
                  'MB'='Brain',
                  'BRCA'='Breast',
                  'UCEC'='Uterus',
                  'COREAD'='Colon',
                  'ESCA'='Esophagus',
                  'KIRC'='Kidney',
                  'LIHC'='Liver',
                  'LUAD'='Lung',
                  'LUSC'='Lung',
                  'SCLC'='Lung',
                  'OV'='Ovary',
                  'PAAD'='Pancreas',
                  'PRAD'='Prostate',
                  'SKCM'='Skin',
                  'STAD'='Stomach',
                  'THCA'='Thyroid')
celllines_annotation$GTEx = 'none'
celllines_annotation$GTEx[ celllines_annotation$Study.Abbreviation %in% names(tissues_map) ] = unlist(tissues_map[ 
  celllines_annotation$Study.Abbreviation[ celllines_annotation$Study.Abbreviation %in% names(tissues_map) ]])
table(celllines_annotation$GTEx)
unique(celllines_annotation[ , c('GTEx', 'Study.Abbreviation') ])

save(celllines_annotation, file = 'data/regulons_QC/B3_cell_lines/cellines2tissues_mapping.rdata')
