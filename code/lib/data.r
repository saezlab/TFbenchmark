
############################################################################################################
## Load data
############################################################################################################


load_proteins_annotation = function(type){
  if( type == 'ensembl' ){
    # library(biomaRt)
    # message('Loading Ensembl-Biomart data')
    # ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
    # ensemblgenes_annot = getBM(mart = ensembl,
    #                            attributes=c('hgnc_symbol', 'uniprotswissprot', 'ensembl_gene_id'),
    #                            filters = 'transcript_biotype', values = c('protein_coding') )
    # annot = unique(ensemblgenes_annot)
    # save(annot, file = 'data/annotations/hsa_ensembl_proteincoding_20180314.rdata')
    # write.table(annot, file = 'data/annotations/hsa_ensembl_proteincoding_20180314.txt', sep = '\t', col.names = T, row.names = F, quote = F)
    load(file = 'data/annotations/hsa_ensembl_proteincoding_20180314.rdata')
  }
  if(type == 'uniprot'){
    annot = read.delim('data/annotations/hsa_uniprot_20180314.txt', stringsAsFactors = F)
    annot$Gene.name.primary = ''
    annot$Gene.name.primary[ annot$Gene.names != '' ] = sapply(strsplit(annot$Gene.names[ annot$Gene.names != '' ], split = ' '), head, 1)
    annot$Gene.name.primary[ annot$Gene.names != '' ] = sapply(strsplit(annot$Gene.name.primary[ annot$Gene.names != '' ], split = ';'), head, 1)
    annot$Gene.names.secondary = sapply(strsplit(annot$Gene.names, split = ' '), function(x) paste(x[-1], collapse = ';') )
  }
  return(annot)
}




## TF properties
load_TFs_census = function(){
  # message('Load TF census')
  # TF_census_vaquerizas = read.delim('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/DoRothEA/DATA/regulons/census_vaquerizas2009/nrg2538-s3_clean.txt', stringsAsFactors = F, header = F)
  # TF_census = unique(TF_census_vaquerizas$V6[TF_census_vaquerizas$V1 %in% c('a', 'b') ])
  # TF_census = sort(setdiff(TF_census, ''))
  # write.table(TF_census, file = 'data/TF_census/vaquerizas/TF_census.txt', col.names = F, row.names = F, quote = F)
  # TF_census = read.delim(file = 'data/TF_census/vaquerizas/TF_census.txt', header = F, stringsAsFactors = F)[,1]
  TF_census = unique(read.delim(file = 'data/TF_census/TFClass/huTF_census.txt', header = F, stringsAsFactors = F)[,1])
}


load_activators = function(){
  act_rep = read.csv('data/TF_info/regulation_type/uniprot/activators_repressors_curated.csv', stringsAsFactors = F)
  activators = unique(subset(act_rep, Veredict_curated == 'A')$Gene.names)
  return(activators)
}

load_dualActRep = function(){
  act_rep = read.csv('data/TF_info/regulation_type/uniprot/activators_repressors_curated.csv', stringsAsFactors = F)
  dual = unique(subset(act_rep, Veredict_curated == 'A,R')$Gene.names)
  return(dual)
}


load_repressors = function(strict = F){
  stric_repressors = c('REST', 'ZEB1', 'ZEB2', 'E2F6', 'NFKB2', 'ZHX2', 'NCOR2', 'IRF2', 'BCL6')
  if( ! strict ){
    act_rep = read.csv('data/TF_info/regulation_type/uniprot/activators_repressors_curated.csv', stringsAsFactors = F)
    repressors = unique(subset(act_rep, Veredict_curated == 'R')$Gene.names)
    return(repressors)
  }else{
    return(stric_repressors)
  }
}

load_complexes = function(){
  complexes = read.delim(file = 'data/TF_info/regulation_type/uniprot/complex_classes.txt', header = T, stringsAsFactors = F)
}


load_CHregulation = function(){
  CHregulation = read.delim(file = 'data/TF_info/epigenetics/Ehsani_2016/chromatin_modulator.txt', header = T, stringsAsFactors = F)
  CHregulation$Chromatin_Opening_Type[ CHregulation$Chromatin_Opening_Type == ''] = "unknown"
  return(CHregulation)
}


load_TFclass_classes = function(){
  TFclassification = unique(read.delim(file = 'data/TF_info/census/TFClass/huTF_classification.txt', header = T, stringsAsFactors = F))
  return(TFclassification)
}



load_TFtissues = function(){
  TF_x_tissue = unique(read.delim(file = 'data/TF_info/expression/GTEx_TF_per_tissue.txt', header = F, stringsAsFactors = F))
  return(TF_x_tissue)
}



## benchmark
load_Amp = function(){
  # amp = (get(load(file = '/Volumes/GoogleDrive/My Drive/datasets/jsr-gdsc/CopyNumerVariation/Gene_level_CN_max.rdata')) >= 8) + 0
  # fpkms = get(load('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/rnaseq/data/expression/processed/rnaseq_cell_lines/fpkms_aggregatedduplicates.RData'))
  # load('/Volumes/GoogleDrive/My Drive/databases/ensemblgenes_annot.Rdata')
  # ensemblgenes_annot = subset(ensemblgenes_annot, hgnc_symbol != '')
  # rownames(fpkms) = ensemblgenes_annot$hgnc_symbol[ match(rownames(fpkms), ensemblgenes_annot$ensembl_gene_id) ]
  # genes = intersect(rownames(amp), rownames(fpkms))
  # samples = intersect(colnames(amp), colnames(fpkms))
  # fpkms = fpkms[ genes, samples]
  # amp = amp[ genes, samples ]
  # amp[ fpkms < 2 ] = 0
  # save(amp, file = 'data/regulons_QC/B3_cell_lines/cna/amplifications.rdata')
  load('data/regulons_QC/B3_cell_lines/cna/amplifications.rdata')
  return(amp)
}


load_homDel = function(){
  # del = (get(load(file = '/Volumes/GoogleDrive/My Drive/datasets/jsr-gdsc/CopyNumerVariation/Gene_level_CN_max.rdata')) == 0 ) + 0
  # fpkms = get(load('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/rnaseq/data/expression/processed/rnaseq_cell_lines/fpkms_aggregatedduplicates.RData'))
  # load('/Volumes/GoogleDrive/My Drive/databases/ensemblgenes_annot.Rdata')
  # ensemblgenes_annot = subset(ensemblgenes_annot, hgnc_symbol != '')
  # rownames(fpkms) = ensemblgenes_annot$hgnc_symbol[ match(rownames(fpkms), ensemblgenes_annot$ensembl_gene_id) ]
  # genes = intersect(rownames(del), rownames(fpkms))
  # samples = intersect(colnames(del), colnames(fpkms))
  # fpkms = fpkms[ genes, samples]
  # del = del[ genes, samples ]
  # del[ fpkms > 2 ] = 0
  # save(del, file = 'data/regulons_QC/B3_cell_lines/cna/deletions.rdata')
  load('data/regulons_QC/B3_cell_lines/cna/deletions.rdata')
  return(del)
}


load_ACHILLES = function(){
  load('data/regulons_QC/B3_cell_lines/essentiality/Achilles_v2.20.2_GeneSolutions.clean.Rdata')
  return(essentiality)
}

load_DRIVEproject = function(){
  # require(stringr)
  # essentiality = readRDS('data/regulons_QC/B3_cell_lines/essentiality/DRIVE_ATARiS_data.rds')
  # # manage names
  # colnames(essentiality) = sapply(strsplit(colnames(essentiality), split = '_'), head, 1)
  # celllines_annotation = get(load('../../datasets/jsr-gdsc/R_objects/cell_lines/COSMIC_covariates_19062013_updated_luz.RData'))
  # celllines_annotation$Cell.line.name_formated = str_replace_all(celllines_annotation$Cell.line.name, "[[:punct:]]", " ") %>% gsub(' ', '', .) %>% tolower(.)
  # essentiality = essentiality[ , colnames(essentiality) %in%  celllines_annotation$Cell.line.name_formated ]
  # colnames(essentiality) = celllines_annotation$CosmicID[ match(colnames(essentiality), celllines_annotation$Cell.line.name_formated) ]
  # # z-transform
  # essentiality = as.matrix(essentiality)
  # save(essentiality, file = 'data/regulons_QC/B3_cell_lines/essentiality/DRIVE_ATARiS_data_clean.rdata')
  load('data/regulons_QC/B3_cell_lines/essentiality/DRIVE_ATARiS_data_clean.rdata')
  return(essentiality)
}





