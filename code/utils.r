library(plyr)
library(reshape2)
library(ggplot2)
source("code/lib_basics.r")
source("code/my_ggplot_theme.r")


############################################################################################################
## Load data
############################################################################################################

# ensemblgenes_annot =  get(load(file = '/Volumes/GoogleDrive/My Drive/databases/ensemblgenes_annot.Rdata'))

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
  # write.table(TF_census, file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_census/vaquerizas/TF_census.txt', col.names = F, row.names = F, quote = F)
  # TF_census = read.delim(file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_census/vaquerizas/TF_census.txt', header = F, stringsAsFactors = F)[,1]
  TF_census = unique(read.delim(file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_classification.txt', header = F, stringsAsFactors = F)[,2])
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
  complexes = read.delim(file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_info/regulation_type/uniprot/complex_classes.txt', header = T, stringsAsFactors = F)
}


load_CHregulation = function(){
  CHregulation = read.delim(file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_info/epigenetics/Ehsani_2016/chromatin_modulator.txt', header = T, stringsAsFactors = F)
  CHregulation$Chromatin_Opening_Type[ CHregulation$Chromatin_Opening_Type == ''] = "unknown"
  return(CHregulation)
}


load_TFclass_classes = function(){
  TFclassification = unique(read.delim(file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_info/census/TFClass/huTF_classification.txt', header = T, stringsAsFactors = F))
  return(TFclassification)
}



load_TFtissues = function(){
  TF_x_tissue = unique(read.delim(file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_info/expression/GTEx_TF_per_tissue.txt', header = F, stringsAsFactors = F))
  return(TF_x_tissue)
}





## benchmark
load_Amp = function(){
  amp = (get(load(file = '/Volumes/GoogleDrive/My Drive/datasets/jsr-gdsc/CopyNumerVariation/Gene_level_CN_max.rdata')) >= 8) + 0
  fpkms = get(load('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/rnaseq/data/expression/processed/rnaseq_cell_lines/fpkms_aggregatedduplicates.RData'))
  load('/Volumes/GoogleDrive/My Drive/databases/ensemblgenes_annot.Rdata')
  ensemblgenes_annot = subset(ensemblgenes_annot, hgnc_symbol != '')
  rownames(fpkms) = ensemblgenes_annot$hgnc_symbol[ match(rownames(fpkms), ensemblgenes_annot$ensembl_gene_id) ]
  genes = intersect(rownames(amp), rownames(fpkms))
  samples = intersect(colnames(amp), colnames(fpkms))
  fpkms = fpkms[ genes, samples]
  amp = amp[ genes, samples ]
  amp[ fpkms < 2 ] = 0
  return(amp)
}


load_homDel = function(){
  del = (get(load(file = '/Volumes/GoogleDrive/My Drive/datasets/jsr-gdsc/CopyNumerVariation/Gene_level_CN_max.rdata')) == 0 ) + 0
  fpkms = get(load('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/rnaseq/data/expression/processed/rnaseq_cell_lines/fpkms_aggregatedduplicates.RData'))
  load('/Volumes/GoogleDrive/My Drive/databases/ensemblgenes_annot.Rdata')
  ensemblgenes_annot = subset(ensemblgenes_annot, hgnc_symbol != '')
  rownames(fpkms) = ensemblgenes_annot$hgnc_symbol[ match(rownames(fpkms), ensemblgenes_annot$ensembl_gene_id) ]
  genes = intersect(rownames(del), rownames(fpkms))
  samples = intersect(colnames(del), colnames(fpkms))
  fpkms = fpkms[ genes, samples]
  del = del[ genes, samples ]
  del[ fpkms > 2 ] = 0
  return(del)
}


load_ACHILLES = function(){
  essentiality = get(load('/Volumes/GoogleDrive/My Drive/datasets/achilles/v2.20.2/Achilles_v2.20.2_GeneSolutions.clean.Rdata'))
}

load_DRIVEproject = function(){
  library(stringr)
  essentiality = readRDS('data/regulons_QC/cell_lines/essentiality/DRIVE_ATARiS_data.rds')
  # manage names
  colnames(essentiality) = sapply(strsplit(colnames(essentiality), split = '_'), head, 1)
  celllines_annotation = get(load('../../datasets/jsr-gdsc/R_objects/cell_lines/COSMIC_covariates_19062013_updated_luz.RData'))
  celllines_annotation$Cell.line.name_formated = str_replace_all(celllines_annotation$Cell.line.name, "[[:punct:]]", " ") %>% gsub(' ', '', .) %>% tolower(.)
  essentiality = essentiality[ , colnames(essentiality) %in%  celllines_annotation$Cell.line.name_formated ]
  colnames(essentiality) = celllines_annotation$CosmicID[ match(colnames(essentiality), celllines_annotation$Cell.line.name_formated) ]
  # z-transform
  essentiality = as.matrix(essentiality)
  return(essentiality)
}


average_duplicated_genes = function(X){
  genes = setdiff(rownames(X), '')
  mat = matrix(NA, ncol = ncol(X), nrow = length(genes), dimnames = list(genes, colnames(X)) )
  for( g in genes ){
    if( sum(rownames(X) == g) > 1  ){
      mat[g, ] = apply( X[ rownames(X) == g,  ], 2, mean)
    }else{
      mat[g, ] = X[ rownames(X) == g,  ]
    }
  }
  return(mat)
}
############################################################################################################








############################################################################################################
## plots
############################################################################################################

summarize_network = function(network){
  network$TF = as.character(network$TF)
  network$target = as.character(network$target)
  network = unique(network)
  
  TF_freq = as.numeric(unname(table(network$TF)))
  targets_freq = as.numeric(unname(table(network$target)))
  message('Number of edges: ', nrow(network))
  message('Number of TFs: ', length(unique(network$TF)))
  message('Number of targets: ', length(unique(network$target)))
  message('Targets per TF')
  print( summary(TF_freq) )
  message('Top TFs')
  print(head(sort(table(network$TF), decreasing = T)))
  message('TFs per targets: ')
  print( summary(targets_freq) )
  message('Top targets')
  print(head(sort(table(network$target), decreasing = T)))
  
  
  plot_df = as.data.frame(rbind(cbind(freq=TF_freq, gene_class='TF'), 
                                cbind(freq=targets_freq, gene_class='targets')), stringsAsFactors = F)
  plot_df$freq = as.numeric(plot_df$freq)
  ggplot(plot_df, aes(x = gene_class, y = freq) ) + geom_boxplot() + facet_wrap(~gene_class, scales = 'free') + 
    theme_light(18) + theme(axis.title.x = element_blank()) + ylab('Freq')
}

############################################################################################################







############################################################################################################
## Other
############################################################################################################

network_similarity =function(network){
  network$TF = as.character(network$TF)
  network$target = as.character(network$target)
  network = unique(network)
  
  tfs = unique(network$TF)
  regulons = lapply(tfs, function(tf) subset(network, TF == tf)$target  )
  names(regulons) = tfs
  regulons = regulons[ sapply(regulons, length) > 2  ]
  overlap = sapply( names(regulons) , function(tf1) sapply( names(regulons) , function(tf1, tf2){
    length(intersect(regulons[[tf1]], regulons[[tf2]])) / length(unique(c(regulons[[tf1]], regulons[[tf2]])))
  } , tf1 ) )
  pheatmap::pheatmap(overlap)
  diag(overlap) = NA
  print(summary(c(overlap)))
  print('Top overlaps')
  df = melt(overlap)
  print( head(df[ order(df$value, decreasing = T), ], 10) )
}


effectsize_cohensD  =  function(x, y) { # Effect size
  lx  =  length(x)- 1
  ly  =  length(y)- 1
  md   =  abs(mean(x) - mean(y)) ## mean difference (numerator)
  csd  =  lx * var(x) + ly * var(y)
  csd  =  csd/(lx + ly)
  csd  =  sqrt(csd)                     ## common sd computation
  cd   =  md/csd                        ## cohen's d
}


manage_datasets_names = function(x, what = 'dataset2group'){
  # message('Options for what are:')
  # message('- regulon_id2TF')
  # message('- regulon_id2regulon_dataset_full')
  # message('- regulon_dataset_full2evidence')
  # message('- regulon_dataset_full2dataset')
  # message('- dataset2group')
  x_out = x
  if(what == 'regulon_id2TF'){
    x_out = sapply(strsplit(x, ' '), head, 1)
    if( length(grep('_[A-E]$', x_out)) > 0 )
      x_out[ grep('_[A-E]$', x_out) ] = sapply(strsplit(x_out[ grep('_[A-E]$', x_out) ], '_'), head, 1)
  } 
  if(what == 'regulon_id2regulon_dataset_full'){
    x_out = sapply(strsplit(x, ' '), tail, 1) %>% gsub('.sif', '', .)
  }
  if(what == 'regulon_dataset_full2evidence'){
    x_out = sapply(strsplit(x, '/'), head, 1)
  }
  if(what == 'regulon_dataset_full2dataset'){
    x_out = sapply(strsplit(x, '/'), function(xx) paste(xx[2:3], collapse = '_') ) %>%
      gsub('network', '', .) %>%
      gsub('_NA', '', .) %>%   gsub('__', '_', .) %>% 
      gsub('_$', '', .)  %>% gsub('^_', 'inferredGTEx_', .) %>% 
      gsub('scanning_output_', '', .) %>% 
      gsub('_signed_signed', '_signed', .) 
  }
  if(what == 'dataset2group'){
    x_out[ grep('^inferredGTEx_t', x_out) ] = 'tissue-specific'
    x_out[ grep('^inferredGTEX_t', x_out) ] = 'tissue-specific'
    x_out[ grep('^inferredGTEX ', x_out) ] = 'tissue-specific'
    x_out[ grep('^inferredGTEx_', x_out) ] = 'tissues consensus'
    x_out[ grep('hocomoco', x_out) ] = 'hocomoco'
    x_out[ grep('ReMap', x_out) ] = 'ReMap'
    x_out[ grep('jaspar', x_out) ] = 'jaspar'
    # x_out[ x_out %in% c(LETTERS[1:5], 'TOP')  ] = 'omnipath_scores'
    x_out = gsub('_signed_signed', '', x_out)
    x_out = gsub('evidence2_', '', x_out)
  }
  if(what == 'evidence2color'){
    x_out = c('coral', my_color_palette$EMBL[c(5,5,8,3,4,2,6,1,7)])[
      match(x_out,
            c('chipSeq', 'consensus', 'consensus_between', 'consensus_within', 'curated_databases', 'inferredGTEx',  'old_consensus', 'TFBS_scanning', 'omnipath_scores') )]
  }
  if(what == 'evidence2shape'){
    x_out = c(15, 6, 6, 6, 16:17, 5, 18, 16, 1)[
      match(x_out,
            c('chipSeq', 'consensus', 'consensus_between', 'consensus_within', 'curated_databases', 'inferredGTEx',  'old_consensus', 'TFBS_scanning', 'omnipath_scores') )]
  }
  return(x_out)
}


get_percentile = function(ob_nes, ra_nes){
  per = ecdf(ra_nes)(ob_nes)
}
