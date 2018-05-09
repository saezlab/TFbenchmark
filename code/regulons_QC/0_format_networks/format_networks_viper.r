rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')



# Set minimum regulon size
Nmax = 10000
Nmin = 4
remove_overrepresented_targets = T
overrepresented_targets_th = 0.1


# Important! considere sing for repressors!
repressors = load_repressors(strict = F)
TFcensus = load_TFs_census()



### --------------------------------------------------------------------------------------------------------------
# Load REAL networks + format for viper input
network_files = list.files('data/TF_target_sources/', recursive = T, pattern = 'network') %>% grep('sif$', ., value=T)  %>% 
  grep('RANDOM', ., value=T, invert= T) %>% grep('omni', ., value=T) 
for (n in network_files){
  viper_regulon = list()
  message(n)
  # Load network and fromat regulons
  net = read.delim(paste('data/TF_target_sources/', n, sep = ''), stringsAsFactors = F, sep = '\t')
  net = unique(net)
  # Remove self regulation
  net = subset(net, TF != target)
  tfs = unique(net$TF)
  message('TFs: ', length(tfs))
  if( length(grep('_', net$TF)) > 0 ){
    net$regulon_id = net$TF
    net$TF = sapply(strsplit(net$TF, split = '_'), head, 1)
  }
  # Remove empty targets
  net = subset(net, target != "")
  # Remove overrepresented targets
  if( remove_overrepresented_targets & length(tfs) > 20 ){
    overrepresented_targets = which( table(net$target) / length(tfs) > overrepresented_targets_th) %>% names(.)
    net = subset(net, ! target %in% overrepresented_targets)
  }
  # Define TFs of interest
  tfs_of_interest = names(which(table(net$TF) >= Nmin))
  message('TFs with (n>=', Nmin,'): ', length(tfs_of_interest))
  if( length(tfs_of_interest) < 2 )
    next()
  for (tf in tfs_of_interest ){
    tf_targets = subset(net, TF == tf)$target
    # Define TF mode
    if( ! is.null(net$effect) ){ # if known, use predefined
      tf_targets_sign = subset(net, TF == tf)$effect
      if( tf %in% repressors ) # unless effect equals 0, check if tf is repressor
        tf_targets_sign[ tf_targets_sign == 0 ] = - 1
      if( ! tf %in% repressors ) # unless effect equals 0, check if tf is repressor
        tf_targets_sign[ tf_targets_sign == 0  ] = 1
    }else{ # if unknown, check if tf is repressor
      tf_targets_sign = rep(1, length(tf_targets))
      if( tf %in% repressors )
        tf_targets_sign = -tf_targets_sign
    }
    if( length(tf_targets_sign) == 0 )
      next()
    # Define regulation weight
    if( ! is.null(net$weight) ){
      tf_targets_weight = subset(net, TF == tf)$weight
    }else{
      tf_targets_weight = rep(1, length(tf_targets))
    }
    si = length(tf_targets)
    # Build VIPER regulon
    tf_name = paste(tf, n, sep = ' - ')
    if( ! is.null(net$regulon_id) )
      tf_name = paste(subset(net, TF == tf)$regulon_id[1], n, sep = ' - ')
    viper_regulon[[tf_name]]$tfmode = tf_targets_sign # NOTE: ttmode has to be the first element of the list. Otherwise msviper will fail.
    names(viper_regulon[[tf_name]]$tfmode) = tf_targets
    viper_regulon[[tf_name]]$likelihood = tf_targets_weight
  }
  save(viper_regulon, file = paste('data/TF_target_sources/', gsub('.sif', '_viperRegulon.rdata', n) , sep = '') )
}







### --------------------------------------------------------------------------------------------------------------
# Format omniphat with scores
net = read.csv('data/TF_target_sources/omnipath_scores/database.csv', stringsAsFactors = F)
viper_regulon = list()
# Remove self regulation
tfs = unique(net$TF)
message('TFs: ', length(tfs))
if( length(grep('_', net$TF)) > 0 ){
  net$regulon_id = net$TF
  net$TF = sapply(strsplit(net$TF, split = '_'), head, 1)
}
net = subset(net, TF != target)
# Remove empty targets
net = subset(net, target != "")
# Define TFs of interest
tfs_of_interest = names(which(table(net$TF) >= Nmin))
message('TFs with (n>=', Nmin,'): ', length(tfs_of_interest))
for (tf in tfs_of_interest ){
  tf_targets = subset(net, TF == tf)$target
  # Define TF mode
  tf_targets_sign = subset(net, TF == tf)$effect
  if( tf %in% repressors ) # unless effect equals 0, check if tf is repressor
    tf_targets_sign[ tf_targets_sign == 0 ] = - 1
  if( ! tf %in% repressors ) # unless effect equals 0, check if tf is repressor
      tf_targets_sign[ tf_targets_sign == 0  ] = 1
  # Define regulation weight
  tf_targets_weight = c(1, 0.8, 0.6, 0.4, 0.2)[ match(subset(net, TF == tf)$score, LETTERS[1:5]) ]
  si = length(tf_targets)
  # Build VIPER regulon
  tf_name = paste(tf, 'omnipath_scores/all_weighted/network.sif', sep = ' - ')
  viper_regulon[[tf_name]]$tfmode = tf_targets_sign # NOTE: ttmode has to be the first element of the list. Otherwise msviper will fail.
  names(viper_regulon[[tf_name]]$tfmode) = tf_targets
  viper_regulon[[tf_name]]$likelihood = tf_targets_weight
}
save(viper_regulon, file = 'data/TF_target_sources/omnipath_scores/all_weighted/network_viperRegulon.rdata')






### --------------------------------------------------------------------------------------------------------------
# Load REAL ARACNE tissue specific networks + format for viper input
aracne2viper_tissueSpecificNet = function(tissue){
  aracne_file = paste('data/TF_target_sources/inferredGTEx/aracne/', tissue, '/network.txt', sep = '')
  expression_file = paste('data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/', tissue, '.txt', sep = '')
  
  A = read.delim(aracne_file, stringsAsFactors = F)
  if( nrow(A) == 0 )
    return(NULL)
  write.table(A, file = '~/tmp/network2viper.txt', sep='\t', quote = F, col.names = F, row.names = F)
  
  E = read.delim(file = expression_file, stringsAsFactors = F)
  rownames(E) = E$gene
  E = as.matrix(E[,-1])
  
  regul = aracne2regulon(afile = '~/tmp/network2viper.txt', eset = E, verbose = F, format = '3col')
}
for (tissue in list.dirs('data/TF_target_sources/inferredGTEx/aracne/', full.names = F)[-1]){
  message(tissue)
  tissue_regulon = aracne2viper_tissueSpecificNet(tissue)
  if( ! is.null(tissue_regulon) ){
    names(tissue_regulon) = paste(names(tissue_regulon), gsub(' ', '_', tissue) )
    save(tissue_regulon, file = paste('data/TF_target_sources/inferredGTEx/aracne/', tissue, '/', tissue, '_viperRegulon.rdata' , sep = '') )
  }
}


