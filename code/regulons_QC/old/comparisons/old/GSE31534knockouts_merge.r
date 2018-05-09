rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
get_percentile = function(ob_nes, ra_nes, value = 'pvalue_1tailed'){
  per = ecdf(ra_nes)(ob_nes)
  pvalue_1tailed = per
  if(  value == 'pvalue_1tailed' )
    return(pvalue_1tailed)
  pvalue_2tailed = 1-abs(per - 0.5)*2
  return(pvalue_2tailed)
}



# Load unsigned activities
tissue = 'Skin'
activity_files = list.files('data/regulons_QC/pertubations/GSE31534/activities', full.names = T) %>% grep(.,  pattern = 'RANDOM', value = T, invert= T) %>% grep(.,  pattern = 'aracne', value = T, invert= T)
# build unsigned randoms
random_activity_files = list.files('data/regulons_QC/pertubations/GSE31534/activities', pattern = 'RANDOM_NETWORKS_coexpressionPriors', full.names = T) %>% grep('tissueSpecific', ., value =T, invert = T)
randoms = vector("list", length = 1000)
for (ra in random_activity_files){
  nes = get(load(ra))
  nes = nes[ nes[,1] <= 1000,  ]
  for( si in unique(nes[,1]) ){
    randoms[[si]] = c(randoms[[si]], nes[ nes[,1] == si , -1])
  }
}
sapply(randoms, length)
names(randoms) = 1:1000
randoms = randoms[ sapply(randoms, length) > 100 ]
# normalize unsigned activities
observed = list()
for (o in activity_files){
  message(o)
  nes = get(load(o))
  nes = nes[ as.character(nes[,'Size']) %in%  names(randoms) , ]
  if( ! is.matrix(nes) )
    next
  if( nrow(nes) == 0)
    next
  if( is.matrix(nes) ) {
    p1t = nes[,-1]
    p1t_rank = nes[,-1]
    p2t = nes[,-1]
    p2t_rank = nes[,-1]
    for( tf in rownames(nes) ){
      tf_size = as.character(nes[tf,'Size'])
      p1t[tf,] = get_percentile(p1t[tf,], randoms[[tf_size]], value = 'pvalue_1tailed')
      p2t[tf,] = get_percentile(p2t[tf,], randoms[[tf_size]], value = 'pvalue_2tailed')
      p1t_rank[tf,] = rank(p1t[tf,]) / max(rank(p1t[tf,])) 
      p2t_rank[tf,] = rank(p2t[tf,]) / max(rank(p2t[tf,])) 
    }
    df = melt(p1t, value.name = 'pvalue_1tailed')
    df$pvalue_2tailed = melt(p2t, value.name = 'val')$val
    df$pvalue1tailed_rank = melt(p1t_rank, value.name = 'val')$val
    df$pvalue2tailed_rank = melt(p2t_rank, value.name = 'val')$val
    df$nes = melt(nes[ , colnames(p1t) ], value.name = 'nes')$nes
    df$regulon = tail(unlist(strsplit(o, split = "/")),1)
    df$regulon_type = head(unlist(strsplit(df$regulon[1], split = "\\.")), 1)
    df$regulon_shortname = head(tail(unlist(strsplit(df$regulon[1], split = "\\.")), 3), 1)
    df$Var1 = as.character(df$Var1)
    df$TF = sapply(strsplit(df$Var1, split = ' '), head, 1)
    if( df$regulon_shortname[1] %in% c('consensus', 'scanning_output') )
      df$regulon_shortname = tail(unlist(strsplit(df$regulon[1], split = "\\.")), 3)[2] %>% gsub('_viperRegulon', '', .)
    if( df$regulon_shortname[1] == 'ReMap' )
      df$regulon_shortname = tail(unlist(strsplit(df$regulon[1], split = "\\.")), 3)[2] %>% gsub('_viperRegulon', '', .) %>% paste('ReMap_', ., sep ='')
    if( df$regulon_shortname[1] == 'inferredGTEx' )
      df$regulon_shortname = tail(unlist(strsplit(df$regulon[1], split = "\\.")), 3)[2] %>% gsub('_viperRegulon', '', .) %>% paste('inferredGTEx_', ., sep ='')
    df$regulon_shortname = gsub('network', '', df$regulon_shortname) %>% gsub('_$', '', .) %>% gsub('^_', '', .) 
    df$perturbed_gene = sapply(strsplit(as.character(df$Var2), split = '_'), tail, 1)
    observed[[o]] = df
  }
}
save(observed, file = 'data/regulons_QC/pertubations/GSE31534/activity_comparison_results.r')







################################################ 
# Load signed activities
random_activity_file = list.files('data/regulons_QC/pertubations/GSE31534/activities', pattern = 'RAN', full.names = T) %>% grep(tissue, ., value = T)
nes = get(load(random_activity_file))
nes = nes[ nes[,1] <= 1000,  ]
nes[,1] = sapply(strsplit(rownames(nes), split = ' - '), head, 1)
randoms = list()
for( si in unique(nes[,1]) ){
  randoms[[si]] = as.numeric(nes[ nes[,1] == si, -1])
}
# normalize unsigned activities
activity_files = list.files('data/regulons_QC/pertubations/GSE31534/activities', full.names = T) %>% grep(.,  pattern = tissue, value = T) %>% grep(.,  pattern = 'aracne', value = T)
for (o in activity_files){
  message(o)
  nes = get(load(o))
  nes = nes[ rownames(nes) %in%  names(randoms) , ]
  if( ! is.matrix(nes) )
    next
  if( nrow(nes) == 0)
    next
  if( is.matrix(nes) ) {
    p1t = nes[,-1]
    p1t_rank = nes[,-1]
    p2t = nes[,-1]
    p2t_rank = nes[,-1]
    for( tf in rownames(nes) ){
      p1t[tf,] = get_percentile(p1t[tf,], randoms[[tf]], value = 'pvalue_1tailed')
      p2t[tf,] = get_percentile(p2t[tf,], randoms[[tf]], value = 'pvalue_2tailed')
      p1t_rank[tf,] = rank(p1t[tf,]) / max(rank(p1t[tf,])) 
      p2t_rank[tf,] = rank(p2t[tf,]) / max(rank(p2t[tf,])) 
    }
    df = melt(p1t, value.name = 'pvalue_1tailed')
    df$Var1 = as.character(df$Var1)
    df$pvalue_2tailed = melt(p2t, value.name = 'val')$val
    df$pvalue1tailed_rank = melt(p1t_rank, value.name = 'val')$val
    df$pvalue2tailed_rank = melt(p2t_rank, value.name = 'val')$val
    df$nes = melt(nes[ , colnames(p1t) ], value.name = 'nes')$nes
    df$regulon = tail(unlist(strsplit(o, split = "/")),1)
    df$regulon_type = head(unlist(strsplit(df$regulon[1], split = "\\.")), 1)
    df$regulon_shortname = head(tail(unlist(strsplit(df$regulon[1], split = "\\.")), 3), 1)
    df$TF = sapply(strsplit(df$Var1, split = ' '), head, 1)
    if( df$regulon_shortname[1] %in% c('consensus', 'scanning_output') )
      df$regulon_shortname = tail(unlist(strsplit(df$regulon[1], split = "\\.")), 3)[2] %>% gsub('_viperRegulon', '', .)
    if( df$regulon_shortname[1] == 'ReMap' )
      df$regulon_shortname = tail(unlist(strsplit(df$regulon[1], split = "\\.")), 3)[2] %>% gsub('_viperRegulon', '', .) %>% paste('ReMap_', ., sep ='')
    if( df$regulon_shortname[1] == 'inferredGTEx' )
      df$regulon_shortname = tail(unlist(strsplit(df$regulon[1], split = "\\.")), 3)[2] %>% gsub('_viperRegulon', '', .) %>% paste('inferredGTEx_', ., sep ='')
    df$regulon_shortname = gsub('network', '', df$regulon_shortname) %>% gsub('_$', '', .) %>% gsub('^_', '', .) 
    df$perturbed_gene = sapply(strsplit(as.character(df$Var2), split = '_'), tail, 1)
    observed[[o]] = df
  }
}
save(observed, file = 'data/regulons_QC/pertubations/GSE31534/activity_comparison_results.rdata')
