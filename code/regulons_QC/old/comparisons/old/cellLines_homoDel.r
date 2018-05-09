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




# Load deletions
del = load_homDel()
del = del[ rowSums(del) > 0, ]



# Load random activities
random_activity_files = list.files('data/regulons_QC/cell_lines/activities/', pattern = 'RANDOM_NETWORKS_coexpressionPriors', full.names = T) %>% grep('tissueSpecific', ., invert = T, value = T)
randoms = vector("list", length = 1000)
for (ra in random_activity_files){
  message(ra)
  nes = get(load(ra))
  nes = nes[ nes[,1] <= 1000,  ]
  for( si in unique(nes[,1]) ){
    randoms[[si]] = c(randoms[[si]], nes[ nes[,1] == si , -1])
  }
}
names(randoms) = 1:1000
randoms = randoms[ sapply(randoms, length) > 100 ] 


# normalize activities
activity_files = list.files('data/regulons_QC/cell_lines/activities/', full.names = T) %>% grep(.,  pattern = 'RANDOM', value = T, invert= T)  %>% grep(.,  pattern = 'aracne', value = T, invert= T)
observed = list()
for (o in activity_files){
  message(o)
  nes = get(load(o))
  nes = nes[ as.character(nes[,'Size']) %in%  names(randoms) , ]
  if( ! is.matrix(nes) )
    next
  if( nrow(nes) == 0)
    next
  rownames(nes) = sapply(strsplit(rownames(nes), ' - '), head, 1)
  nes = nes[ rownames(nes) %in% rownames(del), ]
  if( ! is.matrix(nes) )
    next
  if( nrow(nes) == 0)
    next
  message(nrow(nes))
  if( is.matrix(nes) ) {
    p1t = nes[, colnames(nes) %in% colnames(del) ]
    p1t_rank = p1t
    p2t = nes[, colnames(nes) %in% colnames(del) ]
    p2t_rank = p2t
    for( tf in rownames(nes) ){
      tf_size = as.character(nes[tf,'Size'])
      p1t[tf,] = get_percentile(p1t[tf,], randoms[[tf_size]], value = 'pvalue_1tailed')
      p2t[tf,] = get_percentile(p2t[tf,], randoms[[tf_size]], value = 'pvalue_2tailed')
      p1t_rank[tf,] = rank(p1t[tf,]) / max(rank(p1t[tf,])) 
      p2t_rank[tf,] = rank(p2t[tf,]) / max(rank(p2t[tf,])) 
    }
    message(nrow(p1t))
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
    df$del_status = melt(del[ rownames(p1t) , colnames(p1t) ], value.name = 'val', id.vars = NULL)$val
    observed[[o]] = df
  }
}
save(observed, file = 'data/regulons_QC/cell_lines/comparison_results_homDel.rdata')



