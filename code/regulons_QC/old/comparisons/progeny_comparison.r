rm(list = ls())
home = '~/Google Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


repressors = load_repressors()
progeny_scores = t(get(load('../pathway_activities/PROGENY/data/celllines_progeny14PW.rdata')))



# get correlations
fi = list.files('data/TF_activities/cell_lines/', recursive = T)
correlation_list = list()
for (f in fi){
  message(f)
  x = get(load(paste('data/TF_activities/cell_lines/', f, sep = '')))  %>% apply( ., 2, rank )
  x = x / max(x)
  x[ rownames(x) %in% repressors,  ] = 1 - x[ rownames(x) %in% repressors,  ]
  M = matrix(0, ncol = nrow(x), 
             nrow = nrow(progeny_scores), dimnames = list(rownames(progeny_scores), rownames(x)) )
  sam = intersect(colnames(x), colnames(progeny_scores))
  for( pro in rownames(progeny_scores) )
    for( tf in rownames(x))
      M[pro, tf] = cor(x[tf, sam ], progeny_scores[pro, sam ])
  correlation_list[[f]] = M
}



# rank top TFs
progeny_TFs = list('AR'='Androgen', 
                   'ESR1'='Estrogen', 
                   'HIF1A'='Hypoxia',
                   'STAT1'='JAK.STAT',
                   'STAT2'='JAK.STAT', 
                   'RELA'='NFkB',
                   "NFKB1"='NFkB',
                   'ELK1'='MAPK',
                   'TP53'='p53')  
ranks = list()
for( co in names(correlation_list) ){
  x = correlation_list[[co]]
  x = x[, colnames(x) %in% names(progeny_TFs)  ]
  if( ! is.matrix(x) )
    next
  ranks[[co]] = sapply( colnames(x) , function(tf){
    po = sort(rank(x[,tf])  / nrow(x))
    15 - which(names(po) == progeny_TFs[[tf]])
    })
  ranks[[co]] = data.frame(rank=ranks[[co]], TF=names(ranks[[co]]), stringsAsFactors = F)
  ranks[[co]]$n = ncol(x)
}
ranks = Filter(Negate(function(x) ncol(x)==0),ranks)



# plot
df = melt(ranks, id.vars = names(ranks[[1]]))
df$L1 = gsub('reverse_engenieered.activities.rdata', 'reverse_engenieered.aracne.activities.rdata', df$L1)
df$L1 = gsub('old_dorotheav1.activities.rdata', 'old_dorotheav1.DoRothEAv1.activities.rdata', df$L1)
df$type = sapply(strsplit(as.character(df$L1), split = '\\.'), head, 1)
df$dataset = sapply(strsplit(as.character(df$L1), split = '\\.'), function(x)  paste(x[-1], collapse = '.') ) %>% gsub('.activities.rdata', '', .)
df$name = sapply(strsplit(gsub('.activities.rdata', '', df$dataset), split = '\\.'), tail, 1)




ddf = ddply(df, c('name', 'type'), summarise, median = median(rank))  
df$name = factor(df$name, levels = ddf$name[ order(ddf$median) ])
p2 = ggplot(df, aes(y=rank, x=name, fill = type)) + geom_boxplot() + 
  scale_fill_brewer(palette = 'Set2') +
  geom_hline(yintercept = 7, color = 'red') +
  theme_bw(18) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + ylab('relative activity\' rank')
p1 = ggplot(unique(df[, -(1:2) ]), aes(y=n, x=name, fill = type)) + geom_bar(stat='identity') + 
  scale_fill_brewer(palette = 'Set2') +
  theme_bw(18) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))  + ylab('number of TFs')
plot_grid(p2, p1, ncol = 1, align = 'v')
ggsave('data/regulons_QC/progeny_TFactivity_correlations.png', dpi = 300, width = 12, height = 15)


ggplot(df, aes(y=rank, x=TF)) + geom_boxplot() + 
  geom_hline(yintercept = 7, color = 'red') +
  theme_bw(18) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) 

