rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


repressors = load_repressors()
del = load_homDel()
del = del[ rowSums(del) > 0, colSums(del) > 0 ]
del = del[ rowSums(del) > 0, ]
del = del[ , colSums(del) > 0 ]
dim(del)


fi = list.files('data/regulons_QC/cell_lines/', recursive = T)





activities_list = list()
ranks_list = list()
for (f in fi){
  x = get(load(paste('data/regulons_QC/cell_lines/', f, sep = ''))) 
  tfs = intersect(rownames(x), rownames(del) ) 
  samples = intersect(colnames(x), colnames(del))
  x = x[ tfs , samples ] 
  if( ! is.matrix(x) )
    next
  if( nrow(x) < 2 )
    next
  activities_list[[f]] = t(apply(x, 1, rank ))
  activities_list[[f]] = activities_list[[f]] / max(activities_list[[f]])
  activities_list[[f]][ rownames(activities_list[[f]]) %in% repressors,  ] = 1 - activities_list[[f]][ rownames(activities_list[[f]]) %in% repressors,  ]
  ranks_list[[f]] = activities_list[[f]][ del[tfs, samples] == 1 ]
  ranks_list[[f]] = data.frame(ranks = ranks_list[[f]], TF = unlist(apply(del[tfs, samples], 2, function(idx) rownames(del[tfs, ])[idx==1]  )))
}





# plot
df = melt(ranks_list, id.vars = names(ranks_list[[1]]))
df$L1 = gsub('reverse_engenieered.activities.rdata', 'reverse_engenieered.aracne.activities.rdata', df$L1)
df$L1 = gsub('old_dorotheav1.activities.rdata', 'old_dorotheav1.DoRothEAv1.activities.rdata', df$L1)
df$type = sapply(strsplit(as.character(df$L1), split = '\\.'), head, 1)
df$dataset = sapply(strsplit(as.character(df$L1), split = '\\.'), function(x)  paste(x[-1], collapse = '.') ) %>% gsub('.activities.rdata', '', .)
df$name = sapply(strsplit(gsub('.activities.rdata', '', df$dataset), split = '\\.'), tail, 1)




ddf = ddply(df, c('name', 'type'), summarise, median = median(ranks))  
df$name = factor(df$name, levels = ddf$name[ order(ddf$median) ])
p2 = ggplot(df, aes(y=ranks, x=name, fill = type)) + geom_boxplot() + 
  scale_fill_brewer(palette = 'Set2') +
  geom_hline(yintercept = 0.5, color = 'red') +
  theme_bw(18) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + ylab('relative activity\' rank')
p2
ddf = ddply(df, c('name', 'type'), summarise, n = length(unique(TF)))  
p1 = ggplot(ddf, aes(y=n, x=name, fill = type)) + geom_bar(stat='identity') + 
  scale_fill_brewer(palette = 'Set2') +
  theme_bw(18) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))  + ylab('number of TFs')
p1
plot_grid(p2, p1, ncol = 1, align = 'v')
ggsave('data/regulons_QC/homDel_TFactivity.png', dpi = 300, width = 12, height = 15)


ggplot(subset(df, TF %in% subset(df, type = 'consensus')$TF), aes(y=ranks, x=TF)) + geom_boxplot() + 
  theme_bw(18) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) 

