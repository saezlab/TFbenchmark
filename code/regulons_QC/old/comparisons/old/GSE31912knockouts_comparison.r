rm(list = ls())
home = '~/googledrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


repressors = load_repressors()
nonreliable_perturbations = read.delim('data/TF_target_sources/perturbations/GEO/GSE31912_MCF7_TFknockdown_nonreliableTFperturbation.txt', stringsAsFactors = F, header = F)[,1]


fi = list.files('data/TF_activities/perturbations/', recursive = T, pattern = 'GSE31912')
activities_list = list()
for (f in fi){
  activities_list[[f]] = get(load(paste('data/TF_activities/perturbations/', f, sep = '')))  %>% apply( ., 2, rank )
  activities_list[[f]] = activities_list[[f]] / max(activities_list[[f]])
  activities_list[[f]][ rownames(activities_list[[f]]) %in% repressors,  ] = 1 - activities_list[[f]][ rownames(activities_list[[f]]) %in% repressors,  ]
  
}

df = melt(activities_list)
df$perturbation = sapply(strsplit(as.character(df$Var2), split = '_'), tail, 1)
df$TF =   toupper( as.character(df$Var1) )
df = subset(df, TF == perturbation)
df = subset(df, ! perturbation %in% nonreliable_perturbations )

df$study = sapply(strsplit(as.character(df$L1), split = '\\.'), head, 1)
df$dataset = sapply(strsplit(as.character(df$L1), split = '\\.'), function(x)  paste(x[-1], collapse = '.') ) %>% gsub('.activities.txt', '', .)
df$type = sapply(strsplit(df$dataset, split = '\\.'), head, 1)
df$name = sapply(strsplit(gsub('.activities.rdata', '', df$dataset), split = '\\.'), tail, 1)

frequent_tfs = names(which(table(df$TF)>20))
length(frequent_tfs)
ddf = ddply(subset(df, TF %in%  frequent_tfs), c('dataset', 'name'), summarise, median = median(value))  
df$name = factor(df$name, levels = ddf$name[ order(ddf$median) ])
p2 = ggplot(subset(df, TF %in% frequent_tfs), aes(y = value, x = name, fill = type) ) + geom_boxplot() + 
  scale_fill_brewer(palette = 'Set2') +
  geom_hline(yintercept = 0.5, color = 'red') + ylab('relative activity of the perturbed TF') + xlab('regulon') +
  theme_bw(18) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) 
ddf = ddply(subset(df, TF %in%  frequent_tfs), c('dataset', 'name', 'type'), summarise, n = length(unique(TF)))  
p1 = ggplot(ddf, aes(y=n, x=name, fill = type)) + geom_bar(stat='identity') + 
  scale_fill_brewer(palette = 'Set2') +
  theme_bw(18) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))  + ylab('number of TFs')
p1
plot_grid(p2, p1, ncol = 1, align = 'v')
ggsave('data/regulons_QC/pertubations/GSE31912/relativeActivity.png', dpi = 300, width = 12, height = 15)

ggplot(subset(df, TF %in% frequent_tfs), aes(y = value, x = TF) ) + geom_boxplot() + facet_wrap(~type, scales = 'free_x') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  geom_hline(yintercept = 0.1, color = 'red')
ggsave('data/regulons_QC/pertubations/GSE31912/relativeActivity_perTF.png', dpi = 300, width = 10, height = 7)
