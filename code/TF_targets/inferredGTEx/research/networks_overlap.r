rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)



# Load networks
aracne_path = 'data/TF_target_sources/reverse_engenieered/aracne/'
networks = list()
for (tissue in list.files(path = aracne_path, full.names = T, recursive = T, pattern = 'network')){
  tissue_name = unlist(strsplit(tissue, split = '/'))[6]
  message(tissue_name)
  networks[[tissue_name]] = read.delim(file = tissue, stringsAsFactors = F, header = T)
}
networks = Filter(Negate(function(x){ nrow(x) == 0}), networks)
length(networks)




# Define TFs
overlaping_TFs = Reduce(intersect, lapply(networks, function(x) x$Regulator )[1:3])
length(overlaping_TFs)



# Check global overlap
df = melt(networks, id.vars = names(networks[[1]]))
df$interaction = paste(df$Regulator, df$Target, sep = '-')
interactions_count = table(df$interaction)
# plot network size as num tissues threshold
plot_df = data.frame(interactions_count=1:30, count = sapply(1:30, function(i) sum((interactions_count >= i ) + 0) ) )
plot_df$interactions_count = factor(plot_df$interactions_count)
head(plot_df)
p1 = ggplot(plot_df, aes(y = count, x = interactions_count, label = count) ) + geom_bar(stat= 'identity', fill = '#00AFBB') + theme_light(18) + 
  ggtitle('network size') +
  geom_text(angle = 45, hjust = 0) + ylab('count') + xlab('n tissues threshold')
p1
# plot TFs keept as num tissues threshold
size_summaries_per_filter = list()
targets_summaries_per_filter = list()
TFsize_summaries_per_filter = list()
targetssize_summaries_per_filter = list()
for( i in 1:30){
  ddf = subset(df, interaction %in% names(which(interactions_count>=i)) )
  TFfreqs = table(ddf$Regulator)
  targetfreqs = table(ddf$Target)
  size_summaries_per_filter[[ paste('n tissue', i) ]]$regulon_size_3 = length(which(TFfreqs > 2))
  size_summaries_per_filter[[ paste('n tissue', i) ]]$regulon_size_5 = length(which(TFfreqs > 4))
  size_summaries_per_filter[[ paste('n tissue', i) ]]$regulon_size_10 = length(which(TFfreqs > 9))
  size_summaries_per_filter[[ paste('n tissue', i) ]]$regulon_size_20 = length(which(TFfreqs > 19))
  TFsize_summaries_per_filter[[ paste('n tissue', i) ]] = TFfreqs
  targetssize_summaries_per_filter[[ paste('n tissue', i) ]] = targetfreqs
  targets_summaries_per_filter[[ paste('n tissue', i) ]] = length(unique(ddf$Target))
}
plot_df = melt(size_summaries_per_filter)
plot_df$L2 = factor(plot_df$L2, levels = unique(plot_df$L2), labels =  gsub('regulon_size_', '>=', unique(plot_df$L2)) )
plot_df$L1 = factor(plot_df$L1, levels = unique(plot_df$L1), labels = 1:30)
p2 = ggplot(plot_df, aes(x=L1, y=value, fill = L2) ) + geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_brewer(palette = 'Dark2') + 
  theme_light(18) + ylab('count') + xlab('n tissues threshold') + theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.title = element_blank()) + ggtitle('TFs')
# plot targets as num tissues threshold
plot_df = melt(targets_summaries_per_filter)
plot_df$L1 = factor(plot_df$L1, levels = unique(plot_df$L1), labels = 1:30) 
p3 = ggplot(plot_df, aes(y = value, x = L1) ) + geom_bar(stat= 'identity', fill = '#E7B800') + theme_light(18) + 
  ggtitle('targets')  + ylab('count') + xlab('n tissues threshold')
# plot regulon sizes as num tissues threshold
plot_df = melt(TFsize_summaries_per_filter)
plot_df$L1 = factor(plot_df$L1, levels = c(unique(plot_df$L1), 28:30), labels = 1:30) 
p4 = ggplot(plot_df, aes(x= L1, y= value) ) + geom_boxplot()  + scale_y_log10() +
  theme_light(18) + ylab('targets per TF') + xlab('n tissues threshold') + ggtitle('TF regulon size')
# plot regulon sizes as num tissues threshold
plot_df = melt(targetssize_summaries_per_filter)
plot_df$L1 = factor(plot_df$L1, levels = c(unique(plot_df$L1), 28:30), labels = 1:30) 
p5 = ggplot(plot_df, aes(x= L1, y= value) ) + geom_boxplot()  + scale_y_log10() +
  theme_light(18) + ylab('TFs per target') + xlab('n tissues threshold') + ggtitle('TFs per target')
myplot = plot_grid(p1, p2, p3, p4, p5, ncol = 1, align = 'v')  
save_plot(myplot, filename = 'data/TF_target_sources/reverse_engenieered/aracne/filter_summaries.png', base_height = 15, base_width = 10)


# intersect tissues
df$value = 1
m = dcast(df, formula = interaction~L1, fill = 0)
names(m)[1] = 'Identifier'
upset(m[ which(rowSums(m[,-1]) > 1) , ], sets = colnames(m)[-1], nintersects = 40, main.bar.color = 'gray20', order.by = "freq", number.angles = 0, text.scale = 1.3, att.color = 'black', mb.ratio = c(.7, .3))
ggsave(filename = 'data/TF_target_sources/reverse_engenieered_networks/aracne/overlaps.png')






# Build regulons
regulons = list()
for (tissue in names(networks)){
  message(tissue)
  for (tf in overlaping_TFs) {
    regulons[[tissue]][[tf]] = networks[[tissue]]$Target[ networks[[tissue]]$Regulator == tf ]
  }
}




# compute overlap
overlaps = list()
roverlaps = list()
M = matrix(0, nrow = length(networks), ncol = length(networks), dimnames = list(names(networks), names(networks)))
for (tf in sort(overlaping_TFs)){
  message(tf)
  overlaps[[tf]] = M
  roverlaps[[tf]] = M
  for (tissue1 in names(networks)){
    targets1 = regulons[[tissue1]][[tf]]
    for (tissue2 in names(networks)){
      targets2 = regulons[[tissue2]][[tf]]
      overlaps[[tf]][tissue1, tissue2 ] = length(intersect(targets1, targets2)) / length(unique(c(targets1, targets2))) 
      rtargets = regulons[[tissue2]][[sample(overlaping_TFs, 1)]]
      roverlaps[[tf]][tissue1, tissue2 ] = length(intersect(targets1, rtargets)) / length(unique(c(targets1, rtargets))) 
    }
  }
  diag(overlaps[[tf]]) = NA
  diag(roverlaps[[tf]]) = NA
}







# plots 
# compare same TF vs other TFs similarity
df = melt(list(observed=overlaps, expected = roverlaps))
sub_df = df[ ! is.na(df$value),  ]
ddf_ttest = ddply(sub_df, c('L2'), function(x) t.test(x$value[ x$L1 == 'observed'], x$value[ x$L1 != 'observed'], alternative = 'greater')$p.value  )
ddf_ttest = ddf_ttest[ order(ddf_ttest$V1) , ]
table(p.adjust(ddf_ttest$V1, method = 'fdr') < 0.01)


sub_df$L2 = factor(sub_df$L2, levels = ddf_ttest$L2 )
sub_df$level = 1
sub_df$level[ sub_df$L2 %in% ddf_ttest$L2[1:300] ] = 1
sub_df$level[ sub_df$L2 %in% ddf_ttest$L2[301:600] ] = 2
sub_df$level[ sub_df$L2 %in% ddf_ttest$L2[601:900] ] = 3
sub_df$level[ sub_df$L2 %in% ddf_ttest$L2[901:1205] ] = 4
ggplot(sub_df[ grep('^Z[H-Z]', sub_df$L2, invert = T) , ], aes(x= L2, y = value, color = L1) ) + geom_boxplot(position = 'dodge', outlier.size = .4) + scale_color_brewer(palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =.5)) + facet_wrap(facets = ~ level, ncol = 1, scales = 'free_x') + theme(legend.position = 'none')




# compare similarity when TF is or is not expressed
# 1) Differential expression
deTFs = get(load(file = 'data/TF_target_sources/reverse_engenieered_networks/tf_differential_expression.rdata'))
deTFs = deTFs[ deTFs$adj.P.Val < 0.05 & deTFs$logFC > 0 , ]
df = melt(overlaps)
df = df[ ! is.na(df$value),  ]
df$expressed_in_both_tissues = F
for( tf in overlaping_TFs ){
  tissues_expressed = intersect(gsub(' ', '_', deTFs$L1[ deTFs$TF == tf ]), df$Var1 )
  if( length(tissues_expressed) > 2 )
    df$expressed_in_both_tissues[ df$Var1 %in% tissues_expressed & df$Var2 %in% tissues_expressed & df$L1 == tf ] = T
  if ( length(tissues_expressed) < 3 )
    df = df[ df$L1 != tf, ]
}

ddf_ttest = ddply(df, c('L1'), function(x) t.test(x$value[ x$expressed_in_both_tissues], x$value[ ! x$expressed_in_both_tissues ], alternative = 'greater')$p.value   )
ddf_ttest = ddf_ttest[ order(ddf_ttest$V1) , ]
table(p.adjust(ddf_ttest$V1, method = 'fdr') < 0.05 )
dim(ddf_ttest)
head(ddf_ttest, 20)

df$L1 = factor(df$L1, levels = ddf_ttest$L1[ order(ddf_ttest$V1) ] )
df$level = 6
df$level[ df$L1 %in% ddf_ttest$L1[1:220] ] = 1
df$level[ df$L1 %in% ddf_ttest$L1[221:440] ] = 2
df$level[ df$L1 %in% ddf_ttest$L1[441:623] ] = 3
df$level[ df$L1 %in% ddf_ttest$L1[623:800] ] = 4
df$level[ df$L1 %in% ddf_ttest$L1[801:1000] ] = 5
ggplot(df, aes(x= L1, y = value, color = expressed_in_both_tissues) ) + geom_boxplot(position = 'dodge', outlier.size = .4) + scale_color_brewer(palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =.5)) + facet_wrap(facets = ~ level, ncol = 1, scales = 'free_x') + theme(legend.position = 'none')



# 2) 25th percentile > 5 fpkms
callTFs = get(load(file = 'data/TF_target_sources/reverse_engenieered_networks/tf_tissues_call_percentile.rdata'))
callTFs = subset(callTFs, L2 == 'active')
df = melt(overlaps)
df = df[ ! is.na(df$value),  ]
df$expressed_in_both_tissues = F
for( tf in overlaping_TFs ){
  tissues_expressed = intersect(gsub(' ', '_', callTFs$L1[ callTFs$value == tf ]), df$Var1 )
  df$expressed_in_both_tissues[ df$Var1 %in% tissues_expressed & df$Var2 %in% tissues_expressed & df$L1 == tf ] = T
  if ( length(tissues_expressed) < 3 | length(tissues_expressed) > 25 )
    df = df[ df$L1 != tf, ]
}
table(df$expressed_in_both_tissues)
ddf_ttest = ddply(df, c('L1'), function(x) t.test(x$value[ x$expressed_in_both_tissues], x$value[ ! x$expressed_in_both_tissues ], alternative = 'greater')$p.value   )
ddf_ttest = ddf_ttest[ order(ddf_ttest$V1) , ]
table(p.adjust(ddf_ttest$V1, method = 'fdr') < 0.05 )
dim(ddf_ttest)
head(ddf_ttest, 20)

df$L1 = factor(df$L1, levels = ddf_ttest$L1[ order(ddf_ttest$V1) ] )
df$level = 2
df$level[ df$L1 %in% ddf_ttest$L1[1:200] ] = 1
ggplot(df, aes(x= L1, y = value, color = expressed_in_both_tissues) ) + geom_boxplot(position = 'dodge', outlier.size = .4) + scale_color_brewer(palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =.5)) + facet_wrap(facets = ~ level, ncol = 1, scales = 'free_x') + theme(legend.position = 'none')


