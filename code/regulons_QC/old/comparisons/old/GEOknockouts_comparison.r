rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(plyr)


repressors = load_repressors()


plotNES = function(df, regulon_filter = NULL, TFs_filter= NULL, mycolors){
  if(! is.null(regulon_filter))
    df = subset(df, as.character(regulon_type) %in% regulon_filter)
  if(! is.null(TFs_filter))
    df = subset(df, TF %in% TFs_filter )
  # # get random distribution
  # random_q = quantile(df$NES[ df$is_random ], na.rm = T)
  # plot NES distributions per database
  ddf_ttest = ddply(df , c('regulon_shortname', 'regulon_type'), function(x){
    t.test(x$NES[ ! x$is_random ], x$NES[ x$is_random ], alternative = 'greater')$p.value}
    )
  ddf_ttest$significance = ''
  ddf_ttest$significance[ ddf_ttest$V1 < 0.05 ] = '*'
  ddf_ttest$significance[ ddf_ttest$V1 < 0.01 ] = '**'
  ddf_ttest$significance[ ddf_ttest$V1 < 0.001 ] = '***'
  
  df$regulon_shortname = factor(df$regulon_shortname, levels = ddf_ttest$regulon_shortname[ order(ddf_ttest$V1, decreasing = F) ] )
  p2 = ggplot(df, aes(y = NES, x = regulon_shortname, fill = randomization_type) ) + 
    geom_boxplot(outlier.size = 0.7, alpha = 0.5) + # NES boxplot
    geom_text(data = ddf_ttest, aes(x = regulon_shortname, label = significance), y = Inf, vjust = 1.5, color = 'red3', size = 5, inherit.aes = F) + # add significance
    # geom_hline(yintercept = random_q[2], color = 'red3', linetype = 'dashed') + geom_hline(yintercept = random_q[3], color = 'red3') + geom_hline(yintercept = random_q[4], color = 'red3', linetype = 'dashed') +  # add random quantiles
    ylab('relative activity of the perturbed TF (rank)') + xlab('regulon') +
    mytheme + facet_grid(regulon_type~., scales = 'free') +
    # scale_fill_manual(values = mycolors) + 
    scale_fill_manual(values = c('gray20', 'red', 'red4') ) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  p2
  
  ddf = ddply(subset(df, ! is_random ), c('regulon_shortname', 'regulon_type'), summarise, n = length(unique(experiment)))  
  p1 = ggplot(ddf, aes(y=n, x=regulon_shortname, fill = regulon_type)) + geom_bar(stat='identity', position = 'dodge') + 
    mytheme + scale_fill_manual(values = mycolors) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))  + ylab('TF perturbations covered')
  
  return(list(NES_boxplot=p2, n_barplot=p1))
}





perturbation_files = list.files(path = 'data/regulons_QC/pertubations/GEO/gsea', pattern = 'fgsea', full.names = T)
activities_list = list()
for (f in perturbation_files){
  message(f)
  gsea_results = get(load(f))[, c('pathway', 'NES', 'ES', 'size', 'pval', 'padj') ]
  gsea_results$TF = sapply(strsplit(gsea_results$pathway, split = ' - '), head, 1)
  gsea_results$regulon = sapply(strsplit(gsea_results$pathway, split = ' - '), tail, 1)
  gsea_results$is_random = F
  gsea_results$randomization_type = 'none'
  gsea_results$is_random[ grep('random', gsea_results$pathway) ] = T
  gsea_results$randomization_type[ grep('random_unb', gsea_results$pathway) ] = 'unbalanced'
  gsea_results$randomization_type[ grep('random_bal', gsea_results$pathway) ] = 'balanced'
  gsea_results$regulon_type = sapply(strsplit(gsea_results$regulon, split = '\\/'), head, 1)
  gsea_results$NES[ gsea_results$TF %in% repressors  ] = gsea_results$NES[ gsea_results$TF %in% repressors  ] * (-1)
  gsea_results = gsea_results[ order(gsea_results$NES, decreasing = T) ]
  gsea_results$rank = 1:nrow(gsea_results) / nrow(gsea_results)
  gsea_results$experiment = gsub('data/regulons_QC/pertubations/GEO/gsea/fgsea_', '', f) %>% gsub('.rdata', '', .)
  gsea_results$perturbed_TF = sapply(strsplit(gsea_results$experiment, split = '\\.'), head, 1)
  gsea_results$regulon_shortname = sapply(strsplit(as.character(gsea_results$regulon), split = '\\/'), function(x) x[2])
  gsea_results$regulon_shortname[ gsea_results$regulon_shortname == 'scanning_output' ] = sapply(strsplit(gsea_results$regulon[ gsea_results$regulon_shortname == 'scanning_output' ], split = "\\/"), tail, 1)
  gsea_results$regulon_shortname = gsub('network_', '', gsea_results$regulon_shortname) %>% gsub('.sif', '', .) %>% gsub('_network', '', .) %>% gsub('network', '', .)
  gsea_results$regulon_shortname[ gsea_results$regulon_type == 'reverse_engenieered' ] = paste('>', gsea_results$regulon_shortname[ gsea_results$regulon_type == 'reverse_engenieered' ], 'tissues')
  activities_list[[f]] = gsea_results
}





# take random perturbations pairs
df = melt(activities_list, id.vars = names(activities_list[[1]]))
summary(df$NES)
df$is_TF_perturbed = df$TF == df$perturbed_TF
df$randomization_type = factor(df$randomization_type, levels = c('none', 'balanced', 'unbalanced') )
ggplot(df, aes(y = NES, x= is_TF_perturbed, fill = randomization_type)) + geom_boxplot(alpha = 0.8) + mytheme + scale_fill_manual(values = my_color_palette$clear[3:1])
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_distribution.png', dpi = 300)



# check number of perturbations
length(unique(subset(df, !is.na(NES))$experiment)) # total perturbations
length(unique(subset(df, !is.na(NES))$regulon)) # total regulon sources



# select paired TF - TF pertrubation expriment
df = subset(df, is_TF_perturbed)
hist(sapply(unique(df$experiment), function(e) nrow(subset(df, ! is_random & experiment == e)) ), xlab = 'databases covering the perturbed TF', main = 'perturbations coverage', breaks = 44)
hist(sapply(unique(df$regulon), function(r) nrow(subset(df, ! is_random & regulon == r)) ), xlab = 'TFs per database', main = 'TFs coverage across databases', breaks = 173)






# plot TFBS predictions
myplots = plotNES(df, 
                  regulon_filter = c('TFBS_scanning'), 
                  TFs_filter = intersect(df$TF[ grep('hoco', df$regulon_shortname)], df$TF[ grep('jasp', df$regulon_shortname)] ),
                  mycolors = my_color_palette$clear[c(3)] )
myplots$NES_boxplot
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_TFBS_scanning.png', dpi = 300, width = 9, height = 7)



# plot ARACNe
myplots = plotNES(df, 
                  regulon_filter = c('reverse_engenieered'), 
                  TFs_filter = intersectSeveral(lapply(c('> 2 tissues', '> 5 tissues', '> 10 tissues'), function(f) subset( df, regulon_shortname == f)$TF ) ),
                  mycolors = my_color_palette$clear[c(1)])
myplots$NES_boxplot
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_reverse_engenieered.png', dpi = 300, width = 5, height = 7)


# plot curated
myplots = plotNES(df, 
                  regulon_filter = c("curated_databases"), 
                  mycolors = my_color_palette$clear[c(1)])
myplots$NES_boxplot
plot_grid(myplots$NES_boxplot, myplots$n_barplot, ncol = 1, align = 'v', rel_heights = c(6,4))
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_curated.png', dpi = 300, width = 9, height = 15)



# plot all
myplots = plotNES(df, mycolors = my_color_palette$clear[c(7,1, 4,3,6,5,2)])
myplots$NES_boxplot
plot_grid(myplots$NES_boxplot, myplots$n_barplot, ncol = 1, align = 'v', rel_heights = c(7,3))
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities.png', dpi = 300, width = 12, height = 15)





frequent_tfs = names(which(table(unique(df[, c('TF', 'regulon_shortname') ])$TF) > 15))
ddf = ddply(subset(df, TF %in% frequent_tfs), c('regulon_shortname'), summarise, median = median(NES, na.rm = T))
ddf$regulon_shortname = as.character(ddf$regulon_shortname)
df$regulon_shortname = factor(df$regulon_shortname, levels = ddf$regulon_shortname[ order(ddf$median) ] )
p2 = ggplot(subset(df, TF %in% frequent_tfs), aes(y = NES, x = regulon_shortname, fill = regulon_type) ) + geom_boxplot() + 
  geom_hline(yintercept = random_q[2], color = 'gray') + geom_hline(yintercept = random_q[3], color = 'red2') + geom_hline(yintercept = random_q[4], color = 'gray') + 
  ylab('relative activity of the perturbed TF (rank)') + xlab('regulon') +
  mytheme + scale_fill_manual(values = my_color_palette$clear[c(7,1, 4,3,6,5,2)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
p2
ddf = ddply(subset(df, TF %in% frequent_tfs), c('regulon_shortname', 'regulon_type'), summarise, n = length(unique(experiment)))  
p1 = ggplot(ddf, aes(y=n, x=regulon_shortname, fill = regulon_type)) + geom_bar(stat='identity', position = 'dodge') + 
  mytheme + scale_fill_manual(values = my_color_palette$clear[c(7,1, 4,3,6,5,2)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))  + ylab('TF perturbations covered')
p1
plot_grid(p2, p1, ncol = 1, align = 'v')
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_freq10TFs.png', dpi = 300, width = 12, height = 15)





ggplot(subset(df, TF %in% frequent_tfs), aes(y = NES, x = TF) ) + geom_boxplot() + facet_wrap(~regulon_type, scales = 'free_x') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  geom_hline(yintercept = random_q[2], color = 'gray') + geom_hline(yintercept = random_q[3], color = 'red2') + geom_hline(yintercept = random_q[4], color = 'gray') 
ggsave('data/regulons_QC/pertubations/GEO/relativeActivity_perTF.png', dpi = 300, width = 10, height = 7)

