rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(plyr)



plotNES = function(df, regulon_filter = NULL, TFs_filter= NULL, mycolors){
  if(! is.null(regulon_filter))
    df = subset(df, as.character(regulon_type) %in% regulon_filter)
  if(! is.null(TFs_filter))
    df = subset(df, TF %in% TFs_filter )
  # # get random distribution
  # random_q = quantile(df$NES[ df$is_random ], na.rm = T)
  # plot NES distributions per database
  ddf_ttest = ddply(df , c('regulon_shortname', 'regulon_type'), function(x){
    pval = NULL
    pval = try(t.test(x$NES[ ! x$is_random ], x$NES[ x$randomization_type == 'balanced' ], alternative = 'greater')$p.value, silent = T)
    effectsize = effectsize_cohensD(x$NES[ ! x$is_random ], x$NES[ x$randomization_type == 'balanced' ])
    if( class(pval) == 'try-error' )
      pval = 1
    return(c(pval, effectsize))
    } )
  names(ddf_ttest)[3:4] = c('pval', 'effectsize')
  ddf_ttest = ddf_ttest[ order(ddf_ttest$effectsize, na.last = T, decreasing = T) , ]
  ddf_ttest$significance = ''
  ddf_ttest$significance[ ddf_ttest$pval < 0.05 ] = '*'
  ddf_ttest$significance[ ddf_ttest$pval < 0.01 ] = '**'
  ddf_ttest$significance[ ddf_ttest$pval < 0.001 ] = '***'
  
  df$regulon_shortname = factor(df$regulon_shortname, levels = ddf_ttest$regulon_shortname)
  p2 = ggplot(df, aes(y = NES, x = regulon_shortname, fill = randomization_type) ) + 
    geom_boxplot(outlier.size = 0.7, alpha = 0.5) + # NES boxplot
    geom_text(data = ddf_ttest, aes(x = regulon_shortname, label = significance), y = Inf, vjust = 1.5, color = 'red3', size = 5, inherit.aes = F) + # add significance
    # geom_hline(yintercept = random_q[2], color = 'red3', linetype = 'dashed') + geom_hline(yintercept = random_q[3], color = 'red3') + geom_hline(yintercept = random_q[4], color = 'red3', linetype = 'dashed') +  # add random quantiles
    ylab('relative activity of the perturbed TF\n(VIPER-NES)') + xlab('regulon') +
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





perturbation_files = list.files(path = 'data/regulons_QC/pertubations/GEO/gsea', pattern = 'viper', full.names = T)
activities_list = list()
for (f in perturbation_files){
  message(f)
  gsea_results = get(load(f))
  gsea_results$TF = sapply(strsplit(as.character(gsea_results$Regulon), split = ' - '), head, 1)
  gsea_results$regulon = sapply(strsplit(as.character(gsea_results$Regulon), split = ' - '), tail, 1)
  gsea_results$is_random = F
  gsea_results$randomization_type = 'none'
  gsea_results$is_random[ grep('random', gsea_results$Regulon) ] = T
  # gsea_results$randomization_type[ grep('random_unb', gsea_results$Regulon) ] = 'unbalanced'
  gsea_results$randomization_type[ grep('random_bal', gsea_results$Regulon) ] = 'balanced'
  gsea_results$regulon_type = sapply(strsplit(gsea_results$regulon, split = '\\/'), head, 1)
  gsea_results = gsea_results[ order(gsea_results$NES, decreasing = T) , ]
  gsea_results$experiment = gsub('data/regulons_QC/pertubations/GEO/gsea/viper_', '', f) %>% gsub('.rdata', '', .)
  gsea_results$perturbed_TF = sapply(strsplit(gsea_results$experiment, split = '\\.'), head, 1)
  gsea_results$regulon_shortname = sapply(strsplit(as.character(gsea_results$regulon), split = '\\/'), function(x) x[2])
  gsea_results$regulon_shortname[ gsea_results$regulon_shortname == 'scanning_output' ] = sapply(strsplit(gsea_results$regulon[ gsea_results$regulon_shortname == 'scanning_output' ], split = "\\/"), tail, 1)
  gsea_results$regulon_shortname[ gsea_results$regulon_shortname == 'ReMap' ] = sapply(strsplit(gsea_results$regulon[ gsea_results$regulon_shortname == 'ReMap' ], split = "\\/"), tail, 1)
  gsea_results$regulon_shortname = gsub('network_', '', gsea_results$regulon_shortname) %>% gsub('.sif', '', .) %>% gsub('_network', '', .) %>% gsub('network', '', .)
  gsea_results$regulon_shortname[ gsea_results$regulon_type == 'reverse_engenieered' ] = paste('>', gsea_results$regulon_shortname[ gsea_results$regulon_type == 'reverse_engenieered' ], 'tissues')
  activities_list[[f]] = gsea_results
}





# take random perturbations pairs
df = melt(activities_list, id.vars = names(activities_list[[1]]))
summary(df$NES)
df$is_TF_perturbed = df$TF == df$perturbed_TF
df$randomization_type = factor(df$randomization_type, levels = c('none', 'balanced') )
ggplot(df, aes(y = NES, x= is_TF_perturbed, fill = randomization_type)) + geom_boxplot(alpha = 0.8) + mytheme + scale_fill_manual(values = my_color_palette$clear[3:1])
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_viper_distribution.png', dpi = 300, width = 4, height = 6)



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
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_viper_TFBS_scanning.png', dpi = 300, width = 9, height = 7)


# plot chip-Seq
myplots = plotNES(df, 
                  regulon_filter = c('chipX'), 
                  mycolors = my_color_palette$clear[c(1)])
myplots$NES_boxplot
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_viper_chipX.png', dpi = 300, width = 5, height = 7)



# plot ARACNe
myplots = plotNES(df, 
                  regulon_filter = c('reverse_engenieered'), 
                  TFs_filter = intersectSeveral(lapply(c('> 2 tissues', '> 5 tissues', '> 10 tissues'), function(f) subset( df, regulon_shortname == f)$TF ) ),
                  mycolors = my_color_palette$clear[c(1)])
myplots$NES_boxplot
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_viper_reverse_engenieered.png', dpi = 300, width = 5, height = 7)




# plot curated
myplots = plotNES(df, TFs_filter = NULL,
                  regulon_filter = c("curated_databases"), 
                  mycolors = my_color_palette$clear[c(1)])
myplots$NES_boxplot
plot_grid(myplots$NES_boxplot, myplots$n_barplot, ncol = 1, align = 'v', rel_heights = c(6,4))
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_viper_curated.png', dpi = 300, width = 5, height = 10)




# plot all
myplots = plotNES(df, mycolors = my_color_palette$clear[c(7,1, 4,3,6,5,2)])
myplots$NES_boxplot
plot_grid(myplots$NES_boxplot, myplots$n_barplot, ncol = 1, align = 'v', rel_heights = c(68,32))
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_viper.png', dpi = 300, width = 10, height = 18)





frequent_tfs = names(which(table(unique(df[, c('TF', 'regulon_shortname') ])$TF) > 20))
ggplot(subset(df, TF %in% frequent_tfs & ! is_random), aes(y = NES, x = TF, fill = regulon_type) ) + geom_boxplot() + 
  mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) 
ggsave('data/regulons_QC/pertubations/GEO/relativeActivities_viper_perTF.png', dpi = 300, width = 10, height = 7)

