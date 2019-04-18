rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(fgsea)




# Load TF genesets
load(file = 'data/TF_info/TFrole_genesets.rdata')
TFrole_features = unique(sapply(strsplit(names(TFrole_genesets), '\\.'), head, 1))



# LOAD ##################################################################################################
load(file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors_formated.rdata')
df = subset(df, is_TF_perturbed)




all_results = list()
df$regulon_dataset[ df$regulon_dataset %in% c('jaspar2018_n500', 'hocomocov11_n500') ] = 'TFbindingmotif'
datasets2study = unique(df$regulon_dataset) %>% grep('evidence2', ., invert = T, value = T) %>% grep('omnipath', ., invert = T, value = T) 
for(dataset in datasets2study){
  message(dataset)
  
  perturbation_ranking = subset(df, regulon_dataset == dataset & TF == perturbed_TF)
  mySignature = perturbation_ranking$rank_nes
  names(mySignature) = perturbation_ranking$TF
  mySignature = mySignature[order(mySignature, decreasing = T)]
  
  for ( feature in TFrole_features ) {
    feature_network = TFrole_genesets[[feature]]
    
    results = fgsea(pathways = feature_network,  stats = mySignature[ names(mySignature) %in% unlist(feature_network)  ], minSize=3, 
                    maxSize=2000, nperm=1000)
    results = results[ order(results$padj, decreasing = F) , ]
    results$regulon_evidence = perturbation_ranking$regulon_evidence[1]
    all_results[[feature]][[dataset]] = results
  }
}



mdf = melt(all_results, id.vars = names(all_results[[1]][[1]]))
mdf = mdf[ order(mdf$pval), which(names(mdf) != 'leadingEdge')  ]
mdf$regulon_evidence[ mdf$L2 == 'evidence_1curateddatabases' ] = 'curated_databases'
subset(mdf, pval < 0.05 &  padj < 0.1 & L2 == 'evidence_1curateddatabases')
subset(mdf, pval < 0.05 &  padj < 0.1 & L2 == 'evidence_2curateddatabases')
subset(mdf, pval < 0.05 & L2 == 'ReMap_top_500')
subset(mdf, pval < 0.05 & padj < 0.1 & L2 %in% 'TFbindingmotif')
subset(mdf, pval < 0.05 & padj < 0.1 & L2  == 'inferredGTEx_3' )
subset(mdf, pval < 0.05 &  padj < 0.1 & L2 %in% c('A', 'B', 'C', 'D', 'E'))
x = subset(mdf, pval < 0.05 &  padj < 0.1 &
             L2 %in% c('evidence_1curateddatabases', 'ReMap_top_500',
                       'inferredGTEx_3', 'TFbindingmotif' ))
unique(x$L1)
x$label = paste(x$L1, x$pathway, sep = ' --- ') %>% gsub('_', ' ', .)
gtitle = 'individual perturbations'
gsubtitle = 'TF annotation - Enrichment'
x
p3a = ggplot(x, aes(y = NES, x = label, fill = pval ) ) + geom_bar(stat = 'Identity')  + 
  coord_flip() +
  facet_grid(regulon_evidence~., scales = 'free', space = 'free', as.table = F) +
  scale_fill_gradient(low = brewer.pal(6, 'Blues')[6], high =  brewer.pal(6, 'Blues')[2], name = 'pvalue') +
  scale_y_continuous(limits = c(-2, 2))  +
  labs(title = gtitle,subtitle = gsubtitle) + xlab('')  +
  mytheme + theme(strip.text.y = element_text(angle=0, hjust = 0), 
                  legend.key.width=unit(1.5,"cm"), 
                  plot.subtitle = element_text(face = "italic", hjust = 0.5))
p3a
ggsave('paper/figures/Figure3/Figure3a.png', width = 12, height = 5, dpi = 300)
save(p3a, file = 'paper/figures/Figure3/Figure3a.rdata')
