rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(viper)




# Load TF genesets
load(file = 'data/TF_info/TFrole_genesets.rdata')
TFrole_genesets = unlist(TFrole_genesets, recursive = F)
networks = list()
for( cl in names(TFrole_genesets) ){
  geneset = TFrole_genesets[[cl]]
  networks[[cl]]$tfmode = rep(1, length(geneset))
  networks[[cl]]$likelihood = rep(1, length(geneset))
  names(networks[[cl]]$tfmode) = geneset
}
TFrole_features = unique(sapply(strsplit(names(TFrole_genesets), '\\.'), head, 1))




# LOAD ##################################################################################################
load(file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors_formated.rdata')
df = subset(df, is_TF_perturbed)




all_results = list()
for(dataset  in unique(df$regulon_dataset) ){
  message(dataset)
  
  perturbation_ranking = subset(df, regulon_dataset == dataset & TF == perturbed_TF)
  mySignature = 1- perturbation_ranking$rank_pvalue1tailed
  names(mySignature) = perturbation_ranking$TF
  mySignature = mySignature[order(mySignature, decreasing = T)]
  
  
  sum(names(mySignature) %in% unlist(lapply(unlist(networks, recursive = F), names)))
  
  for ( feature in TFrole_features ) {
    feature_network = networks[ grep(feature, names(networks) ) ]
    mrs = NULL
    try(mrs <- msviper(mySignature, feature_network, verbose = F, minsize = 4, pleiotropy = F, ges.filter = F , adaptive.size = F), silent = T)
    if( is.null(mrs) )
      next()
    results = data.frame(Feature = names(mrs$es$nes),
                         Size = mrs$es$size[ names(mrs$es$nes) ], 
                         NES = mrs$es$nes, 
                         p.value = mrs$es$p.value, 
                         FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
    results = results[ order(results$p.value, decreasing = F) , ]
    results$regulon_evidence = perturbation_ranking$regulon_evidence[1]
    all_results[[feature]][[dataset]] = results
  
  }
}


mdf = melt(all_results, id.vars = names(all_results[[1]][[1]]))
mdf = mdf[ order(mdf$p.value), ]
head(mdf, 20)


mdf$regulon_evidence[ mdf$L2 == 'evidence_1curateddatabases' ] = 'curated_databases'
subset(mdf, p.value < 0.05 &  FDR < 0.1 & L2 == 'evidence_1curateddatabases')
subset(mdf, p.value < 0.05 &  FDR < 0.1 & L2 == 'evidence_2curateddatabases')
subset(mdf, p.value < 0.05 & L2 == 'ReMap_top_500')
subset(mdf, p.value < 0.05 & FDR < 0.1 & L2 %in% 'TFbindingmotif')
subset(mdf, p.value < 0.05 & FDR < 0.1 & L2  == 'inferredGTEx_3' )
subset(mdf, p.value < 0.05 &  FDR < 0.1 & L2 %in% c('A', 'B', 'C', 'D', 'E'))
x = subset(mdf, p.value < 0.05 &  FDR < 0.25 &
             L2 %in% c('evidence_1curateddatabases', 'ReMap_top_500',
                       'inferredGTEx_3', 'TFbindingmotif' ))
unique(x$L1)
x$label = gsub('\\.', ' --- ', x$Feature) %>% gsub('_', ' ', .)
gtitle = 'individual perturbations'
gsubtitle = 'TF annotation - Enrichment'
p3a = ggplot(x, aes(y = NES, x = label, fill = p.value ) ) + geom_bar(stat = 'Identity')  + 
  coord_flip() +
  facet_grid(regulon_evidence~., scales = 'free', space = 'free', as.table = F) +
  scale_fill_gradient(low = brewer.pal(6, 'Blues')[6], high =  brewer.pal(6, 'Blues')[2], name = 'pvalue') +
  labs(title = gtitle,subtitle = gsubtitle) + xlab('')  +
  mytheme + theme(strip.text.y = element_text(angle=0, hjust = 0), 
                  legend.key.width=unit(1.5,"cm"), 
                  plot.subtitle = element_text(face = "italic", hjust = 0.5))
p3a
ggsave('paper/figures/Figure3a.png', width = 12, height = 5, dpi = 300)
save(p3a, file = 'paper/figures/Figure3a.rdata')


