rm(list = ls())
source('code/utils.r')



# Load network
network = read.csv(file = 'data/TF_target_sources/omnipath_scores/database.csv', stringsAsFactors = F)
names(network)[5:8] = c('curated', 'ChIPseq', 'TF binding motif', 'inferred GTEx')
nrow(unique(network[, 1:2]))


# retrieve TF
x = lapply(names(network)[5:8], function(evidence){
  unique(network[  network[,evidence] , ]$TF)
})
names(x) = gsub('is_evidence_', '', names(network)[5:8])
sapply(x, length)
df = melt(x)
names(df) = c('TF', 'dataset')
df$value = 1
m = dcast(df, TF~dataset, fill = 0)
names(m)[1] = 'Identifier'
png('paper/figures/Figure1/Figure1b.png', res = 300, width = 1600, height = 1200)
upset(m, sets = colnames(m)[-1], main.bar.color = 'gray20', sets.bar.color = my_color_palette$EMBL[4], 
      empty.intersections = F, set_size.angles = 90, number.angles = 25, 
      order.by = "freq", point.size = 2.5, line.size = 0.5, mb.ratio = c(.6, .4), text.scale = c(1.5, 1.2, 1.5, 1.2, 1.5, 1),
      mainbar.y.label = 'shared TFs', sets.x.label = 'total TFs') # 9000x3100
dev.off()
TFlist = list()
TFlist$only_inferred = as.character(m[ which(rowSums(m[,-1]) == 1  & m[,'inferred GTEx'] == 1), "Identifier" ])
TFlist$inferred = as.character(m[ which(m[,'inferred GTEx'] == 1), "Identifier" ])
TFlist$curated = as.character(m[ which(m[,'curated'] == 1), "Identifier" ])
TFlist$chip = as.character(m[ which(m[,'ChIPseq'] == 1), "Identifier" ])
TFlist$TFBS = as.character(m[ which(m[,'TF binding motif'] == 1), "Identifier" ])
TFlist$atleast3 = as.character(m[ which(rowSums(m[,-1]) > 2), "Identifier" ])



# retrieve TFTG
network$TFtarget = paste(network$TF, network$target)
x = lapply(names(network)[5:8], function(evidence){
  unique(network[  network[,evidence] , ]$TFtarget)
})
names(x) = gsub('is_evidence_', '', names(network)[5:8])
sapply(x, length)
df = melt(x)
names(df) = c('TFtarget', 'dataset')
df$value = 1
m = dcast(df, TFtarget~dataset, fill = 0)
names(m)[1] = 'Identifier'
png('paper/figures/Figure1/Figure1d.png', res = 300, width = 2000, height = 1200)
upset(m, sets = colnames(m)[-c(1)], main.bar.color = 'gray20', sets.bar.color = my_color_palette$EMBL[2], scale.intersections = 'log2',
      empty.intersections = F, set_size.angles = 90, number.angles = 21, 
      order.by = "freq", point.size = 2.5, line.size = 0.5, mb.ratio = c(.6, .4), text.scale = c(1.5, 1.2, 1.5, 1.2, 1.5, 1),
      mainbar.y.label = 'shared TF-TG', sets.x.label = 'total TF-TG') # 9000x3100
dev.off()





# enrichment
load(file = 'data/TF_info/TFrole_genesets.rdata')
load('../pathway_activities/DoRothEA/DATA/pathways/SLAPE.MSigDB_KEGG_hugoUpdated.rdata')
TFrole_genesets$regulatory_effect$unknown = setdiff(unlist(TFrole_genesets$TF_class), unlist(TFrole_genesets$regulatory_effect))
TFrole_genesets$tissue_of_expression$intermediate = setdiff(unlist(TFrole_genesets$TF_class), unlist(TFrole_genesets$tissue_of_expression))
source('code/paper_figures/1/lib_enrichment.r')
TFrole_features = unique(sapply(strsplit(names(TFrole_genesets), '\\.'), head, 1))


re3 = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$atleast3)
re3
plot_enrichment(re3,feature = 'regulatory_effect') + ggtitle('TFs covered by > 2 evidences')
plot_enrichment(re3) + ggtitle('TFs covered by > 2 evidences')
ggsave(filename = 'paper/figures/Figure1/Figure1c.png',  dpi=300, width = 7, height = 3.5)

re4 = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$atleast3)
plot_enrichment(re4, feature = NULL) + ggtitle('TFs covered by > 2 evidences')
re4
ggsave(filename = 'paper/figures/supplementary/FigureS1.png',  dpi=300, width = 10, height = 10)


re1 = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$only_inferred)
plot_enrichment(re1) + ggtitle('only_inferred')
plot_enrichment(re1, feature = 'regulatory_effect') + ggtitle('only_inferred')
re1 = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$only_inferred)
plot_enrichment(re1, feature = NULL) + ggtitle('only_inferred')

re2 = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$chip)
plot_enrichment(re2, feature = 'regulatory_effect') + ggtitle('chip')
plot_enrichment(re2) + ggtitle('chip')
re2 = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$chip)
plot_enrichment(re2, feature = NULL) + ggtitle('chip')


re2 = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$curated)
plot_enrichment(re2, feature = 'regulatory_effect') + ggtitle('curated')
plot_enrichment(re2) + ggtitle('curated')
re2 = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$curated)
plot_enrichment(re2, feature = NULL) + ggtitle('curated')


re2 = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$TFBS)
plot_enrichment(re2, feature = 'regulatory_effect') + ggtitle('TFBS')
plot_enrichment(re2) + ggtitle('TFBS')
re2 = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$TFBS)
plot_enrichment(re2, feature = NULL) + ggtitle('TFBS')



re2 = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$inferred)
plot_enrichment(re2, feature = 'regulatory_effect') + ggtitle('inferred')
plot_enrichment(re2) + ggtitle('inferred')
re2 = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$inferred)
plot_enrichment(re2, feature = NULL) + ggtitle('inferred')

