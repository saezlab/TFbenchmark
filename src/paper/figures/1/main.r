rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')



# Load network
network = read.csv(file = 'data/TF_target_sources/omnipath_scores/database_20180915.csv', stringsAsFactors = F)
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
TFlist$only_inferred_GTEx = as.character(m[ which(rowSums(m[,-1]) == 1  & m[,'inferred GTEx'] == 1), "Identifier" ])
TFlist$only_chip = as.character(m[ which(rowSums(m[,-1]) == 1  & m[,'ChIPseq'] == 1), "Identifier" ])
TFlist$only_curated = as.character(m[ which(rowSums(m[,-1]) == 1  & m[,'curated'] == 1), "Identifier" ])
TFlist$only_tbfs = as.character(m[ which(rowSums(m[,-1]) == 1  & m[,'TF binding motif'] == 1), "Identifier" ])
TFlist$inferred = as.character(m[ which(m[,'inferred GTEx'] == 1), "Identifier" ])
TFlist$curated = as.character(m[ which(m[,'curated'] == 1), "Identifier" ])
TFlist$ChIPseq = as.character(m[ which(m[,'ChIPseq'] == 1), "Identifier" ])
TFlist$TFBS = as.character(m[ which(m[,'TF binding motif'] == 1), "Identifier" ])
TFlist$at_least_3_evidences = as.character(m[ which(rowSums(m[,-1]) > 2), "Identifier" ])
TFlist$at_least_2_evidences = as.character(m[ which(rowSums(m[,-1]) > 1), "Identifier" ])
TFlist$at_least_4_evidences = as.character(m[ which(rowSums(m[,-1]) > 3), "Identifier" ])



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
load(file = 'data/annotations/KEGGpathways_SLAPE_MSigDB.rdata')
TFrole_genesets$regulatory_effect$unknown = setdiff(unlist(TFrole_genesets$TF_class), unlist(TFrole_genesets$regulatory_effect))
TFrole_genesets$tissue_of_expression$intermediate = setdiff(unlist(TFrole_genesets$TF_class), unlist(TFrole_genesets$tissue_of_expression))
source('code/paper/figures/1/lib_enrichment.r')
TFrole_features = unique(sapply(strsplit(names(TFrole_genesets), '\\.'), head, 1))


# Figure 1C
re3 = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$at_least_3_evidences)
re3
plot_enrichment(re3,feature = 'regulatory_effect') + ggtitle('TFs covered by > 2 evidences')
plot_enrichment(re3) + ggtitle('TFs covered by > 2 evidences')
ggsave(filename = 'paper/figures/Figure1/Figure1c.png',  dpi=300, width = 7, height = 3.5)



# Figure S1
# A
re = list()
re$tissues = list()
TFrole_genesets$tissue_of_expression$expressed_in_most_tissues = TFrole_genesets$tissue_of_expression$`no_tissue-specific`
TFrole_genesets$tissue_of_expression$`no_tissue-specific` = NULL
re$tissues$at_least_3_evidences = analyse_genesets(geneset =  TFrole_genesets$tissue_of_expression, genes = TFlist$at_least_3_evidences)
re$tissues$at_least_2_evidences = analyse_genesets(geneset =  TFrole_genesets$tissue_of_expression, genes = TFlist$at_least_2_evidences)
re$tissues$at_least_4_evidences = analyse_genesets(geneset =  TFrole_genesets$tissue_of_expression, genes = TFlist$at_least_4_evidences)
re$tissues$only_inferred_GTEx = analyse_genesets(geneset =  TFrole_genesets$tissue_of_expression, genes = TFlist$only_inferred_GTEx)
re$tissues$ChIPseq = analyse_genesets(geneset =  TFrole_genesets$tissue_of_expression, genes = TFlist$ChIPseq)
re$tissues$curated = analyse_genesets(geneset =  TFrole_genesets$tissue_of_expression, genes = TFlist$curated)
re$tissues$TFBS = analyse_genesets(geneset =  TFrole_genesets$tissue_of_expression, genes = TFlist$TFBS)



# B
re$kegg = list()
re$kegg$at_least_3_evidences = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$at_least_3_evidences)
re$kegg$at_least_2_evidences = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$at_least_2_evidences)
re$kegg$at_least_4_evidences = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$at_least_4_evidences)
re$kegg$only_inferred_GTEx = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$only_inferred_GTEx)
re$kegg$ChIPseq = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$ChIPseq)
re$kegg$curated = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$curated)
re$kegg$TFBS = analyse_genesets(geneset =  KEGG_PATH$HGNC_SYMBOL, genes = TFlist$TFBS)




# Figure S1
# C
re$domains = list()
re$domains$at_least_3_evidences = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$at_least_3_evidences)
re$domains$at_least_2_evidences = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$at_least_2_evidences)
re$domains$at_least_4_evidences = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$at_least_4_evidences)
re$domains$only_inferred_GTEx = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$only_inferred_GTEx)
re$domains$ChIPseq = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$ChIPseq)
re$domains$curated = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$curated)
re$domains$TFBS = lapply(TFrole_genesets, analyse_genesets, genes = TFlist$TFBS)


plot_enrichment_grid(re$tissues)
ggsave(filename = 'paper/figures/supplementary/FigureS1a.png', width = 15, height = 2.3)
plot_enrichment_grid(re$kegg)
ggsave(filename = 'paper/figures/supplementary/FigureS1b.png', width = 15, height = 10)
plot_enrichment_grid(lapply(re$domains, function(x) x$TF_class ))
ggsave(filename = 'paper/figures/supplementary/FigureS1c.png', width = 15, height = 4)


