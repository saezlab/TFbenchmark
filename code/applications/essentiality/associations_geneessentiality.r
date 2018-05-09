rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
source('code/applications/associations_geneessentiality_lib.r')



#---------------------- TF activities
message('Loading activities')
activities = get(load('data/regulons_QC/cell_lines/activities/omnipath_scores.TOP.network.rdata'))
rownames(activities) = sapply(strsplit(rownames(activities), split = ' - '), head, 1)

# Remove low confidence regulons
activities = activities[ grep('_E$', rownames(activities), invert = T) , ]
activities = activities[ grep('_D$', rownames(activities), invert = T) , ]
# activities = activities[ grep('_C$', rownames(activities), invert = T) , ]
# activities = activities[ grep('_B$', rownames(activities), invert = T) , ]
# activities = activities[ grep('_A$', rownames(activities), invert = T) , ]
dim(activities)
#----------------------



#---------------------- Essentiality profiles
message('Loading essentiality profiles')
essentiality = load_ACHILLES()
essentiality = t(scale(t(essentiality)))

# invert DEMETER zscore (the larger the more essential)
essentiality = (-1) * essentiality

# Remove core essential genes from the essentiality features
core_essential = sapply(list.files(path = 'data/applications/essentiality/core_essential_genes/', full.names = T), function(x) get(load(x))) %>% unlist(.) %>% unique(.)
essentiality = essentiality[ ! rownames(essentiality) %in% core_essential, ]
dim(essentiality)

# Remove nondruggable genes
druggable_genes = read.csv('data/applications/druggability/drugable_genome_ Finan2017_Science_aag1166_Table S1.csv', header = T, stringsAsFactors = F)
essentiality = essentiality[   rownames(essentiality) %in% druggable_genes$hgnc_names, ]
dim(essentiality)

# Remove TFs
TFs = load_TFs_census()
essentiality = essentiality[ ! rownames(essentiality) %in% TFs , ]
dim(essentiality)

# # Filter essentiality scores
# max_zscores = apply(essentiality, 1, function(x) max(abs(x), na.rm = T) )
# essentiality = essentiality[ max_zscores > 3, ]
# dim(essentiality)
#----------------------


# # 
# # #---------------------- Pancancer ASSOS
# associations_df = analyse_TF_essentiality_associations(activities, essentiality, tissue = 'pancancer')
# plot_volcano(associations_df, th = 0.01, main = 'TF-gene essentiality interactions\n(pancancer)')
# ggsave(filename = 'results/applications/gene_essentiality_assos/pancancer.png', dpi = 300, width = 6, height = 6)
# save(associations_df, file = 'results/applications/gene_essentiality_assos/pancancer.rdata')
# # #----------------------



#---------------------- Cancer-type ASSOS
tissues_of_interest = names(which(table(covariates.annot$gdsc_desc_2[ covariates.annot$CosmicID %in% colnames(essentiality) ])  > 10) )
for (tissue in rev(tissues_of_interest)){
  message(tissue)
  associations_df = analyse_TF_essentiality_associations(activities, essentiality, tissue = tissue)
  plot_volcano(associations_df)
  ggsave(filename = paste('results/applications/gene_essentiality_assos/', tissue, '.png', sep = ''), dpi = 300, width = 10, height = 10)
  save(associations_df, file = paste('results/applications/gene_essentiality_assos/', tissue, '.rdata', sep = ''))
}
#----------------------