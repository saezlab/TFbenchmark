rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/regulons_QC/5_TFproperties_enrichment/gsealm/lib_gsealm.r')




# Load TF genesets
load(file = 'data/TF_info/TFrole_genesets.rdata')
TFrole_features = unique(sapply(strsplit(names(TFrole_genesets), '\\.'), head, 1))




# Load and merge activity comparisons
load(file = 'data/regulons_QC/cell_lines/activity_comparison_results_formated.rdata' )
df = subset(df, is_TF_perturbed)
df$perturbed_TF = df$TF


# Enrichment test
df$regulon_dataset[ df$regulon_dataset %in% c('jaspar2018_n200', 'hocomocov11_n200')  ] = 'TFBS'
df$regulon_dataset[ df$regulon_dataset %in% c('jaspar2018_n100', 'hocomocov11_n100', 'jaspar2018_n200', 'hocomocov11_n200','jaspar2018_n500', 'hocomocov11_n500', 'jaspar2018_n1000', 'hocomocov11_n1000')  ] = 'TFBS'
df$regulon_dataset[ df$regulon_dataset %in% c('evidence_1curateddatabases', 'evidence_2curateddatabases', 'evidence_3curateddatabases')  ] = 'curated databases'
df$regulon_dataset[ df$regulon_dataset %in% c('inferredGTEx_2', 'inferredGTEx_3', 'inferredGTEx_5', 'inferredGTEx_10')  ] = 'consensus inferred GTEx'
df$regulon_dataset[ df$regulon_dataset %in% c('ReMap_top_100', 'ReMap_top_200', 'ReMap_top_500', 'ReMap_top_1000')  ] = 'ChIP-seq'
df$signature = df$NES #df$rank_nes - 0.5 # 
mdf = TF_properties_enrichment(df, genesets = TFrole_genesets)

# save
save(mdf, file = 'paper/figures/Figure3/enrichment_b3.rdata')

