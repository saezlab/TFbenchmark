rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(viper)



TFs = load_TFs_census()


# Load network_files
EXP = get(load('~/data/GTEx_expressionatlas/raw/GTex_v6_atlas.genes.voom.batchcor.merged.rdata'))


# Get coexpression
coEXP = cor(t(EXP[, sample(1:ncol(EXP), 1000)]))
coEXP[ coEXP < 0 ] = abs( min(coEXP) )
diag(coEXP) = min(coEXP)
relcoEXP = coEXP / apply(coEXP, 1, max)


# Get gene names
gene_symbol_annot = get(load('../../databases/ensemblgenes_annot.Rdata'))
genes = gene_symbol_annot$hgnc_symbol[ match(colnames(relcoEXP), gene_symbol_annot$ensembl_gene_id) ]
colnames(relcoEXP) = genes
rownames(relcoEXP) = genes


# Define background targets
targets = colnames(relcoEXP)
sizes = 1300:3
n_randoms = 1:100


# Build random network
for (si in sizes){
  message('size ', si)
  randoms = list()
  for (ri in n_randoms ){
    randomGene = sample(targets, 1)
    randoms[[paste('r', ri, '_', si, sep = '')]] = sample(targets, replace = F, prob = relcoEXP[randomGene, ], size = si)
  }
  networks_random = melt(randoms)[, 2:1]
  names(networks_random) = c("TF", "target")
  # write.table(networks_random, file= paste('data/TF_target_sources/RANDOM_NETWORKS_coexpressionPriors/', si, '_network.sif', sep = ''), col.names = T, row.names = F, quote = F, sep = '\t')
  viper_regulon = list()
  for (tf in unique(networks_random$TF) ){
    tf_targets = subset(networks_random, TF == tf)$target
    # Define TF mode
    tf_targets_sign = rep(1, length(tf_targets))
    # Define regulation weight
    tf_targets_weight = rep(1, length(tf_targets))
    # Build VIPER regulon
    viper_regulon[[tf]]$tfmode = tf_targets_sign # NOTE: ttmode has to be the first element of the list. Otherwise msviper will fail.
    names(viper_regulon[[tf]]$tfmode) = tf_targets
    viper_regulon[[tf]]$likelihood = tf_targets_weight
  }
  save(viper_regulon, file = paste('data/TF_target_sources/RANDOM_NETWORKS_coexpressionPriors/', si, '_network_viperRegulon.rdata', sep = '') )
}



