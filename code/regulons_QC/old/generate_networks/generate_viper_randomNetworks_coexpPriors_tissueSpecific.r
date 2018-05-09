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
coEXP = abs(coEXP)
diag(coEXP) = min(coEXP)
relcoEXP = coEXP / apply(coEXP, 1, max)


# Get gene names
gene_symbol_annot = get(load('../../databases/ensemblgenes_annot.Rdata'))
genes = gene_symbol_annot$hgnc_symbol[ match(colnames(relcoEXP), gene_symbol_annot$ensembl_gene_id) ]
colnames(relcoEXP) = genes
rownames(relcoEXP) = genes
targets = colnames(relcoEXP)



# Load REAL GTEx tissue-specific networks 
network_files = list.files('data/TF_target_sources', recursive = T, pattern = 'viperRegulon.rdata')
network_files = grep("aracne", network_files, value = T)

randoms = 1:100

for(net in network_files){
  regulon = get(load(paste('data/TF_target_sources/', net, sep = '')))
  random_regulon = list()
  for( reg in names(regulon) ){
    message(reg)
    si = length(regulon[[reg]]$tfmode)
    for(i in randoms){
      random_regulon[[paste(reg, i, sep = ' - r')]] = regulon[[reg]]
      randomGene = sample(targets, 1)
      names(random_regulon[[paste(reg, i, sep = ' - r')]]$tfmode) = sample(targets, replace = F, prob = relcoEXP[randomGene, ], size = si)
    }
  }
  filename = paste('data/TF_target_sources/RANDOM_NETWORKS_coexpressionPriors_tissueSpecific/', tail(unlist(strsplit(net, split = '/')), 1), sep = '')
  save(random_regulon, file = filename )
}


