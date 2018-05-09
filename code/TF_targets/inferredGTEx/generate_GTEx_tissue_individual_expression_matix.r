rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)



# Load expression gene symbols
exp = get(load('~/data/GTEx_expressionatlas/raw/GTex_v6_atlas.genes.voom.batchcor.merged.rdata'))


# Format gene names
ensemblgenes_annot =  get(load(file = 'data/annotations/hsa_ensembl_proteincoding_20180314.rdata'))
rownames(exp) = ensemblgenes_annot$hgnc_symbol[ match(rownames(exp), ensemblgenes_annot$ensembl_gene_id) ]



# Check if TFs are there
TF_census = read.delim('data/TF_census/TFClass/huTF_census.txt', stringsAsFactors = F)[,1]
table(TF_census %in% rownames(exp))
table(TF_census %in% ensemblgenes_annot$hgnc_symbol)
setdiff(TF_census, ensemblgenes_annot$hgnc_symbol)


# Load samples annotation
load('~/data/GTEx_expressionatlas/annot/E-MTAB-5214.sdrf.rdata')



tissues = sample_annotation$Comment.histological.type.[ sample_annotation$Source.Name %in% colnames(exp) ]
table(tissues)
dim(exp)
tissues = sort(unique(tissues))
for (ti in tissues){
  message(ti)
  e = exp[ , colnames(exp) %in% sample_annotation$Source.Name[ sample_annotation$Comment.histological.type. == ti ] ]
  colnames(e)[1] = paste('gene\t', colnames(e)[1], sep = '')
  e = e[ , ! is.na(colSums(e)) ]
  write.table(e, file = paste('data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/', gsub(' ', '_', ti), '.txt', sep=''), col.names = T, row.names = T, sep = '\t', quote = F)
}
