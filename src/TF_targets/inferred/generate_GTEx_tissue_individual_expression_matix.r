rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')



# Load expression gene symbols
exp = get(load('~/data/GTEx_expressionatlas/raw/GTex_v6_atlas.genes.voom.nonbatchcor.merged.rdata'))


# Format gene names
ensemblgenes_annot =  get(load(file = '../../databases/ensemblgenes_annot.Rdata'))
# ensemblgenes_annot =  get(load(file = 'data/annotations/hsa_ensembl_proteincoding_20180314.rdata'))
table(rownames(exp) %in% ensemblgenes_annot$ensembl_gene_id)
setdiff(rownames(exp), ensemblgenes_annot$ensembl_gene_id)
ensemblgenes_annot = subset(ensemblgenes_annot, hgnc_symbol != '')
table(rownames(exp) %in% ensemblgenes_annot$ensembl_gene_id)
exp = exp[ rownames(exp) %in% ensemblgenes_annot$ensembl_gene_id, ]
rownames(exp) = ensemblgenes_annot$hgnc_symbol[ match(rownames(exp), ensemblgenes_annot$ensembl_gene_id) ]
any(is.na(rownames(exp) ))


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
tissues_size = list()
for (ti in tissues){
  message(ti)
  e = exp[ , colnames(exp) %in% sample_annotation$Source.Name[ sample_annotation$Comment.histological.type. == ti ] ]
  colnames(e)[1] = paste('gene\t', colnames(e)[1], sep = '')
  e = e[ , ! is.na(colSums(e)) ]
  tissues_size[[ti]] = ncol(e)
  write.table(e, file = paste('data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/', gsub(' ', '_', ti), '.nonbatchcor.txt', sep=''), col.names = T, row.names = T, sep = '\t', quote = F)
}


df = melt(tissues_size)
names(df) = c('n_samples', 'GTEx_tissue')
write.csv(df, file = 'data/TF_target_sources/inferred/GTEx/samples_x_tissue.csv')
ggplot(df, aes(x=GTEx_tissue, y=n_samples, label = n_samples)) + geom_bar(stat= 'identity') + geom_text(hjust = -0.25, size = 3) +
  scale_y_continuous(limits = c(0, 1450) ) +
  coord_flip() + xlab('tissues') + ylab('#samples') +
  mytheme 
ggsave(filename = 'paper/figures/GR_revision_GTEx_samples_per_tissue.png', width = 7, height = 5)
