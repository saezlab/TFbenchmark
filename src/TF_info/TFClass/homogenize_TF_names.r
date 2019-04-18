#' This code updates the gene symbols associated to the TFs in the census
#' Gene symbols in the original TFClass file are updated 
#' to match uniprot primary symbols




rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')


TFclass = read.delim('data/TF_census/TFClass/huTF_classification.txt', stringsAsFactors = F)
TFclass$name = sapply(strsplit(TFclass$name, split = "; "), head, 1)
TFcensus = unique(sort(TFclass$name))
length(TFcensus)



# Load uniprot gene symbols (20180314)
prot_annot = load_proteins_annotation(type = 'uniprot')
gene_annot = load_proteins_annotation(type = 'ensembl')


# Check problematic uniprot ids
table(TFclass$uniprot %in% prot_annot$Entry)
setdiff(TFclass$uniprot, prot_annot$Entry)
subset(TFclass, ! uniprot %in%  prot_annot$Entry)


# Check manually the problematic uniprots
mannual_mapping = data.frame(old.Entry = c('P0C7V5', 'Q14591', 'O75505', 'Q6VB85'),
                             entry = c(NA, NA, NA, 'Q8WXT5'), # NA means obsolete accession
                             stringsAsFactors = F)

# Update uniprots
TFclass$uniprot[ TFclass$uniprot == 'Q6VB85' ] = 'Q8WXT5'
subset(TFclass, ! uniprot %in%  prot_annot$Entry)$name %in% prot_annot$Gene.name.primary
subset(TFclass, ! uniprot %in%  prot_annot$Entry)$name %in% gene_annot$hgnc_symbol




# update names
TFclass$name_updated = prot_annot$Gene.name.primary[ match(TFclass$uniprot, prot_annot$Entry) ]
TFclass$name_synonyms = prot_annot$Gene.names[ match(TFclass$uniprot, prot_annot$Entry) ]
subset(TFclass, name != name_updated)


# save census
TFcensus = unique(sort(TFclass$name_updated))
length(TFcensus)
write.table(TFcensus, file = 'data/TF_census/TFClass/huTF_census.txt', row.names = F, col.names = F, quote = F)

