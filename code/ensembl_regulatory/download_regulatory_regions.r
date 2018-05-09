rm(list = ls())
home = '~/Google Drive/projects/TFbenchmark/'
setwd(home)

library(biomaRt)
library(reshape2)


cat('Downloading gene annotation from biomaRt\n')
ensembl = useMart(biomart='ENSEMBL_MART_FUNCGEN', dataset = 'hsapiens_regulatory_feature', host="www.ensembl.org")


regulatory_regions = getBM(attributes=c('regulatory_stable_id', 'chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), mart = ensembl)
regulatory_regions = regulatory_regions[ order(regulatory_regions$chromosome_start) ,]
regulatory_regions = regulatory_regions[ order(regulatory_regions$chromosome_name) ,]
head(regulatory_regions)
dim(regulatory_regions)
length(unique(regulatory_regions$regulatory_stable_id))
write.table(regulatory_regions, file = 'data/genomic_regulatory_regions/ensembl/human_regulatory_features_GRCh38.p10.txt', sep = '\t', quote = F, col.names = T, row.names = F)


regulatory_regions_celltype = list()
for (ch in unique(regulatory_regions$chromosome_name) ){
  message(ch)
  regions = getBM(attributes=c('regulatory_stable_id',  'activity', 'epigenome_name'), mart = ensembl, filters = list(chromosome_name=ch)) 
  regulatory_regions_celltype[[ as.character(ch) ]] =  regions
}
regulatory_regions_celltype_df = melt(regulatory_regions_celltype, id.vars = names(regulatory_regions_celltype[[1]]))
head(regulatory_regions_celltype_df)
dim(regulatory_regions_celltype_df)
length(unique(regulatory_regions_celltype_df$regulatory_stable_id))
write.table(regulatory_regions_celltype_df[, -ncol(regulatory_regions_celltype_df)], file = 'data/genomic_regulatory_regions/ensembl/human_regulatory_features_activity_GRCh38.p10.txt', sep = '\t', quote = F, col.names = T, row.names = F)






# Find active regions
table(regulatory_regions_celltype_df$activity)
regulatory_regions_celltype_df = regulatory_regions_celltype_df[ ! is.na(regulatory_regions_celltype_df$activity) , ]

freq = table(regulatory_regions_celltype_df$regulatory_stable_id)
df = unique(regulatory_regions_celltype_df[,c('regulatory_stable_id',  'activity')])
homogeneus_regions = names(which(table(df$regulatory_stable_id) == 1))
length(homogeneus_regions)
table(df$activity[ df$regulatory_stable_id %in% homogeneus_regions ])
freq[ homogeneus_regions ]
