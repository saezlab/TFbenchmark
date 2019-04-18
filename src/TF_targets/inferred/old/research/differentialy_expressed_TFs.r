rm(list = ls())
home = '~/Google Drive/projects/TFbenchmark/'
setwd(home)




# Load expression
exp = get(load('~/data/GTEx_expressionatlas/raw/GTex_v6_atlas.genes.voom.batchcor.merged.rdata'))
fpkms = get(load('~/data/GTEx_expressionatlas/raw/GTex_v6_atlas.genes.fpkms.merged.rdata'))



# Load samples annotation
load('~/data/GTEx_expressionatlas/annot/E-MTAB-5214.sdrf.rdata')



# Load TFs
TFs = read.delim(file = 'data/TF_census/vaquerizas/TF_census.txt', header = F, stringsAsFactors = F)[,1]



# Load gene annot and convert gene names
exp =  exp[ rownames(exp) %in% ensemblgenes_annot$ensembl_gene_id, ]
rownames(exp) = ensemblgenes_annot$hgnc_symbol[ match(rownames(exp), ensemblgenes_annot$ensembl_gene_id) ]
exp = exp[ ! duplicated(rownames(exp)), ]
fpkms =  fpkms[ rownames(fpkms) %in% ensemblgenes_annot$ensembl_gene_id, ]
rownames(fpkms) = ensemblgenes_annot$hgnc_symbol[ match(rownames(fpkms), ensemblgenes_annot$ensembl_gene_id) ]
fpkms = fpkms[ ! duplicated(rownames(fpkms)), ]



# Filter exp for TFs
table(rownames(exp) %in% TFs)
table(TFs %in% rownames(exp))
exp = exp[ rownames(exp) %in% TFs, ]
fpkms = fpkms[ rownames(fpkms) %in% TFs, ]




# plot profiles
tissues = sample_annotation$Comment.histological.type.[ match(colnames(exp), sample_annotation$Source.Name)] 
idx = sample(1:ncol(exp), 1000)
tissue_cluster = tissues[idx] 
library(Rtsne)
library(RColorBrewer)
colors = sample(c(brewer.pal(n = 10, name = "Paired"), brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 4, name = "Set1")[3:4])[1:length(unique(tissue_cluster))])
names(colors) = unique(tissue_cluster)
tsne = Rtsne(scale(t(exp[, idx ])), dims = 2, perplexity=50, verbose=TRUE, max_iter = 1000)
plot(tsne$Y, t='n', main="tSNE")
text(tsne$Y, labels=tissue_cluster, col=colors[tissue_cluster], cex = .7)




tissues = sample_annotation$Comment.histological.type.[ sample_annotation$Source.Name %in% colnames(exp) ]



##########################################################################################################################################################################################################################################
# differential expression per tissue type
##########################################################################################################################################################################################################################################
ggplot( sample_annotation[ sample_annotation$Source.Name %in% colnames(exp) , ], aes (x = Comment.histological.type.) ) + geom_bar() + theme_classic(18) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + xlab('tissue')
table(tissues)
results = list()
for (ti in sort(unique(tissues))){
  message(ti)
  design = model.matrix(~factor((sample_annotation$Comment.histological.type.[ match(colnames(exp),  sample_annotation$Source.Name) ] == ti) + 0))
  ## fit the same linear model now to the enrichment scores
  fit = lmFit(exp, design)
  ## estimate moderated t-statistics
  fit = eBayes(fit, robust=TRUE, trend=TRUE)
  ## set1 is differentially expressed
  results[[ti]] = topTable(fit, coef = 2, p.value = 1, number = nrow(exp)) 
  results[[ti]]$TF = rownames(results[[ti]])
}
df = melt(results, id.vars = names(results[[1]]))
df$adj.P.Val_global = p.adjust(df$P.Value, method = 'fdr')
save(df, file = 'data/TF_target_sources/reverse_engenieered_networks/tf_differential_expression.rdata')


# pvalue heatmap
df$value = df$adj.P.Val +1E-322
M = -log10(acast(df, formula = L1~TF))
df$value = df$logFC
M = M * sign(acast(df, formula = L1~TF))
pheatmap::pheatmap(M)


# Counts
df$veredict = 'neutral'
df$veredict[ df$adj.P.Val < 0.01 & df$logFC < 0 ] = '-'
df$veredict[ df$adj.P.Val < 0.01 & df$logFC > 0 ] = '+'
ggplot(df, aes(x= veredict, fill = veredict)) + geom_bar() + 
  theme_bw(18) + xlab('TF status in tissue') + 
  scale_fill_manual(values = c('steelblue3', 'coral3', 'grey')) + theme(legend.position = 'none') 


ddf = melt(ddply(df, 'TF', function(x) table(x$veredict) ), id.vars = 'TF')
ddf$value = as.factor(ddf$value)
ggplot(ddf[ ddf$variable != 'neutral', ], aes(x= value, fill = variable) ) + geom_bar(position = 'dodge') + 
  scale_fill_manual(values = c('steelblue3', 'coral3'), name = 'TF status in tissue') + 
  theme_bw(18) + ylab('number of TFs') + xlab('number of tissues') + ggtitle('TF expression in 30 tissues')
##########################################################################################################################################################################################################################################







##########################################################################################################################################################################################################################################
# 25th > 5 fpkms
##########################################################################################################################################################################################################################################
results = list()
for (ti in sort(unique(tissues))){
  message(ti)
  results[[ti]]$active =  names(which(apply(fpkms[, colnames(fpkms) %in% sample_annotation$Source.Name[ sample_annotation$Comment.histological.type. == ti  ] ], 1, quantile, 0.25) > 1))
  results[[ti]]$inactive = names(which(apply(fpkms[, colnames(fpkms) %in% sample_annotation$Source.Name[ sample_annotation$Comment.histological.type. == ti  ] ], 1, quantile, 0.75) < 1))
}
df = melt(results)
ddf = as.data.frame(rbind( cbind(table(df$value[ df$L2 == 'active']), 'active'), cbind(table(df$value[ df$L2 == 'inactive']), 'inactive')), stringsAsFactors = F)
ddf$V1 = as.factor(as.integer(ddf$V1))
ggplot(ddf, aes(x= V1, fill = V2) ) + geom_bar(position = 'dodge') + 
  scale_fill_manual(values = c('steelblue3', 'coral3'), name = 'TF status in tissue') + 
  theme_bw(18) + ylab('number of TFs') + xlab('number of tissues') + ggtitle('TF expression in 30 tissues')
save(df, file = 'data/TF_target_sources/reverse_engenieered_networks/tf_tissues_call_percentile.rdata')
