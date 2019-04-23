#' This code generates an annotation file for indicating which TF is expressed in wich GTEx tissue
#' GTEx gene expression data (fpkm) 
#' is transformed into a binary expression matrix
#'    0: no expression
#'    1: expression
#' output file: data/TF_info/expression/GTEx_TF_per_tissue.txt
#' and reproduces the plots from Vaquerizas et al 2009





rm(list = ls())
home = '~/data/GTEx_expressionatlas/'
setwd(home)
source('/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/code/lib/utils.r')


# Load expression gene symbols
exp = get(load('raw/GTex_v6_atlas.genes.fpkms.merged.rdata'))
ensemblgenes_annot =  get(load(file = '/Volumes/GoogleDrive/My Drive/databases/ensemblgenes_annot.Rdata'))
rownames(exp) = ensemblgenes_annot$hgnc_symbol[ match(rownames(exp), ensemblgenes_annot$ensembl_gene_id) ]



# Load samples annotation
load('annot/E-MTAB-5214.sdrf.rdata')


# Define tissues
tissues = sample_annotation$Comment.histological.type.[ sample_annotation$Source.Name %in% colnames(exp) ]
tissues = sort(unique(tissues))


# Define Tfs
TFcesus = load_TFs_census()
tf_idx = rownames(exp) %in% TFcesus


tissue_TFs = list()
exp_dist = list()
for (ti in tissues){
  message(ti)
  e = exp[ tf_idx , colnames(exp) %in% sample_annotation$Source.Name[ sample_annotation$Comment.histological.type. == ti ] ]
  tissue_TFs[[ti]] = names(which(apply(e, 1, mean) >= 2)) # or 1st quantile > 1fpkm
  
  exp_dist[[ti]]$nonTF = apply(exp[ ! tf_idx , colnames(exp) %in% sample_annotation$Source.Name[ sample_annotation$Comment.histological.type. == ti ] ], 1, median)
  exp_dist[[ti]]$TF = apply(e, 1, median)
}
df = melt(tissue_TFs)
df_dist = melt(exp_dist)
df_dist = subset(df_dist, value != 0 )
df_dist$value = log2(df_dist$value)


# Save annotation
write.table(df, file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_info/expression/GTEx_TF_per_tissue.txt', col.names = F, row.names = F, sep = '\t', quote = F)



# PLOTS as in Vaquerizas
p1 = ggplot(df_dist, aes(x = L1, y = value, fill = L2) ) + geom_boxplot(outlier.size = 0) +
  scale_fill_manual(values = my_color_palette$shiny[c(6,4)], name = 'class') + 
  mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylab('log2 fpkms') + xlab('tissues')


df_freq = as.data.frame(table(df$L1), stringsAsFactors = F)
df_freq = df_freq[ order(df_freq$Freq, decreasing = T), ]
df_freq$Var1 = factor(df_freq$Var1, levels = df_freq$Var1 )
p2 = ggplot(df_freq, aes(x= Var1, y = Freq)) + geom_bar(stat = 'identity', fill = my_color_palette$clear[1]) + 
  geom_hline(yintercept = mean(df_freq$Freq), color = 'red' ) +
  mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylab('# TFs') + xlab('tissue')


df_freq = as.data.frame(table(table(df$value)), stringsAsFactors = F)
df_freq$Var1 = factor(df_freq$Var1, levels = df_freq$Var1 )
p3 = ggplot(df_freq, aes(x= Var1, y = Freq)) + geom_bar(stat = 'identity', fill = my_color_palette$clear[1]) + 
  geom_hline(yintercept = mean(df_freq$Freq), color = 'red' ) +
  mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylab('# TF') + xlab('# tissues')

P = plot_grid(p1, p2, p3, align = 'v', ncol = 1)
save_plot(P, filename = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_info/expression/GTEx_TF_expression_distribution.png', base_width = 6, base_height = 11)
