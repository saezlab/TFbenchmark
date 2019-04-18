rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/regulons_QC/4_plot_benchmark_results/lib_plots_regulonsQC.r')




# Load and merge activity comparisons
load('data/regulons_QC/pertubations/GSE_activity_comparison_results_formated.rdata')



# Counts
nrow(unique(df[ df$is_TF_perturbed, c('expreiment_GEOid', 'TF') ]))
num_TFs = length(unique(df$perturbed_TF[ df$TF == df$perturbed_TF ]))
num_TFs



# Define test and plot features
balanced_accuracy = F
plot_title = 'B2'
outdir = 'results/regulons_QC/GSE/'


# Load samples x GTEx tissue
samples_x_tissue = read.csv('data/TF_target_sources/inferred/GTEx/samples_x_tissue.csv')
samples_x_tissue$GTEx_tissue = gsub(' ', '_', samples_x_tissue$GTEx_tissue)
samples_x_tissue$GTEx_tissue = paste('inferred', samples_x_tissue$GTEx_tissue)





idx = which(df$regulon_group == 'GTEx tissue-specific')
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, legend_position = 'bottom')
P$aucs$sample_x_tissue = samples_x_tissue$n_samples[ match(P$aucs$regulon_dataset, samples_x_tissue$GTEx_tissue) ]
P$aucs
cortest = cor.test(P$aucs$AUC, P$aucs$sample_x_tissue)
label = paste(paste('cor =', round(cortest$estimate, 2)), 
              paste('p-value =', signif(cortest$p.value, 2)), sep = '\n')
P = ggplot(P$aucs, aes(x=sample_x_tissue, y=AUC, size = coverage)) + 
  geom_smooth(method = "lm", show.legend = F, color = 'red4', fill = 'gray90') +
  geom_point() + 
  # annotate('text', label = label, x= Inf, y=Inf, vjust = 1.2, hjust = 1.2) +
  ylab('Area under the PR curve (AUC)') + xlab('#samples') +
  scale_size(name = '', limits = c(1, 40), breaks = c(1, 5, 10, 20, 40)) +
  # scale_color_continuous(name = 'perturbed TFs coverage') +
  ggtitle(plot_title) +
  mytheme 
P
save(P, file = paste('paper/figures/GR_revision_AUC_vs_sample', plot_title,'.rdata', sep = ''))
