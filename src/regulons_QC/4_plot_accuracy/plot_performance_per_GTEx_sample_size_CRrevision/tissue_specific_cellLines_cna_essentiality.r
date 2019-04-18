rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/regulons_QC/4_plot_benchmark_results/lib_plots_regulonsQC.r')




# Load and merge activity comparisons
load(file = 'data/regulons_QC/cell_lines/activity_comparison_results_formated.rdata' )
df = subset(df, ! regulon_dataset %in% c('all_weighted', 'cnio', 'inferredGTEx_dynamic', 'evidence4', 'TOP') )
df$rank_nes = 1 - df$rank_nes
df$NES = df$NES * -(1)



# Counts
nrow(unique(df[ df$is_TF_active, c('experiment', 'TF') ]))
num_TFs = length(unique(df$TF[  df$is_TF_active ]))
num_TFs


# Define test and plot features
balanced_accuracy = F
plot_title = 'B3'
outdir = 'results/regulons_QC/cell_lines/'





# Load samples x GTEx tissue
samples_x_tissue = read.csv('data/TF_target_sources/inferred/GTEx/samples_x_tissue.csv')
samples_x_tissue$GTEx_tissue = gsub(' ', '_', samples_x_tissue$GTEx_tissue)
samples_x_tissue$GTEx_tissue = paste('inferred', samples_x_tissue$GTEx_tissue)



# PLOTS MAIN ##################################################################################################
x_scale = c(0.47, 1)
y_scale = c(0, num_TFs)




idx = which(df$regulon_group == 'GTEx tissue-specific')
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale, legend_position = 'bottom')
P$aucs$regulon_dataset = gsub('_COMBAT', '', P$aucs$regulon_dataset)
P$aucs$sample_x_tissue = samples_x_tissue$n_samples[ match(P$aucs$regulon_dataset, samples_x_tissue$GTEx_tissue) ]
P$aucs
cortest = cor.test(P$aucs$AUC, P$aucs$sample_x_tissue)
label = paste(paste('cor =', round(cortest$estimate, 2)), 
              paste('p-value =', signif(cortest$p.value, 2)), sep = '\n')
P = ggplot(P$aucs, aes(x=sample_x_tissue, y=AUC, size = coverage)) + 
  geom_smooth(method = "lm", show.legend = F, color = 'red4', fill = 'gray90') +
  geom_point() + 
  annotate('text', label = label, x= Inf, y=Inf, vjust = 1.2, hjust = 1.2) +
  ylab('Area under the PR curve (AUC)') + xlab('#samples') +
  scale_size(name = 'perturbed TFs coverage', limits = c(1, 40), breaks = c(1, 5, 10, 20, 40)) +
  # scale_color_continuous(name = 'perturbed TFs coverage') +
  ggtitle(plot_title) +
  mytheme 
P
save(P, file = paste('paper/figures/GR_revision_AUC_vs_sample', plot_title,'.rdata', sep = ''))
