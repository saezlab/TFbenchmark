rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/regulons_QC/4_plot_benchmark_results/lib_plots_regulonsQC.r')




# LOAD ##################################################################################################
load(file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors_formated.rdata')
df$rank_nes = 1 - df$rank_nes
df$NES = df$NES * -(1)
df = df[ df$perturbed_TF %in% df$TF, ]
df = subset(df, ! regulon_dataset %in% c('all_weighted', 'cnio', 'inferred_dynamic', 'evidence4', 'TOP') )






# Define test and plot features
balanced_accuracy = F
plot_title = 'B1'
outdir = 'results/regulons_QC/GEOperturbations/'





# Counts
df = subset(df, is_TF_DE_expected | ! perturbation_treatment %in% c("shRNA", "siRNA", "overexpression" ) )



# Load samples x GTEx tissue
samples_x_tissue = read.csv('data/TF_target_sources/inferred/GTEx/samples_x_tissue.csv')
samples_x_tissue$GTEx_tissue = gsub(' ', '_', samples_x_tissue$GTEx_tissue)
samples_x_tissue$GTEx_tissue = paste('inferred', samples_x_tissue$GTEx_tissue)




# Confidence vs coverage data ##################################################################################################
P = auc2coverage_plot(df, balanced_accuracy = balanced_accuracy)
x_scale = c(min(P$aucs$AUC), signif(max(P$aucs$AUC),1))
if( min(P$aucs$AUC) > 0.475 )
  x_scale = c(0.47, max(P$aucs$AUC)+0.01)
y_scale = c(0, max(P$aucs$coverage))





idx = which(df$regulon_group == 'GTEx tissue-specific')
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale, legend_position = 'bottom')
P$aucs
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
  scale_size(name = '', limits = c(1, 40), breaks = c(1, 5, 10, 20, 40)) +
  # scale_color_continuous(name = 'perturbed TFs coverage') +
  ggtitle(plot_title) +
  mytheme 
P
save(P, file = paste('paper/figures/GR_revision_AUC_vs_sample', plot_title,'.rdata', sep = ''))
