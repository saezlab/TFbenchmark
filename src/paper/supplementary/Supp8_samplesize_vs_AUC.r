rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/regulons_QC/4_plot_benchmark_results/lib_plots_regulonsQC.r')


plot_files = list.files(path = home, pattern = 'GR_revision_AUC_vs_sample', recursive = T) %>% grep('rdata', ., value = T)
my_plots = lapply(plot_files, function(x) get(load(x))+ theme(legend.position = 'none') )
P_nolegend = plot_grid(plotlist = my_plots, nrow = 3)
P = plot_grid( P_nolegend, get_legend(get(load(plot_files[[1]]))), rel_heights = c(3, .3), nrow = 2)
save_plot(plot = P, filename = 'paper/figures/GR_revision_samplesize_vs_AUC.png', base_height = 17, base_width = 5.5)


