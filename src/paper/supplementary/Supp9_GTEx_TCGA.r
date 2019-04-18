rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/plots.r')



pB1 = get(load('results/regulons_QC/B1_perturbations/PR_inferred_cancer.rdata'))
pB3 = get(load('results/regulons_QC/B3_cell_lines/PR_inferred_cancer.rdata'))
P = plot_grid(pB1, pB3, labels = 'AUTO')
save_plot(P, filename = 'paper/figures/supplementary/Supplemental_Fig_S9.png', base_height = 5, base_width = 12)
