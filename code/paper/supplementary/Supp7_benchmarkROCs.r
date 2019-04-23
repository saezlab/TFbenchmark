rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


pA = list.files(path = 'results/regulons_QC/',  recursive = T, full.names = T, pattern = 'chip') %>% 
  grep(., pattern = 'rdata', value = T) %>% 
  grep(., pattern = 'ROC_', value = T)
pA = lapply(pA, function(x) get(load(x)) + theme(axis.title.x = element_blank())  )
pA[[2]] = pA[[2]] + ggtitle('b1')
pA[[3]] = pA[[3]] + ggtitle('b2')
pA[[1]] = pA[[1]] + ggtitle('b3')

# pA = plot_grid(pA[[2]], pA[[3]], pA[[1]], align = 'h', nrow = 1)



pB = list.files(path = 'results/regulons_QC/',  recursive = T, full.names = T, pattern = 'TFBS_size') %>% 
  grep(., pattern = 'rdata', value = T) %>% 
  grep(., pattern = 'ROC_', value = T)
pB = lapply(pB, function(x) get(load(x)) + theme(title = element_blank()))
# pB = plot_grid(pB[[2]], pB[[3]], pB[[1]], align = 'h', nrow = 1)


pC = list.files(path = 'results/regulons_QC/',  recursive = T, full.names = T, pattern = 'TFBS_2') %>% 
  grep(., pattern = 'rdata', value = T) %>% 
  grep(., pattern = 'ROC_', value = T)
pC = lapply(pC, function(x) get(load(x)) + theme(title = element_blank()) )
# pC = plot_grid(pC[[2]], pC[[3]], pC[[1]], align = 'h', nrow = 1)


pD = list.files(path = 'results/regulons_QC/',  recursive = T, full.names = T, pattern = 'inferred') %>% 
  grep(., pattern = 'rdata', value = T) %>% 
  grep(., pattern = 'ROC_', value = T)
pD = lapply(pD, function(x) get(load(x)) + theme(title = element_blank()) )
# pD = plot_grid(pD[[2]], pD[[3]], pD[[1]], align = 'h', nrow = 1)



pE = list.files(path = 'results/regulons_QC/',  recursive = T, full.names = T, pattern = 'curated') %>% 
  grep(., pattern = 'rdata', value = T) %>% 
  grep(., pattern = 'ROC_', value = T)
pE = lapply(pE, function(x) get(load(x)) + theme(plot.title = element_blank()))
# pE = plot_grid(pE[[2]], pE[[3]], pE[[1]], align = 'h', nrow = 1)




# P = plot_grid(pA, pB, pC, pD, pE, ncol = 1, labels = 'AUTO', label_size = 25, rel_heights = c(1,1.5,1.5,1,2.5), align = 'v')
P = plot_grid(pA[[2]], pA[[3]], pA[[1]], 
              pB[[2]], pB[[3]], pB[[1]], 
              pC[[2]], pC[[3]], pC[[1]], 
              pD[[2]], pD[[3]], pD[[1]], 
              pE[[2]], pE[[3]], pE[[1]], 
              ncol = 3, 
              rel_heights = c(1.1,1.5,1.5,1,2.4), 
              align = 'hv', 
              labels = c('A', '', '', 'B', '', '', 'C', '', '', 'D', '', '', 'E', '', ''), 
              label_size = 25)
save_plot(P, filename = 'paper/figures/supplementary/GarciaAlonso_FigureS3.png', dpi = 300, base_height = 15, base_width = 16)
