


pA = list.files(path = 'results/regulons_QC/',  recursive = T, full.names = T, pattern = 'chip') %>% 
  grep(., pattern = 'rdata', value = T)
pA = lapply(pA, function(x) get(load(x)) )
pA = plot_grid(pA[[2]], pA[[3]], pA[[1]], align = 'h', nrow = 1)



pB = list.files(path = 'results/regulons_QC/',  recursive = T, full.names = T, pattern = 'TFBS_size') %>% 
  grep(., pattern = 'rdata', value = T)
pB = lapply(pB, function(x) get(load(x)) )
pB = plot_grid(pB[[2]], pB[[3]], pB[[1]], align = 'h', nrow = 1)


pC = list.files(path = 'results/regulons_QC/',  recursive = T, full.names = T, pattern = 'inferred') %>% 
  grep(., pattern = 'rdata', value = T)
pC = lapply(pC, function(x) get(load(x)) )
pC = plot_grid(pC[[2]], pC[[3]], pC[[1]], align = 'h', nrow = 1)



P = plot_grid(pA, pB, pC, ncol = 1, labels = 'AUTO', label_size = 25)
save_plot(P, filename = 'paper/figures/supplementary/FigureS_single_evidence_ROCs.png', dpi = 300, base_height = 12, base_width = 13)
