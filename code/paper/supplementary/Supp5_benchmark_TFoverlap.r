rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/plots.r')
library(UpSetR)



x = list()
x$benchmark_1 = read.table('data/regulons_QC/benchmark1_TFs.txt', stringsAsFactors = F, header = F)[,1]
x$benchmark_2 = read.table('data/regulons_QC/benchmark2_TFs.txt', stringsAsFactors = F, header = F)[,1]
x$benchmark_3 = read.table('data/regulons_QC/benchmark3_TFs.txt', stringsAsFactors = F, header = F)[,1]
df = melt(x)
names(df) = c('TF', 'dataset')
df$value = 1
m = dcast(df, TF~dataset, fill = 0)
names(m)[1] = 'Identifier'

png('paper/figures/supplementary/Supplemental_Fig_S5.png', res = 300, width = 2600, height = 1200)
upset(m, sets = colnames(m)[-1], main.bar.color = 'gray20', sets.bar.color = my_color_palette$EMBL[4], 
      empty.intersections = F, set_size.angles = 90, number.angles = 25, 
      order.by = "freq", point.size = 2.5, line.size = 0.5, mb.ratio = c(.6, .4), text.scale = c(1.5, 1.2, 1.5, 1.2, 1.5, 1),
      mainbar.y.label = 'shared TFs', sets.x.label = 'total TFs') # 9000x3100
dev.off()
