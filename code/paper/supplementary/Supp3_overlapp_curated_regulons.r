rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')



# Load network
network_files = list.files('data/TF_target_sources/curated_databases/', recursive = T, pattern = 'network', full.names = T) %>% 
  grep('sif$', ., value=T) %>% 
  grep('signed', ., invert =T,  value=T)  %>% 
  grep('cnio', ., invert =T,  value=T)  
networks = lapply(network_files, read.delim, stringsAsFactors = F, sep = '\t')


# retrieve TF
TFs = load_TFs_census()
x = lapply(networks, function(xx) intersect(xx[,1], TFs) )
names(x) = gsub('data/TF_target_sources/curated_databases/', '', network_files) 
names(x) = gsub('/network.sif', '', names(x))
sapply(x, length)
length(unique(unlist(x)))
df = melt(x)
names(df) = c('TF', 'dataset')
df$value = 1
m = dcast(df, TF~dataset, fill = 0)
names(m)[1] = 'Identifier'
png('paper/figures/supplementary/FigureS3.png', res = 300, width = 3000, height = 2200)
upset(m, sets = colnames(m)[-1], main.bar.color = 'gray20', sets.bar.color = my_color_palette$EMBL[4], nintersects = 59,
      empty.intersections = F, set_size.angles = 90, number.angles = 25, 
      order.by = "freq", point.size = 2.5, line.size = 0.5, mb.ratio = c(.6, .4), text.scale = c(1.5, 1.2, 1.5, 1.2, 1.5, 1),
      mainbar.y.label = 'shared TFs', sets.x.label = 'total TFs') # 9000x3100
dev.off()
