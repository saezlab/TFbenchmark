rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')




path = 'data/TF_target_sources/curated_databases/'
filesf = list.files(path, pattern = 'network.sif', recursive = T, full.names = T)
files = list.files(path, pattern = 'network.sif', recursive = T, full.names = F)
l = lapply(filesf, read.delim, header = T, stringsAsFactors = F)
names(l) = sapply(strsplit(files, split = "\\/"), head, 1)



df = melt(l, id.vars = names(l[[1]]))
df$id = paste(df$TF, df$target)


# TFs coverage
tf_df = unique(df[,c('TF', 'L1')])
tf_df$value = 1
m = dcast(tf_df, TF~L1, fill = 0)
names(m)[1] = 'Identifier'
upset(m, sets = colnames(m)[-1], main.bar.color = 'gray20', sets.bar.color = "#56B4E9", order.by = "freq", show.numbers = F, point.size = 5, line.size = 1, mb.ratio = c(.7, .3), text.scale = 2, nintersects = 50) # 9000x3100
ggsave(filename = paste(path, 'TFs_coverage.png', sep=''), dpi = 300, width = 20, height = 20)



# interactions coverage
tf_df = unique(df[,c('id', 'L1')])
tf_df$value = 1
m = dcast(tf_df, id~L1, fill = 0)
names(m)[1] = 'Identifier'
upset(m, sets = colnames(m)[-c(1,8)], main.bar.color = 'gray20', sets.bar.color = "#56B4E9", order.by = "freq", show.numbers = F, point.size = 5, line.size = 1, mb.ratio = c(.5, .5), text.scale = 2, nintersects = 50) # 9000x3100
ggsave(filename = paste(path, 'interactions_coverage.png', sep=''), dpi = 300, width = 20, height = 20)
