rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')



# Load network
network = read.csv(file = 'data/TF_target_sources/omnipath_scores/database_20180915.csv', stringsAsFactors = F)
length(unique(network$TF))
length(unique(network$target))
pancancer_network = read.csv(file = 'data/TF_target_sources/omnipath_scores/database_PANCANCER_20181030.csv', stringsAsFactors = F)
length(unique(pancancer_network$TF))
length(unique(pancancer_network$target))



# retrieve TF
x = list()
x$normal = lapply(unique(network$score), function(sc){
  unique(subset(network, score %in% LETTERS[1:which(LETTERS == sc)] )$TF)
})
names(x$normal) = unique(network$score)
x$pancancer = lapply(unique(pancancer_network$score), function(sc){
  unique(subset(pancancer_network, score %in% LETTERS[1:which(LETTERS == sc)] )$TF)
})
names(x$pancancer) = unique(pancancer_network$score)

sapply(x$normal, length)
sapply(x$pancancer, length)

df = melt(x)
head(df)
names(df) = c('TF', 'score', 'dataset')
mdf = ddply(df, 'dataset', function(xx)  data.frame(table(xx$score)))
p1 = ggplot(mdf, aes(x = Var1, fill = dataset, y = Freq, label = Freq, group = dataset) ) + 
  geom_bar(stat = 'Identity', position = 'dodge', color = 'white', alpha = 0.8) + 
  # scale_fill_manual(values = brewer.pal(5, 'Paired')[c(2,1,4,3,5)]) + 
  scale_fill_manual(values = c(my_color_palette$EMBL[1], brewer.pal(6, 'Paired')[6])) + 
  # geom_text(aes(y = 100), position = position_dodge(width = 1)) +
  mytheme + theme(legend.position = 'none', axis.title.x = element_blank()) + coord_flip() + 
  ggtitle('TFs coverage') + xlab('score threshold')
p1





# retrieve TFTG
network$TFtarget = paste(network$TF, network$target)
length(unique(network$TFtarget))
pancancer_network$TFtarget = paste(pancancer_network$TF, pancancer_network$target)
length(unique(pancancer_network$TFtarget))
x$normal = lapply(unique(network$score), function(sc){
  unique(subset(network, score %in% LETTERS[1:which(LETTERS == sc)] )$TFtarget)
})
names(x$normal) = unique(network$score)
x$pancancer = lapply(unique(pancancer_network$score), function(sc){
  unique(subset(pancancer_network, score %in% LETTERS[1:which(LETTERS == sc)] )$TFtarget)
})
names(x$pancancer) = unique(pancancer_network$score)
sapply(x$normal, length)
sapply(x$pancancer, length)
df = melt(x)
names(df) = c('TFtarget', 'score', 'dataset')
mdf = ddply(df, 'dataset', function(xx)  data.frame(table(xx$score)))
mdf$Freq_lab = mdf$Freq
mdf$Freq_lab[ mdf$Freq > 100000] = ''
p2 = ggplot(mdf, aes(x = Var1, fill = dataset, y = Freq, label = Freq_lab, group = dataset) ) + 
  geom_bar(stat = 'Identity', position = 'dodge', color = 'white', alpha = 0.8) + 
  # scale_fill_manual(values = brewer.pal(5, 'Paired')[c(2,1,4,3,5)]) + 
  scale_fill_manual(values = c(my_color_palette$EMBL[1], brewer.pal(6, 'Paired')[6])) + 
  geom_text(position = position_dodge(width = 1), size = 3.5, hjust = -0.2) +
  mytheme + theme(legend.position = 'bottom', axis.title.x = element_blank()) + coord_flip() + 
  ggtitle('TF-targets coverage') + xlab('score threshold')
p2



P = plot_grid(p1, p2, ncol = 1, rel_heights = c(45,55)) 
P
ggsave(P, filename = 'paper/figures/Figure5b.png', dpi = 300, width = 6, height = 6)
