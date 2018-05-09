rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')



# Load network
network = read.csv(file = 'data/TF_target_sources/omnipath_scores/database.csv', stringsAsFactors = F)



# retrieve TF
x = lapply(unique(network$score), function(sc){
  unique(subset(network, score %in% LETTERS[1:which(LETTERS == sc)] )$TF)
})
names(x) = unique(network$score)
sapply(x, length)
df = melt(x)
names(df) = c('TF', 'score')
mdf = data.frame(table(df$score))
p1 = ggplot(mdf, aes(x = Var1, fill = Var1, y = Freq, label = Freq) ) + geom_bar(stat = 'Identity') + 
  scale_fill_manual(values = brewer.pal(5, 'Paired')[c(2,1,4,3,5)]) + 
  geom_text(y = 100) +
  mytheme + theme(legend.position = 'none', axis.title.x = element_blank()) + coord_flip() + 
  ggtitle('TFs coverage') + xlab('score threshold')
p1





# retrieve TFTG
network$TFtarget = paste(network$TF, network$target)
x = lapply(unique(network$score), function(sc){
  unique(subset(network, score %in% LETTERS[1:which(LETTERS == sc)] )$TFtarget)
})
names(x) = unique(network$score)
sapply(x, length)
df = melt(x)
names(df) = c('TFtarget', 'score')
mdf = data.frame(table(df$score))
p2 = ggplot(mdf, aes(x = Var1, fill = Var1, y = Freq, label = Freq) ) + geom_bar(stat = 'Identity') + 
  scale_fill_manual(values = brewer.pal(5, 'Paired')[c(2,1,4,3,5)]) + 
  geom_text(y=100000) +
  mytheme + theme(legend.position = 'none', axis.title.x = element_blank()) + coord_flip() + 
  ggtitle('TF-targets coverage') + xlab('score threshold')
p2



P = plot_grid(p1, p2, ncol = 1) 
P
ggsave(P, filename = 'paper/figures/Figure4b.png', dpi = 300, width = 6, height = 6)
