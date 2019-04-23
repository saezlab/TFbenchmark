rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/plots.r')



mdf = rbind(data.frame(get(load(file = 'paper/figures/Figure4/enrichment_b1.rdata')), benchmark = 'b1'), 
            data.frame(get(load(file = 'paper/figures/Figure4/enrichment_b2.rdata')), benchmark = 'b2'),
            data.frame(get(load(file = 'paper/figures/Figure4/enrichment_b3.rdata')), benchmark = 'b3'))


x = subset(mdf, padj < 0.05 & abs(ES) > 0.5 &
             L2 %in% c('curated databases', 'inferred GTEx\nconsensus', 'TFBS', 'ChIP-seq'))
x$regulon_evidence = x$L2
x = x[ order(x$ES), ]

x$label = paste(x$L1, x$feature, sep = ' --- ') %>% gsub('_', ' ', .)
gtitle = 'TF annotation - Enrichment'
P = ggplot(x[ grep('class', x$label, invert = T) , ], aes(y = ES, x = label, fill = pval ) ) + geom_bar(stat = 'Identity')  + 
  coord_flip() + geom_hline(yintercept = 0) +
  facet_grid(regulon_evidence~benchmark, scales = 'free_y', space = 'free', as.table = F) +
  scale_fill_gradient(low = brewer.pal(6, 'Blues')[6], high =  brewer.pal(6, 'Blues')[2], name = 'pvalue') +
  # scale_y_continuous(limits = c(-2, 2))  +
  labs(title = gtitle) + xlab('')  +
  mytheme + theme(strip.text.y = element_text(angle=0, hjust = 0), 
                  legend.key.width=unit(1.5,"cm"),
                  plot.subtitle = element_text(face = "italic", hjust = 0.5))
P
ggsave(filename = 'paper/figures/Figure4/Figure4_b.png', dpi = 300, height = 10, width = 10)


P = ggplot(x[ grep('class', x$label, invert = F) , ], aes(y = ES, x = label, fill = pval ) ) + geom_bar(stat = 'Identity')  + 
  coord_flip() + geom_hline(yintercept = 0) +
  facet_grid(regulon_evidence~benchmark, scales = 'free_y', space = 'free', as.table = F) +
  scale_fill_gradient(low = brewer.pal(6, 'Blues')[6], high =  brewer.pal(6, 'Blues')[2], name = 'pvalue') +
  # scale_y_continuous(limits = c(-2, 2))  +
  labs(title = gtitle) + xlab('')  +
  mytheme + theme(strip.text.y = element_text(angle=0, hjust = 0), 
                  legend.key.width=unit(1.5,"cm"),
                  plot.subtitle = element_text(face = "italic", hjust = 0.5))
P
ggsave(filename = 'paper/figures/Figure4/Figure4_c.png', dpi = 300, height = 15, width = 10)






# P = plot_grid(p3a+theme(legend.position = 'none', axis.title.x = element_blank()), 
#               p3b+theme(legend.position = 'none', title = element_blank()), 
#               p3c+theme(plot.title = element_blank(), plot.subtitle= element_blank()), 
#               ncol = 1, align = 'v', labels = LETTERS[2:4], 
#               rel_heights = c(1.75, 1.5, 2.25), label_size = 24)
# P
# save_plot(P, filename = 'paper/figures/Figure3/Figure3_2ndHalf.png', dpi = 300, base_height = 10, base_width = 9)
