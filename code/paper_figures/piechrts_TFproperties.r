rm(list = ls())
source('code/utils.r')




# Properties
plot_TFrole_pie = function(TFs, var = 'regulatory_effect', colors = my_color_palette$EMBL){
  df = melt(TFrole_genesets[[var]])
  df = subset(df, value %in% TFs)
  if( any( ! TFs %in% df$value ) )
    df = rbind(df, data.frame(value = setdiff(TFs, df$value), L1='unknown', stringsAsFactors = F))
  mdf = ddply(df, 'L1', nrow)
  mdf$percent = mdf$V1 / sum(mdf$V1)
  if(nrow(mdf)>5){
    mdf = mdf[ order(mdf$percent, decreasing = T) , ]
    mdf$L1[ 6:nrow(mdf) ] = 'other'
    mdf$V1[ 6:nrow(mdf) ] = sum(mdf$V1[ 6:nrow(mdf) ])
    mdf$percent[ 6:nrow(mdf) ] = sum(mdf$percent[ 6:nrow(mdf) ])
    mdf = unique(mdf)
  }
  mdf = mdf[ order(mdf$L1) , ]
  mdf$labels = paste(round(mdf$percent, 2), '%')
  mdf$midpoint = sum(mdf$V1) -  (( c(0, cumsum(mdf$V1[-nrow(mdf)])) ) + mdf$V1/2)
  colors[ mdf$L1 %in% c('other', 'unknown') ] = 'grey'
  mdf$L1 = gsub('_', ' ', mdf$L1)
  pie = ggplot(mdf, aes(x="", y=V1, fill=L1)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)  + 
    scale_fill_manual(values = colors, name = var) +
    # geom_text(aes(x = 1.2, y = midpoint , label = labels), color="black") + 
    theme_minimal(base_size = 18) +
    theme(axis.text=element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) 
  pie
}

load(file = 'data/TF_info/TFrole_genesets.rdata')
names(TFrole_genesets)
TFs = load_TFs_census()
plot(1:14, 1:14, col = my_color_palette$EMBL, cex = 2, pch = 16)
p1 = plot_TFrole_pie(TFs, 'mode_of_regulation', colors = my_color_palette$EMBL[c(10,2,6)])
p2 = plot_TFrole_pie(TFs, 'DNA_binding_mode', colors = my_color_palette$EMBL[c(14,4)])
p3 = plot_TFrole_pie(TFs, 'DNA_binding_specificity', colors = my_color_palette$EMBL[c(7, 3, 12)])
p4 = plot_TFrole_pie(TFs, 'chromatin_regulation_mode', colors = my_color_palette$EMBL[c(8,5,9) ])
# p5 = plot_TFrole_pie(TFs, 'TF_class', colors = c(my_color_palette$EMBL[10:13], my_color_palette$clear))
P = cowplot::plot_grid(p1+theme(legend.position = 'none'), 
                       p2+theme(legend.position = 'none'), 
                       p3+theme(legend.position = 'none'), 
                       p4+theme(legend.position = 'none'), 
                       # p5+theme(legend.position = 'none'), 
                       rel_widths = c(1, 1, 1, 1), nrow = 2)
P
ggsave(P, filename = 'paper/figures/Figure3/Figure3_pies.png', width = 10, height = 10, dpi = 300)
P = cowplot::plot_grid(get_legend(p1), 
                       get_legend(p2), 
                       get_legend(p3), 
                       get_legend(p4),
                       # get_legend(p5),
                       rel_widths = c(1, 1, 1, 1), nrow = 2)
ggsave(P, filename = 'paper/figures/Figure3/Figure3_pies_legend.png', width = 10, height = 10, dpi = 300)


