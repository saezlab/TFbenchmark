library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(pheatmap)
library(ggbeeswarm)


mytheme =   theme_light(15) + 
  theme(legend.position = 'bottom', legend.key = element_blank(), 
        legend.background = element_rect(fill=alpha('white', 0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill=alpha('white', 0)),
        text = element_text(family = 'sans'),
        strip.background =  element_blank(), 
        strip.text.x = element_text(color ='black', family = 'sans'), 
        strip.text.y = element_text(color ='black', family = 'sans'),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.text=element_text(size=12, color ='black'),
        axis.title=element_text(size=12, color ='black'),
        plot.title = element_text(size=15, hjust = 0.5))


my_color_palette = list(EMBL = c('#00777D', '#A2C012', '#b78c00',  '#0081AD' , '#E94949',  brewer.pal(name = 'Dark2', n=8), 'cornflowerblue'),
                        shiny = c('#2B8C7F', '#BFC42C','#C4712D', '#C42D44', '#363636', '#789895',  '#0081AD'),
                        clear = c('#789895' , '#676d7e', '#E94949', '#A2C012', '#DCD6B4', '#E0AD00',  '#0081AD'))
  
  

my_venn_plot = function(l, filename, main = ''){
  require(VennDiagram)
  venn.diagram(l, 
               fill = c("steelblue2", "mediumseagreen"), alpha = 0.5, cex = 2.75, lty =0, label.col = c("steelblue4", "white", "springgreen4"), fontfamily = "sans", 
               main = main, main.fontfamily = "sans", main.cex = 2.6,
               cat.pos = c(-30,30), cat.dist = 0.05, cat.col = c("steelblue3", "mediumseagreen"), cat.cex = 2.5, cat.fontfamily = "sans",
               filename = filename)
  filepath = unlist(strsplit(filename, '/'))
  filepath = paste(filepath[-length(filepath)], collapse = '/')
  unlink(list.files(path =  filepath, pattern = '[0-9]\\.log$', full.names = T))
}


integer_breaks = function(x){
  newx = seq(floor(min(x)), ceiling(max(x)))
  if( length(newx) > 10 )
    newx = seq(floor(min(x)), ceiling(max(x)), by = 2)
  if( length(newx) > 10 )
    newx = seq(floor(min(x)), ceiling(max(x)), by = 5)
  newx = unique(round(newx))
}




load_tissues_colors = function(TCGA = F){
   # COSMIC
    source('CODE/TF_scores_methods_v2/lib_manage_data.r')
    covariates = load_sample_variables_gdsc()
    covariates$gdsc_desc_1[ grep('lung', covariates$gdsc_desc_1) ]= 'lung'
    tissue_hash = data.frame(tissue = sort(unique(covariates$gdsc_desc_1)), stringsAsFactors = F)
    tissue_hash$val = 1
    tissue_hash$color = c("#FB9A99", "#E31A1C", "#7570B3", "#E7298A", "#B15928", 
                          "#A6761D", "#E6AB02", "cornflowerblue", "#66A61E", "#A6CEE3", 
                          "#1F78B4", "#B2DF8A", "#1B9E77", "#666666", "#FFFF99", 
                          "#6A3D9A", "#CAB2D6", "#A6D854", "#E78AC3") 
    if(TCGA){
      # TCGA
      tissueTCGA_hash = unique(covariates[ order(covariates$Study.Abbreviation) , c('gdsc_desc_1', 'Study.Abbreviation') ])
      tissueTCGA_hash$gdsc_desc_1[ grep('lung', tissueTCGA_hash$gdsc_desc_1) ]= 'lung'
      tissueTCGA_hash$color = tissue_hash$color[ match(tissueTCGA_hash$gdsc_desc_1, tissue_hash$tissue)]
      names(tissueTCGA_hash)[2] = 'tissue'
      tissueTCGA_hash$color[ tissueTCGA_hash$tissue %in% c('', 'UNABLE TO CLASSIFY') ] = 'gray'
      tissueTCGA_hash = unique(tissueTCGA_hash[, 2:3])
      tissueTCGA_hash = tissueTCGA_hash[ ! duplicated(tissueTCGA_hash$tissue) , ]
      tissueTCGA_hash = rbind(tissueTCGA_hash, tissue_hash[ , c(1,3)])
      ggplot(tissueTCGA_hash, aes(x = tissue, y = 10, color = color)) + geom_point(size=10) + scale_color_identity() + theme_bw(10)
      return(tissueTCGA_hash)
    }
    ggplot(tissue_hash, aes(x = tissue, y = 10, color = color)) + geom_point(size=10) + scale_color_identity() + theme_bw(10)
    return(tissue_hash)
}


library(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


