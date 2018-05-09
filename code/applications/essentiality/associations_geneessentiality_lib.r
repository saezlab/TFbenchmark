
library(reshape2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
source('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/DoRothEA/CODE/assos/lib.r')
source('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/DoRothEA/CODE/my_ggplot_theme.r')
source('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/project_data_loads.r')




# -------------------------------------------------------------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------------------------------------------------------------
covariates.annot = load_covariates()
new_annot = read.csv('/Volumes/GoogleDrive/My Drive/datasets/jsr-gdsc/CanApps_Cell_lines_core_27Nov2017.csv', stringsAsFactors = F)
covariates.annot$gdsc_desc_2 = ''
covariates.annot$gdsc_desc_2[ covariates.annot$CosmicID %in% new_annot$COSMIC.ID ] = new_annot$Cancer.Type[ match(covariates.annot$CosmicID[ covariates.annot$CosmicID %in% new_annot$COSMIC.ID  ], new_annot$COSMIC.ID) ]
covariates.annot = subset(covariates.annot, gdsc_desc_2 != '')
# -------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------




# -------------------------------------------------------------------------------------------------------------------------------
# Assos
# -------------------------------------------------------------------------------------------------------------------------------
analyse_TF_essentiality_associations = function(activities, essentiality, tissue = 'PANCANCER'){
  samples = colnames(activities)
  message('Defining: ', tissue, ' samples')
  if( tissue %in% covariates.annot$gdsc_desc_2 )
    samples = covariates.annot$CosmicID[ covariates.annot$gdsc_desc_2 == tissue]
  if( tissue %in% covariates.annot$gdsc_desc_1)
    samples = covariates.annot$CosmicID[ covariates.annot$gdsc_desc_1 == tissue]
  if( tissue %in% covariates.annot$Study.Abbreviation)
    samples = covariates.annot$CosmicID[ covariates.annot$Study.Abbreviation == tissue]
  samples = intersect(samples, colnames(activities))
  samples = intersect(samples, colnames(essentiality))
  
  message(length(samples), ' samples considered')
  message(nrow(essentiality), ' essential genes considered')
  message(nrow(activities), ' TFs considered')
  
  if ( length(samples) > 3 ){ 
    ASSOS_RESULTS = list()
    for (TF in rownames(activities) ){
      message(' - ', TF)
        X = activities
        X = scale(X[ TF, samples])[,1] # Centering and scaling ES for linear regression
        names(X) = samples
        ASSOS = apply(essentiality, 1, function(Y){
          # i = which(apply(essentiality, 1, function(ic50) all(Y[ ! is.na(Y) ]==ic50[ ! is.na(ic50) ]) )) # Debug
          # message(i)
          names(Y) = colnames(essentiality)
          Y = Y[samples]
          Y = Y[ ! is.na(Y) ]
          if( length(Y) > 4 )
            IC50_gdscANOVA(X = X, Y = Y[samples], covariates.annot = covariates.annot, lm = T)
        } )
        ASSOS = Filter(Negate(is.null), ASSOS)
        if( ! is.null(ASSOS) ){
          mASSOS = melt(ASSOS, id.vars = names(ASSOS[[1]]))
          names(mASSOS)[ names(mASSOS) == 'L1' ] = 'Essential_gene'
          mASSOS$TF_fdr = NA
          mASSOS$TF_fdr <- p.adjust(mASSOS$X_lm_pval, method = 'fdr')
          mASSOS = mASSOS[, order(colnames(mASSOS)) ]
          ASSOS_RESULTS[[TF]] = mASSOS
        }
      }
    ASSOS_RESULTS = Filter(Negate(is.null), ASSOS_RESULTS)
    if( ! is.null(ASSOS_RESULTS) & length(ASSOS_RESULTS) > 1 ){
      mdf = melt(ASSOS_RESULTS, id.vars = names(ASSOS_RESULTS[[1]]))
      names(mdf)[ names(mdf) == 'L1' ] = 'TF'
      mdf = mdf[, order(colnames(mdf)) ]  
      mdf$global_FDR <- p.adjust(mdf$X_lm_pval, method = 'fdr')
      mdf = mdf[ order( mdf$global_FDR), ]
      return(mdf)
    }
  }
}
# -------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------




# -------------------------------------------------------------------------------------------------------------------------------
# PLOTS
# -------------------------------------------------------------------------------------------------------------------------------
plot_volcano = function(associations_df, main = '', th = .05, repel_n = 1:10){
  if( ! is.null(associations_df) ) {
    fdrcol = grep('fdr', names(associations_df), ignore.case = F) 
    associations_df$effect = associations_df$X_lm_coefficient
    associations_df$log10pval = - log10(associations_df[, fdrcol])
    associations_df$pval = associations_df$X_lm_pval
    
    # associations_df = associations_df[ order(associations_df$pval), ]
    associations_df$color = 'gray'
    associations_df$color[ associations_df[,fdrcol] < th & associations_df$effect < 0 ] = 'green4'
    associations_df$color[ associations_df[,fdrcol] < th & associations_df$effect > 0 ] = 'red'
    associations_df$labels = ''
    associations_df$labels[repel_n] = paste(associations_df$TF[repel_n], " - ", 
                                                 associations_df$Essential_gene[repel_n],
                                                 sep='')
    p = ggplot(associations_df, aes(x=effect, y=log10pval, color = color, size = (log10pval^2)*.25)) + geom_point(alpha = .65) + 
      scale_colour_identity() + 
      xlab('effect size') + ylab(expression(-log[10]~adj.pvalue)) +
      mytheme + theme(legend.position = 'none', text = element_text(size = 14)) + ggtitle(main)
    p = p + geom_text_repel(data=associations_df[ as.numeric(repel_n), ], aes(label = labels), color = 'black',
                            size = 2.7, 
                            fill = rgb(red = 0, blue = 0, green = 0, alpha = 0, maxColorValue = 250), 
                            segment.size = 0.15,
                            parse = T, 
                            label.size = 0, 
                            min.segment.length = unit(0, "lines"),
                            point.padding = unit(0.5, "lines"),
                            box.padding = unit(0.5, "lines")) 
  return(p)
  }
}

plot_volcano_tissues = function(associations_df, main = '', path = NULL, lm = T, th = .05, repel_n = 10){
  tissueTCGA_hash = load_tissues_colors(TCGA = T)
  if( ! is.null(associations_df) ) {
    associations_df$effect = associations_df$X_lm_coefficient
    fdrcol = grep('fdr', names(associations_df), ignore.case = T) 
    associations_df$log10pval = - log10(associations_df[, fdrcol])
    associations_df$pval = associations_df[, fdrcol]
    associations_df = associations_df[ order(associations_df$pval), ]
    associations_df = associations_df[ ! is.na(associations_df$pval), ]
    associations_df$color = 'gray90'
    idx = which( associations_df[,fdrcol] < th & associations_df$L1 %in% tissueTCGA_hash$tissue)
    associations_df$color[ idx ] = tissueTCGA_hash$color[ match(associations_df$L1[ idx  ], tissueTCGA_hash$tissue) ]
    associations_df$labels = ''
    idx = 1:repel_n[1]
    if( length(repel_n) > 1 )
      idx = repel_n
    associations_df$labels[idx] = paste("", associations_df$L1[idx], "[", associations_df$Drug_name[idx], "]^",  gsub(' ', '', associations_df$TF[idx]),  sep='')
    p = ggplot(associations_df[ associations_df$X_lm_pval < 0.75 & abs(associations_df$effect) < 7.5 , ], aes(x=effect, y=log10pval, color = color, size = (log10pval^2)*.01)) + 
      geom_point(alpha = .65) + 
      geom_point(data = associations_df[ associations_df[, fdrcol] < th & abs(associations_df$effect) < 7.5 & ! associations_df$L1 %in% c('SCLC', 'OV', 'COREAD', 'BRCA'), ], alpha = .65) + 
      geom_hline(yintercept = -log10(th), colour="gray20", linetype = "longdash") +
      annotate('text', x =-Inf,  y = -log10(th), label = paste('FDR <', th), colour="gray20", hjust = -.1, vjust = -.5, size = 4.5) +
      scale_colour_identity()+
      guides(size=FALSE) +
      scale_y_continuous(breaks=integer_breaks) +
      xlab('effect size') + ylab(expression(-log[10]~adj.pvalue)) +
      mytheme + 
      ggtitle(main)
    p = p + geom_label_repel(data=associations_df[ idx, ], aes(label = labels), color = 'black', force = 0.5, 
                             size = 3, 
                             fill = rgb(red = 0, blue = 0, green = 0, alpha = 0, maxColorValue = 250), 
                             segment.size = 0.18,
                             parse = T, 
                             label.size = 0, 
                             min.segment.length = unit(0, "lines"),
                             point.padding = unit(0.5, "lines"),
                             box.padding = unit(0.5, "lines")) 
    p
    if( ! is.null(path) )
      ggsave(filename = paste(path, '/', gsub(' ', '_', main), '_volcano.png', sep = ''), dpi = 300, width = 6, height = 10)
  }
  return(p)
}

pretty_drug_labels2volcano = function(values){
  values = sapply(values, function(x) paste(setdiff(unlist(strsplit(x,  split = ',')), c('', '.', 'receptor')), collapse = ','))
  values = gsub('MAP2K1 \\(MEK1), MAP2K2 \\(MEK2)', 'MEK1 and MEK2', values)
}

plot_heatmap = function(associations_df, main = '', path = '', lm = T, th = .01){
  if( ! is.null(associations_df) ) {
    associations_df = associations_df[associations_df[ fdrcol ] < th, ]
    if( nrow(associations_df) > 1){
      if (lm)
        associations_df$value = associations_df$X_lm_coefficient
      if (! lm)
        associations_df$value = associations_df$FEATURE_IC50_effect_size
      associations_df$Drug_name = paste(associations_df$Drug_name, ' [', associations_df$Drug_targets, ']', sep = '')
      X = acast(associations_df, Drug_name~TF, fill = 0, value.var='value')
      paletteLength = 101
      # myColor = colorRampPalette(c(brewer.pal(5, 'Greens')[5], 'white', brewer.pal(5, 'Reds')[5]))(paletteLength)
      myColor =  colorRampPalette(c(my_color_palette[1], my_color_palette[2]))(101)
      myBreaks = c(seq(min(X), 0, length.out=ceiling(paletteLength/2)), 
                   seq(max(X)/paletteLength, max(X), length.out=floor(paletteLength/2)))
      myColor[ myBreaks == 0] = 'white'
      pheatmap(X, color = myColor, breaks = myBreaks, main = main, filename = paste(path, '/', gsub(' ', '_', main), '_heatmap.png', sep = ''), cellwidth = 10, cellheight = 7, border_color = 'gray90', fontsize_row = 8)
      #pheatmap(X, color = c(rev(brewer.pal(5, 'Reds')), 'white', brewer.pal(5, 'Greens')), main = main, filename = paste(path, '/', main, '_heatmap.png', sep = ''), cellwidth = 10, cellheight = 7, border_color = 'gray90', fontsize_row = 8)
    }
  }
}



get_heatmap_TF_clusters = function(associations_df, main = '',  path = '', k, lm = T){
  if( ! is.null(associations_df) ) {
    associations_df = associations_df[associations_df[ fdrcol ] < 0.01, ]
    if( nrow(associations_df) > 1){
      if (lm)
        associations_df$value = associations_df$X_lm_coefficient
      if (! lm)
        associations_df$value = associations_df$FEATURE_IC50_effect_size
      associations_df$Drug_name = paste(associations_df$Drug_name, ' [', associations_df$Drug_targets, ']', sep = '')
      X = acast(associations_df, Drug_name~TF, fill = 0, value.var='value')
      paletteLength = 101
      myColor = colorRampPalette(c(brewer.pal(5, 'Greens')[5], 'white', brewer.pal(5, 'Reds')[5]))(paletteLength)
      myBreaks = c(seq(min(X), 0, length.out=ceiling(paletteLength/2)), 
                   seq(max(X)/paletteLength, max(X), length.out=floor(paletteLength/2)))
      myColor[ myBreaks == 0] = 'white'
      res = pheatmap(X, color = myColor, breaks = myBreaks, main = main, cellwidth = 10, cellheight = 7, border_color = 'gray90', fontsize_row = 8)
      clust = cutree(res$tree_col,  k = k)
      clust = lapply(unique(clust), function(x) names(clust[ clust == x] ))
      return(clust)
    }
  }
}


plot_ic50_tf_corr = function(TF, drug_id, tissue, rcoef, pcoef){
  slea_results$ES = t(scale(t(slea_results$ES)))
  samples = covariates.annot$CosmicID
  if( tissue %in% covariates.annot$Study.Abbreviation )
    samples = covariates.annot$CosmicID[ covariates.annot$Study.Abbreviation == tissue ]
  if( tissue %in% covariates.annot$gdsc_desc_1 )
    samples = covariates.annot$CosmicID[ covariates.annot$gdsc_desc_1 == tissue ]
  samples = Reduce(intersect, list(colnames(IC50s[ , ! is.na(IC50s[drug_id,]) ]), samples, colnames(slea_results$ES)))
  ic50 = IC50s[ drug_id , samples]
  tf = slea_results$ES[TF, samples]
  df = data.frame(TF_activity=tf, IC50=ic50, stringsAsFactors = F)
  col = my_color_palette[4]
  if( rcoef > 0)
    col =  'red'
  drug_name = drug.properties$DRUG_NAME[ drug.properties$DRUG_ID == drug_id ]
  drug_name = gsub('screen', '', drug_name)
  p = ggplot(df, aes(x=TF_activity, y=IC50)) + 
    geom_point(color = col, size=1.5, alpha=.3) +
    stat_smooth(method = 'lm', color = col) + 
    ylab(bquote(.(drug_name)[' log IC50']))  + 
    xlab(paste(TF, 'activity')) +  mytheme +
    annotate('text', label = paste('coef = ', round(rcoef,3), '; p = ', signif(pcoef, 3), sep =''), x = Inf, y=Inf, vjust=1.2, hjust=1, size=3.5) +
    annotate('text', label = tissue, x =-Inf, y=-Inf, vjust=-1, hjust=-1, size=4.5, color = 'red3')
  return(p)
}
# -------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------



