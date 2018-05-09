
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
drug.properties = load_drugpops()
IC50s = load_IC50s()
# -------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------




# -------------------------------------------------------------------------------------------------------------------------------
# Assos
# -------------------------------------------------------------------------------------------------------------------------------
analyse_TF_drug_associations = function(tf_activities, tissue = 'pancancer'){
  
  message('Defining: ', tissue, ' samples')
  samples = covariates.annot$CosmicID[ covariates.annot$gdsc_desc_1 != 'not_specified' ]
  if( tissue != 'pancancer')
    samples = covariates.annot$CosmicID[ covariates.annot$gdsc_desc_2 == tissue]
  if( tissue %in% covariates.annot$gdsc_desc_1)
    samples = covariates.annot$CosmicID[ covariates.annot$gdsc_desc_1 == tissue]
  if( tissue %in% covariates.annot$Study.Abbreviation)
    samples = covariates.annot$CosmicID[ covariates.annot$Study.Abbreviation == tissue]
  samples = intersect(samples, colnames(tf_activities))
  samples = intersect(samples, colnames(IC50s))
  
  message(length(samples), ' samples considered')
  message(nrow(IC50s), ' drugs considered')
  message(nrow(tf_activities), ' TFs considered')
  if ( length(samples) <= 3 )
    return(NULL)
   
  # analyze each TF searately
  ASSOS_RESULTS = list()
  for (TF in rownames(tf_activities) ){
    message(' - ', TF)
    X = tf_activities
    X = scale(X[ TF, samples])[,1] # Centering and scaling ES for linear regression
    names(X) = samples
    ASSOS = apply(IC50s, 1, function(Y){
      names(Y) = colnames(IC50s)
      Y = Y[samples]
      Y = Y[ ! is.na(Y) ]
      if( length(Y) > 4 )
        IC50_gdscANOVA(X = X, Y = Y[samples], covariates.annot = covariates.annot, lm = T)
    })
    ASSOS = Filter(Negate(is.null), ASSOS)
    if( ! is.null(ASSOS) ){
      mASSOS = melt(ASSOS, id.vars = names(ASSOS[[1]]))
      names(mASSOS)[ names(mASSOS) == 'L1' ] = 'Drug_id'
      mASSOS$Drug_name = drug.properties$DRUG_NAME[ match(mASSOS$Drug_id, drug.properties$DRUG_ID) ]
      mASSOS$Drug_targets = drug.properties$PUTATIVE_TARGET[ match(mASSOS$Drug_id, drug.properties$DRUG_ID) ]
      mASSOS = mASSOS[, order(colnames(mASSOS)) ]
      ASSOS_RESULTS[[TF]] = mASSOS
    }
  }
  
  # merge data
  ASSOS_RESULTS = Filter(Negate(is.null), ASSOS_RESULTS)
  if( ! is.null(ASSOS_RESULTS) & length(ASSOS_RESULTS) > 1 ){
    TF_drug_associations = melt(ASSOS_RESULTS, id.vars = names(ASSOS_RESULTS[[1]]))
    names(TF_drug_associations)[ names(TF_drug_associations) == 'L1' ] = 'TF'
    TF_drug_associations = TF_drug_associations[, order(colnames(TF_drug_associations)) ]  
    TF_drug_associations$X_lm_fdr <- p.adjust(TF_drug_associations$X_lm_pval, method = 'fdr')
    TF_drug_associations = TF_drug_associations[ order( TF_drug_associations$X_lm_fdr), ]
    return(TF_drug_associations)
  }
}
# -------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------




# -------------------------------------------------------------------------------------------------------------------------------
# PLOTS
# -------------------------------------------------------------------------------------------------------------------------------
plot_volcano = function(TF_drug_associations, main = '', path = NULL, th = .05, repel_n = 1:10){
  if( ! is.null(TF_drug_associations) ) {
    TF_drug_associations$effect = TF_drug_associations$X_lm_coefficient
    TF_drug_associations$log10pval = - log10(TF_drug_associations$X_lm_fdr)
    TF_drug_associations$pval = TF_drug_associations$X_lm_pval
    TF_drug_associations$color = 'gray'
    TF_drug_associations$color[ which(TF_drug_associations$X_lm_fdr < th & TF_drug_associations$effect < 0) ] = 'green4'
    TF_drug_associations$color[ which(TF_drug_associations$X_lm_fdr < th & TF_drug_associations$effect > 0) ] = 'red'
    TF_drug_associations$labels = ''
    TF_drug_associations$labels[repel_n] = paste(TF_drug_associations$TF[repel_n], "~'\n", 
                                                 TF_drug_associations$Drug_name[repel_n], 
                                                 " '['", 
                                                 gsub(' ', '', pretty_drug_labels2volcano(as.character(TF_drug_associations$Drug_targets[repel_n]))), "']", 
                                                 sep='')
    TF_drug_associations$labels = gsub('rescreen', 're', TF_drug_associations$labels)
    p = ggplot(TF_drug_associations, aes(x=effect, y=log10pval, color = color, size = (log10pval^2)*.25)) + geom_point(alpha = .65) + 
      scale_colour_identity() + scale_y_continuous(breaks=integer_breaks) +
      xlab('effect size') + ylab(expression(-log[10]~adj.pvalue)) +
      mytheme + theme(legend.position = 'none', text = element_text(size = 14)) + ggtitle(main)
    if( any(TF_drug_associations$X_lm_fdr < th) )    
      p = p + geom_hline(yintercept = -log10(th) , colour="grey30", linetype = "longdash") + 
      annotate('text', x =Inf,  y = -log10(th), label = 'FDR < 0.05', colour="grey30", hjust = 1, vjust = -.5, size = 3) 
    p = p + geom_text_repel(data=TF_drug_associations[ as.numeric(repel_n), ], aes(label = labels), color = 'black',
                            size = 2.7, 
                            fill = rgb(red = 0, blue = 0, green = 0, alpha = 0, maxColorValue = 250), 
                            segment.size = 0.15,
                            parse = T, 
                            label.size = 0, 
                            min.segment.length = unit(0, "lines"),
                            point.padding = unit(0.5, "lines"),
                            box.padding = unit(0.5, "lines")) 
    if( ! is.null(path) )
      ggsave(filename = paste(path, '/', gsub(' ', '_', main), '_volcano.png', sep = ''), dpi = 300, width = 6, height = 10)
  }
  return(p)
}

plot_volcano_tissues = function(TF_drug_associations, main = '', path = NULL, lm = T, th = .05, repel_n = 10){
  tissueTCGA_hash = load_tissues_colors(TCGA = T)
  if( ! is.null(TF_drug_associations) ) {
    TF_drug_associations$effect = TF_drug_associations$X_lm_coefficient
    fdrcol = grep('fdr', names(TF_drug_associations), ignore.case = T) 
    TF_drug_associations$log10pval = - log10(TF_drug_associations[, fdrcol])
    TF_drug_associations$pval = TF_drug_associations[, fdrcol]
    TF_drug_associations = TF_drug_associations[ order(TF_drug_associations$pval), ]
    TF_drug_associations = TF_drug_associations[ ! is.na(TF_drug_associations$pval), ]
    TF_drug_associations$color = 'gray90'
    idx = which( TF_drug_associations[,fdrcol] < th & TF_drug_associations$L1 %in% tissueTCGA_hash$tissue)
    TF_drug_associations$color[ idx ] = tissueTCGA_hash$color[ match(TF_drug_associations$L1[ idx  ], tissueTCGA_hash$tissue) ]
    TF_drug_associations$labels = ''
    idx = 1:repel_n[1]
    if( length(repel_n) > 1 )
      idx = repel_n
    TF_drug_associations$labels[idx] = paste("", TF_drug_associations$L1[idx], "[", TF_drug_associations$Drug_name[idx], "]^",  gsub(' ', '', TF_drug_associations$TF[idx]),  sep='')
    p = ggplot(TF_drug_associations[ TF_drug_associations$X_lm_pval < 0.75 & abs(TF_drug_associations$effect) < 7.5 , ], aes(x=effect, y=log10pval, color = color, size = (log10pval^2)*.01)) + 
      geom_point(alpha = .65) + 
      geom_point(data = TF_drug_associations[ TF_drug_associations[, fdrcol] < th & abs(TF_drug_associations$effect) < 7.5 & ! TF_drug_associations$L1 %in% c('SCLC', 'OV', 'COREAD', 'BRCA'), ], alpha = .65) + 
      geom_hline(yintercept = -log10(th), colour="gray20", linetype = "longdash") +
      annotate('text', x =-Inf,  y = -log10(th), label = paste('FDR <', th), colour="gray20", hjust = -.1, vjust = -.5, size = 4.5) +
      scale_colour_identity()+
      guides(size=FALSE) +
      scale_y_continuous(breaks=integer_breaks) +
      xlab('effect size') + ylab(expression(-log[10]~adj.pvalue)) +
      mytheme + 
      ggtitle(main)
    p = p + geom_label_repel(data=TF_drug_associations[ idx, ], aes(label = labels), color = 'black', force = 0.5, 
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

plot_heatmap = function(TF_drug_associations, main = '', path = '', lm = T, th = .01){
  if( ! is.null(TF_drug_associations) ) {
    TF_drug_associations = TF_drug_associations[TF_drug_associations[ fdrcol ] < th, ]
    if( nrow(TF_drug_associations) > 1){
      if (lm)
        TF_drug_associations$value = TF_drug_associations$X_lm_coefficient
      if (! lm)
        TF_drug_associations$value = TF_drug_associations$FEATURE_IC50_effect_size
      TF_drug_associations$Drug_name = paste(TF_drug_associations$Drug_name, ' [', TF_drug_associations$Drug_targets, ']', sep = '')
      X = acast(TF_drug_associations, Drug_name~TF, fill = 0, value.var='value')
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



get_heatmap_TF_clusters = function(TF_drug_associations, main = '',  path = '', k, lm = T){
  if( ! is.null(TF_drug_associations) ) {
    TF_drug_associations = TF_drug_associations[TF_drug_associations[ fdrcol ] < 0.01, ]
    if( nrow(TF_drug_associations) > 1){
      if (lm)
        TF_drug_associations$value = TF_drug_associations$X_lm_coefficient
      if (! lm)
        TF_drug_associations$value = TF_drug_associations$FEATURE_IC50_effect_size
      TF_drug_associations$Drug_name = paste(TF_drug_associations$Drug_name, ' [', TF_drug_associations$Drug_targets, ']', sep = '')
      X = acast(TF_drug_associations, Drug_name~TF, fill = 0, value.var='value')
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



