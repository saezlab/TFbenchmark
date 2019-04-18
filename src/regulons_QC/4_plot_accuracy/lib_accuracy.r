
require(pROC)
require(ROCR)
require(PRROC)




compute_accuracy = function(observed, expected, method = 'aucPR', balanced = F){
  if(balanced){
    if( method == 'aucROC'){
      performance_value = get_aucROC(observed, expected)
    }
    if( method == 'aucPR'){
      performance_value = get_aucPR(observed, expected)
    }
  }else{
    n_positives = sum(expected==1) 
    n_negatives = sum(expected==0) 
    positives = observed[ expected == 1 ]
    negatives = observed[ expected == 0 ]
    n = min(n_positives, n_negatives)
    r_negatives = lapply(1:100, function(ra) sample(negatives, n, replace = F)  ) # down-sample the negatives to balance 
    r_positives = lapply(1:100, function(ra) sample(positives, n, replace = F)  ) # down-sample the positives to balance 
    if( method == 'aucROC'){
      ac = mapply(function(ne, po) get_aucROC(c(ne, po),  c(rep(0, n), rep(1, n) ) ), r_negatives, r_positives)
    }
    if( method == 'aucPR'){
      # ac = sapply(r_negatives, function(ne) get_aucPR(c(ne, positives),  c(rep(0, n_positives), rep(1, n_positives) ) )  )
      ac = mapply(function(ne, po) get_aucPR(c(ne, po),  c(rep(0, n), rep(1, n) ) ), r_negatives, r_positives)
    }
    performance_value = mean(unlist(ac))
  }
  return(performance_value)
}


get_aucPR = function(observed, expected){
  aucPR = pr.curve(scores.class0 = observed[ expected == 0], scores.class1 = observed[ expected == 1])$auc.integral
  return(aucPR)
}



get_aucROC = function(observed, expected){
  myroc = roc(predictor = observed, response = expected, smooth=F)
  performance_value = myroc$auc[1]
  return(performance_value)
}



activities2accuracy = function(ddf, y='NES', wrap_variable='regulon_dataset', 
                             balanced_accuracy = F, performance_method = 'aucPR',
                             regulon_evidence_filter=NULL){
  ddf$wrap_variable = ddf[, wrap_variable]
  ddf$y = ddf[, y]
  if( ! is.null(regulon_evidence_filter) )
    ddf = subset(ddf, regulon_evidence %in% regulon_evidence_filter)
  
  wrap_variable_filter = names(which(table(unique(ddf[, c('wrap_variable', 'is_TF_perturbed') ])$wrap_variable) == 2))
  ddf = subset(ddf, wrap_variable %in% wrap_variable_filter)
  
  # Compute performance ~ AUC
  df_coverage2accuracy = ddply(ddf, c(wrap_variable, 'regulon_evidence', 'regulon_dataset'), function(x) {
    observed = x$y
    expected = (x$is_TF_perturbed + 0)
    performance_value = compute_accuracy(observed, expected, method = performance_method, balanced = balanced_accuracy)
    return(performance_value)
  })
  # Format dataframe
  names(df_coverage2accuracy)[ names(df_coverage2accuracy) == 'V1' ]  = 'AUC'
  df_coverage2accuracy$coverage = ddply(ddf, c(wrap_variable, 'regulon_evidence', 'regulon_dataset'), function(x) {
    length(unique(x$TF[ x$is_TF_perturbed ]))
  })$V1
  # df_coverage2accuracy$coverage[ grep('inferredGTEx ', df_coverage2accuracy$regulon_dataset, ignore.case = T) ] = length(unique(ddf$TF[ intersect(which(ddf$is_TF_perturbed), grep('inferredGTEX ', ddf$regulon_dataset, ignore.case = T) ) ])) # aggregate coverage from tissue-specific networks
  df_coverage2accuracy = df_coverage2accuracy[ order(df_coverage2accuracy$AUC, decreasing = T), ]
  return(df_coverage2accuracy)
}



activities2accuracy_format = function(df_coverage2accuracy, filter_datasets =T){
  # Select which datasets to plot
  if( filter_datasets ){
    df_coverage2accuracy = subset(df_coverage2accuracy, ! regulon_dataset %in% c('TFe', 'tfact', 'trrust', 'oreganno', 'trrd_via_tfact')  )
    # df_coverage2accuracy = df_coverage2accuracy[ grep('evidence2$', df_coverage2accuracy$regulon_dataset, invert = T) , ]
    # df_coverage2accuracy = df_coverage2accuracy[ grep('^evidence[3-4]_', df_coverage2accuracy$regulon_dataset, invert = T) , ]
    # df_coverage2accuracy = df_coverage2accuracy[ ! 1:nrow(df_coverage2accuracy) %in% intersect(grep('evidence', df_coverage2accuracy$regulon_dataset), grep('signed', df_coverage2accuracy$regulon_dataset))  , ]
  }
  # group datasets for the plot
  df_coverage2accuracy$group = manage_datasets_names(df_coverage2accuracy$regulon_dataset, what = 'dataset2group')
  df_coverage2accuracy$regulon_evidence[ df_coverage2accuracy$regulon_evidence == 'consensus' ] = 'consensus_all'
  df_coverage2accuracy$regulon_evidence[ grep('[1-3]curateddatabases', df_coverage2accuracy$regulon_dataset) ] = 'consensus_curated'
  # prettify labels for publication purpose
  df_coverage2accuracy$group[ grep('[1-3]curateddatabases', df_coverage2accuracy$regulon_dataset) ] = gsub('curateddatabases', '', df_coverage2accuracy$group[ grep('[1-3]curateddatabases', df_coverage2accuracy$regulon_dataset) ])
  df_coverage2accuracy$group = gsub('e3', 'any_3', df_coverage2accuracy$group) %>% 
    gsub('e2', 'any_2',  .)  %>% 
    gsub('e4', 'any_4',  .)
  df_coverage2accuracy$group = gsub('_$', '', df_coverage2accuracy$group) %>% 
    gsub('_signed', '*', .) %>% 
    gsub('databases', '', .) %>% 
    gsub('_scanning', '', .)  %>% 
    gsub('scanning', '', .) 
  df_coverage2accuracy$group = sapply(strsplit(df_coverage2accuracy$group, '_via_'), head, 1)
  df_coverage2accuracy$group = gsub('_PANCANCER$', ' ', df_coverage2accuracy$group)
  return(df_coverage2accuracy)
}



auc2coverage_plot = function(df_coverage2accuracy, 
                             y_scale = NULL, x_scale = NULL, legend_position = 'inside'){
  P = ggplot(df_coverage2accuracy, aes(y=coverage, x=AUC, shape = regulon_evidence, color = regulon_evidence, group = group)) +
    geom_vline(xintercept = 0.5) + geom_vline(xintercept = 0.475, linetype = 'dashed') + geom_vline(xintercept = 0.525, linetype = 'dashed')+
    geom_point(size = 3) + geom_line(show.legend = F) +
    geom_label_repel(data = df_coverage2accuracy[ ! duplicated(df_coverage2accuracy$group) , ], aes(label = group ), show.legend = F, segment.colour = 'gray30', size = 3.5, force = 10, parse = F) + 
    scale_colour_manual(values = manage_datasets_names(sort(unique(df_coverage2accuracy$regulon_evidence)), what = 'evidence2color'), name = '') + 
    scale_shape_manual(values = manage_datasets_names(sort(unique(df_coverage2accuracy$regulon_evidence)), what = 'evidence2shape'), name = '') +
    mytheme + ylab('Coverage\n(number of TFs in the benchmark)') + xlab('Area under the PR curve (AUC)')
  if (legend_position == 'inside')  
    P = P + theme(legend.position = c(1,1), legend.justification = c(1.05, 1.05), legend.background = element_rect(size=0.5, linetype="solid", colour ="gray"), legend.text = element_text(size = 9), legend.title = element_blank())
  if (legend_position == 'bottom')  
    P = P +theme(legend.position = 'bottom', legend.background = element_rect(size=0.5, linetype="solid", colour ="gray"), legend.text = element_text(size = 9), legend.title = element_text(size = 10, face = 'bold'))
  if( ! is.null(y_scale) )
    P = P + scale_y_continuous(limits = y_scale)
  if( ! is.null(x_scale) )
    P = P + scale_x_continuous(limits = x_scale)
  return(P)
}




plot_PR = function(ddf, y='NES', wrap_variable='regulon_dataset', 
                     balanced_accuracy = F,
                     regulon_evidence_filter=NULL, TFs_filter = NULL, 
                     line_colors = my_color_palette$EMBL){
  ddf$wrap_variable = ddf[, wrap_variable]
  ddf$y = ddf[, y]
  if( ! is.null(regulon_evidence_filter) )
    ddf = subset(ddf, regulon_evidence %in% regulon_evidence_filter)
  if( ! is.null(TFs_filter) ){
    ddf = subset(ddf, TF %in% TFs_filter)
  }
  if( length(line_colors) == 1)
    line_colors = rep(line_colors, 100)
  wrap_variable_filter = names(which(table(unique(ddf[, c('wrap_variable', 'is_TF_perturbed') ])$wrap_variable) == 2))
  ddf = subset(ddf, wrap_variable %in% wrap_variable_filter)
  
  df_accuracy = dlply(ddf, c(wrap_variable, 'regulon_evidence'), function(x) {
    observed = x$y
    expected = (x$is_TF_perturbed + 0)
    if(balanced_accuracy){
        auc = get_aucPR(observed, expected)
    }else{
      n_positives = sum(expected==1) 
      n_negatives = sum(expected==0) 
      positives = observed[ expected == 1 ]
      negatives = observed[ expected == 0 ]
      n = min(n_positives, n_negatives)
      r_negatives = lapply(1:100, function(ra) sample(negatives, n, replace = F)  ) # down-sample the negatives to balance 
      r_positives = lapply(1:100, function(ra) sample(positives, n, replace = F)  ) # down-sample the positives to balance 
      auc = mapply(function(ne, po) get_aucPR(c(ne, po),  c(rep(0, n), rep(1, n) ) ), r_negatives, r_positives)
    }
    return(auc)
  })
  df_accuracy = melt(df_accuracy)
  df_accuracy[[wrap_variable]] = sapply(strsplit(df_accuracy$L1, split = "\\."), head, 1)
  df_accuracy$regulon_evidence = sapply(strsplit(df_accuracy$L1, split = "\\."), tail, 1)
  df_accuracy$wrap_variable = df_accuracy[, wrap_variable]
  
  P = ggplot(df_accuracy, aes(x=wrap_variable, y=value, color = wrap_variable)) + 
    geom_hline(yintercept = 0.5) +
    geom_boxplot(fill = 'white', alpha = 0, outlier.size = NA, outlier.color = NA, outlier.shape = NA) + 
    geom_quasirandom(size = 0.75, alpha =.3) +
    scale_y_continuous(limits = c(0.4, 1) ) +
    coord_flip() +
    scale_color_manual(values = line_colors, name = '') +
    geom_abline(intercept=1, slope=1, linetype="dashed") + xlab('') + ylab("Area under the PR curve (AUC)") + 
    mytheme + theme(legend.position = 'none')
  if( ! is.null(TFs_filter) )
    P = P + annotate('text', x = Inf, y = Inf, label = paste('n =', length(unique(ddf$TF))), vjust = 1.2, hjust = 1.2)
  return(P)
}



plot_rocs = function(ddf, y='NES', wrap_variable='regulon_dataset', 
                   balanced_accuracy = F,
                   regulon_evidence_filter=NULL, TFs_filter = NULL, 
                   line_colors = my_color_palette$EMBL){
  ddf$wrap_variable = ddf[, wrap_variable]
  ddf$y = ddf[, y]
  if( ! is.null(regulon_evidence_filter) )
    ddf = subset(ddf, regulon_evidence %in% regulon_evidence_filter)
  if( ! is.null(TFs_filter) )
    ddf = subset(ddf, TF %in% TFs_filter)
  if( length(line_colors) == 1)
    line_colors = rep(line_colors, 100)
  wrap_variable_filter = names(which(table(unique(ddf[, c('wrap_variable', 'is_TF_perturbed') ])$wrap_variable) == 2))
  ddf = subset(ddf, wrap_variable %in% wrap_variable_filter)
  
  df_accuracy = dlply(ddf, c(wrap_variable, 'regulon_evidence'), function(x) {
    observed = x$y
    expected = (x$is_TF_perturbed + 0)
    if(balanced_accuracy){
      auc = get_aucROC(observed, expected)
    }else{
      n_positives = sum(expected==1) 
      n_negatives = sum(expected==0) 
      positives = observed[ expected == 1 ]
      negatives = observed[ expected == 0 ]
      n = min(n_positives, n_negatives)
      r_negatives = lapply(1:100, function(ra) sample(negatives, n, replace = F)  ) # down-sample the negatives to balance 
      r_positives = lapply(1:100, function(ra) sample(positives, n, replace = F)  ) # down-sample the positives to balance 
      auc = mapply(function(ne, po) get_aucPR(c(ne, po),  c(rep(0, n), rep(1, n) ) ), r_negatives, r_positives)
    }
    return(auc)
  })
  df_accuracy = melt(df_accuracy)
  df_accuracy[[wrap_variable]] = sapply(strsplit(df_accuracy$L1, split = "\\."), head, 1)
  df_accuracy$regulon_evidence = sapply(strsplit(df_accuracy$L1, split = "\\."), tail, 1)
  df_accuracy$wrap_variable = df_accuracy[, wrap_variable]
  
  P = ggplot(df_accuracy, aes(x=wrap_variable, y=value, color = wrap_variable)) + 
    geom_hline(yintercept = 0.5) +
    geom_boxplot(fill = 'white', alpha = 0, outlier.size = NA, outlier.color = NA, outlier.shape = NA) + 
    geom_quasirandom(size = 0.75, alpha =.3) +
    scale_y_continuous(limits = c(0.4, 1) ) +
    coord_flip() +
    scale_color_manual(values = line_colors, name = '') +
    geom_abline(intercept=1, slope=1, linetype="dashed") + xlab('') + ylab("Area under the ROC curve (AUC)") + 
    mytheme + theme(legend.position = 'none')
  if( ! is.null(TFs_filter) )
    P = P + annotate('text', x = Inf, y = Inf, label = paste('n =', length(unique(ddf$TF))), vjust = 1.2, hjust = 1.2)
  P
  return(P)
}






############ REVISSION specific plots

auc2coverage_plot_oncology = function(df_accuracy,
                             y_scale = NULL, x_scale = NULL, legend_position = 'inside'){
  
  P = ggplot(df_accuracy, aes(y=coverage, x=AUC, shape = group, color = group, group = group)) +
    geom_vline(xintercept = 0.5) + geom_vline(xintercept = 0.475, linetype = 'dashed') + geom_vline(xintercept = 0.525, linetype = 'dashed')+
    geom_point(size = 3) + 
    geom_line(show.legend = F) +
    geom_label_repel(data = df_accuracy[ ! duplicated(df_accuracy$group) , ], aes(label = group ), show.legend = F, segment.colour = 'gray30', size = 3.5, force = 10) + 
    scale_colour_manual(values = brewer.pal(n = 6, name = 'Paired')[c(5,3,4,6)]) + 
    mytheme + ylab('Coverage\n(number of TFs in benchmark)') + xlab('Area under the PR curve (AUC)')
  if (legend_position == 'inside')  
    P = P + theme(legend.position = c(1,1), legend.justification = c(1.05, 1.05), legend.background = element_rect(size=0.5, linetype="solid", colour ="gray"), legend.text = element_text(size = 9), legend.title = element_blank())
  if (legend_position == 'bottom')  
    P = P +theme(legend.position = 'bottom', legend.background = element_rect(size=0.5, linetype="solid", colour ="gray"), legend.text = element_text(size = 9), legend.title = element_text(size = 10, face = 'bold'))
  if( ! is.null(y_scale) )
    P = P + scale_y_continuous(limits = y_scale)
  if( ! is.null(x_scale) )
    P = P + scale_x_continuous(limits = x_scale)
  return(P)
}




plot_PR_compareCOMBAT = function(ddf, y='NES', wrap_variable='regulon_dataset', 
                                 balanced_accuracy = F,
                                 TFs_filter = NULL, 
                                 line_colors = my_color_palette$EMBL){
  ddf$wrap_variable = ddf[, wrap_variable]
  ddf$y = ddf[, y]
  regulon_evidence_filter = 'inferred'
  ddf = subset(ddf, regulon_evidence %in% regulon_evidence_filter)
  if( ! is.null(TFs_filter) ){
    ddf = subset(ddf, TF %in% TFs_filter)
    TF_perturbations = names(which(table(ddf$perturbed_TF[ ddf$is_TF_perturbed ]) >= length(unique(ddf$wrap_variable))))
    ddf = subset(ddf, perturbed_TF %in% TF_perturbations)
  }
  if( length(line_colors) == 1)
    line_colors = rep(line_colors, 100)
  wrap_variable_filter = names(which(table(unique(ddf[, c('wrap_variable', 'is_TF_perturbed') ])$wrap_variable) == 2))
  ddf = subset(ddf, wrap_variable %in% wrap_variable_filter)
  
  df_accuracy = dlply(ddf, c(wrap_variable, 'regulon_evidence'), function(x) {
    observed = x$y
    expected = (x$is_TF_perturbed + 0)
    if(balanced_accuracy){
      auc = get_aucPR(observed, expected)
    }else{
      n_positives = sum(expected==1) 
      n_negatives = sum(expected==0) 
      positives = observed[ expected == 1 ]
      negatives = observed[ expected == 0 ]
      n = min(n_positives, n_negatives)
      r_negatives = lapply(1:100, function(ra) sample(negatives, n, replace = F)  ) # down-sample the negatives to balance 
      r_positives = lapply(1:100, function(ra) sample(positives, n, replace = F)  ) # down-sample the positives to balance 
      auc = mapply(function(ne, po) get_aucPR(c(ne, po),  c(rep(0, n), rep(1, n) ) ), r_negatives, r_positives)
    }
    return(auc)
  })
  df_accuracy = melt(df_accuracy)
  df_accuracy[[wrap_variable]] = sapply(strsplit(df_accuracy$L1, split = "\\."), head, 1)
  df_accuracy$regulon_evidence = sapply(strsplit(df_accuracy$L1, split = "\\."), tail, 1)
  df_accuracy$wrap_variable = df_accuracy[, wrap_variable]
  
  # run ttest
  df_accuracy$ttest_group = gsub('_noCOMBAT', '', df_accuracy$wrap_variable)
  df_accuracy$ttest_group = gsub('_COMBAT', '', df_accuracy$ttest_group)
  df_accuracy$ttest_group = gsub(' noCOMBAT', '', df_accuracy$ttest_group)
  df_accuracy$ttest_group = gsub(' COMBAT', '', df_accuracy$ttest_group)
  df_accuracy$ttest_group = gsub('GTEx_', '', df_accuracy$ttest_group)
  df_accuracy$ttest_group = gsub('GTEx ', '', df_accuracy$ttest_group)
  df_accuracy_ttest = ddply(df_accuracy, 'ttest_group', function(x) {
    x1 = x$value[  grep('noCOMBAT', x$wrap_variable, invert = F) ]
    x2 = x$value[  grep('noCOMBAT', x$wrap_variable, invert = T) ]
    ttest_pvalue = t.test(x1, x2)$p.value
  })
  df_accuracy_ttest$label = signif(df_accuracy_ttest$V1, digits = 2)
  df_accuracy_ttest$label = paste('p =', df_accuracy_ttest$label)
  df_accuracy_ttest$label[ df_accuracy_ttest$V1 < 0.001 ] = 'p < 0.001'
  
  
  df_accuracy$wrap_variable[  grep('noCOMBAT', df_accuracy$wrap_variable) ] = 'noCOMBAT'
  df_accuracy$wrap_variable[  grep('noCOMBAT', df_accuracy$wrap_variable, invert = T) ] = 'COMBAT'
  P = ggplot(df_accuracy, aes(x=wrap_variable, y=value, color = wrap_variable)) + 
    geom_hline(yintercept = 0.5) +
    geom_boxplot(fill = 'white', alpha = 0, outlier.size = NA, outlier.color = NA, outlier.shape = NA) + 
    geom_quasirandom(size = 0.75, alpha =.3) +
    scale_y_continuous(limits = c(0.4, 1) ) +
    coord_flip() +
    scale_color_manual(values = line_colors, name = '') +
    facet_grid(ttest_group~.) + 
    geom_abline(intercept=1, slope=1, linetype="dashed") + xlab('') + ylab("Area under the PR curve (AUC)") + 
    mytheme + theme(legend.position = 'none')
  P = P + geom_text(data = df_accuracy_ttest, aes(label = label), color = 'red', y = 0.9, x = 1.5)
  return(P)
}




plot_PR_compareCANCER = function(ddf, y='NES', wrap_variable='regulon_dataset', 
                                 balanced_accuracy = F,
                                 TFs_filter = NULL, 
                                 line_colors = my_color_palette$EMBL){
  ddf$wrap_variable = ddf[, wrap_variable]
  ddf$y = ddf[, y]
  if( ! is.null(TFs_filter) ){
    ddf = subset(ddf, TF %in% TFs_filter)
    TF_perturbations = names(which(table(ddf$perturbed_TF[ ddf$is_TF_perturbed ]) >= length(unique(ddf$wrap_variable))))
    ddf = subset(ddf, perturbed_TF %in% TF_perturbations)
  }
  if( length(line_colors) == 1)
    line_colors = rep(line_colors, 100)
  wrap_variable_filter = names(which(table(unique(ddf[, c('wrap_variable', 'is_TF_perturbed') ])$wrap_variable) == 2))
  ddf = subset(ddf, wrap_variable %in% wrap_variable_filter)
  
  df_accuracy = dlply(ddf, c(wrap_variable, 'regulon_evidence'), function(x) {
    observed = x$y
    expected = (x$is_TF_perturbed + 0)
    if(balanced_accuracy){
      auc = get_aucPR(observed, expected)
    }else{
      n_positives = sum(expected==1) 
      n_negatives = sum(expected==0) 
      positives = observed[ expected == 1 ]
      negatives = observed[ expected == 0 ]
      n = min(n_positives, n_negatives)
      r_negatives = lapply(1:100, function(ra) sample(negatives, n, replace = F)  ) # down-sample the negatives to balance 
      r_positives = lapply(1:100, function(ra) sample(positives, n, replace = F)  ) # down-sample the positives to balance 
      auc = mapply(function(ne, po) get_aucPR(c(ne, po),  c(rep(0, n), rep(1, n) ) ), r_negatives, r_positives)
    }
    return(auc)
  })
  df_accuracy = melt(df_accuracy)
  df_accuracy[[wrap_variable]] = sapply(strsplit(df_accuracy$L1, split = "\\."), head, 1)
  df_accuracy$regulon_evidence = sapply(strsplit(df_accuracy$L1, split = "\\."), tail, 1)
  df_accuracy$wrap_variable = df_accuracy[, wrap_variable]
  
  # run ttest
  df_accuracy$ttest_group = 'tissue-specific'
  df_accuracy$ttest_group[ grep('i10', df_accuracy$L1) ] = 'i10'
  df_accuracy$ttest_group[ grep('i5', df_accuracy$L1) ] = 'i5'
  df_accuracy$ttest_group[ grep('i3', df_accuracy$L1) ] = 'i3'
  df_accuracy$ttest_group[ grep('i2', df_accuracy$L1) ] = 'i2'
  df_accuracy_ttest = ddply(df_accuracy, 'ttest_group', function(x) {
    x1 = x$value[  grep('GTEx', x$wrap_variable, invert = F) ]
    x2 = x$value[  grep('GTEx', x$wrap_variable, invert = T) ]
    ttest_pvalue = t.test(x1, x2)$p.value
  })
  df_accuracy_ttest$label = signif(df_accuracy_ttest$V1, digits = 2)
  df_accuracy_ttest$label = paste('p =', df_accuracy_ttest$label)
  df_accuracy_ttest$label[ df_accuracy_ttest$V1 < 0.001 ] = 'p < 0.001'
  
  
  df_accuracy$wrap_variable[  grep('GTEx', df_accuracy$wrap_variable) ] = 'GTEx'
  df_accuracy$wrap_variable[  grep('cancer', df_accuracy$wrap_variable) ] = 'TCGA'
  P = ggplot(df_accuracy, aes(x=wrap_variable, y=value, color = wrap_variable)) + 
    geom_hline(yintercept = 0.5) +
    geom_boxplot(fill = 'white', alpha = 0, outlier.size = NA, outlier.color = NA, outlier.shape = NA) + 
    geom_quasirandom(size = 0.75, alpha =.3) +
    scale_y_continuous(limits = c(0.4, 1) ) +
    coord_flip() +
    scale_color_manual(values = line_colors, name = '') +
    facet_grid(ttest_group~.) + 
    geom_abline(intercept=1, slope=1, linetype="dashed") + xlab('') + ylab("Area under the PR curve (AUC)") + 
    mytheme + theme(legend.position = 'none')
  P = P + geom_text(data = df_accuracy_ttest, aes(label = label), color = 'red', y = 0.45, x = 1.5)
  return(P)
}

