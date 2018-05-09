
source('code/utils.r')
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
    if( method == 'PR50'){
      observed = 1 - observed # revert order since the function assumes that the higher the value, the better prediction
      performance_value = get_PR_at_50(observed, expected)
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
    if( method == 'PR50'){
      # ac = sapply(r_negatives, function(ne) get_PR_at_50(c(ne, positives),  c(rep(0, n_positives), rep(1, n_positives) ) )  )
      ac = mapply(function(ne, po) get_PR_at_50(c(ne, po),  c(rep(0, n), rep(1, n) ) ), r_negatives, r_positives)
    }
    performance_value = mean(unlist(ac))
  }
  return(performance_value)
}


get_aucPR = function(observed, expected){
  aucPR = pr.curve(scores.class0 = observed[ expected == 0], scores.class1 = observed[ expected == 1])$auc.integral
  return(aucPR)
}

get_PR_at_50 = function(observed, expected){
  mypred = prediction(observed, expected==1)
  myperf = performance(mypred, 'prec', 'rec')
  # plot(myperf)
  precision = myperf@y.values[[1]]
  recall = myperf@x.values[[1]]
  pre_at_05 = max(precision[ recall == 0.5 ])
  return(pre_at_05)
}


get_aucROC = function(observed, expected){
  myroc = roc(predictor = observed, response = expected, smooth=F)
  performance_value = myroc$auc[1]
  return(performance_value)
}





auc2coverage_plot = function(ddf, y='NES', wrap_variable='regulon_dataset', 
                             balanced_accuracy = F,
                             regulon_evidence_filter=NULL, filter_datasets = T, 
                             y_scale = NULL, x_scale = NULL, legend_position = 'inside'){
  ddf$wrap_variable = ddf[, wrap_variable]
  ddf$y = ddf[, y]
  if( ! is.null(regulon_evidence_filter) )
    ddf = subset(ddf, regulon_evidence %in% regulon_evidence_filter)
  
  wrap_variable_filter = names(which(table(unique(ddf[, c('wrap_variable', 'is_TF_perturbed') ])$wrap_variable) == 2))
  ddf = subset(ddf, wrap_variable %in% wrap_variable_filter)
  
  df_plot = ddply(ddf, c(wrap_variable, 'regulon_evidence', 'regulon_dataset'), function(x) {
    observed = x$y
    expected = (x$is_TF_perturbed + 0)
    performance_value = compute_accuracy(observed, expected, method = 'aucPR', balanced = balanced_accuracy)
    # compute_accuracy(observed, expected, method = 'PR')
    return(performance_value)
  })
  names(df_plot)[ names(df_plot) == 'V1' ]  = 'AUC'
  df_plot$coverage = ddply(ddf, c(wrap_variable, 'regulon_evidence', 'regulon_dataset'), function(x) {
    length(unique(x$TF[ x$is_TF_perturbed ]))
  })$V1
  df_plot$coverage[ grep('inferredGTEx ', df_plot$regulon_dataset, ignore.case = T) ] = length(unique(ddf$TF[ intersect(which(ddf$is_TF_perturbed), grep('inferredGTEX ', ddf$regulon_dataset, ignore.case = T) ) ])) # aggregate coverage from tissue-specific networks
  df_plot = df_plot[ order(df_plot$AUC, decreasing = F), ]
  
  # select which datasets to plot
  if( filter_datasets ){
      # df_plot = subset(df_plot, regulon_evidence != 'old_consensus')
      df_plot = subset(df_plot, ! regulon_dataset %in% c('TFe', 'tfact', 'trrust', 'oreganno', 'trrd_via_tfact', 'evidence3_signed', 'evidence4_signed')  )
      df_plot = df_plot[ grep('evidence2$', df_plot$regulon_dataset, invert = T) , ]
      df_plot = df_plot[ grep('^evidence[3-4]_', df_plot$regulon_dataset, invert = T) , ]
      # df_plot = df_plot[ grep('_curated', df_plot$regulon_dataset, invert = T) , ]
      df_plot = df_plot[ ! 1:nrow(df_plot) %in% intersect(grep('evidence', df_plot$regulon_dataset), grep('signed', df_plot$regulon_dataset))  , ]
  }
  # group datasets for the plot
  df_plot$group = manage_datasets_names(df_plot$regulon_dataset, what = 'dataset2group')
  df_plot$regulon_evidence[ df_plot$regulon_evidence == 'consensus' ] = 'consensus_between'
  df_plot$regulon_evidence[ grep('[1-3]curateddatabases', df_plot$regulon_dataset) ] = 'consensus_within'
  df_plot$group = gsub('evidence_', '', df_plot$group)
  # PLOT
  P = ggplot(df_plot, aes(y=coverage, x=AUC, shape = regulon_evidence, color = regulon_evidence, group = group)) +
    geom_vline(xintercept = 0.5) + geom_vline(xintercept = 0.475, linetype = 'dashed') + geom_vline(xintercept = 0.525, linetype = 'dashed')+
    geom_label_repel(data = df_plot[ ! duplicated(df_plot$group) , ], aes(label = group ), show.legend = F, segment.colour = 'gray30', size = 3) + 
    geom_point(size = 3) + geom_line(show.legend = F) +
    scale_colour_manual(values = manage_datasets_names(sort(unique(df_plot$regulon_evidence)), what = 'evidence2color'), name = '') + 
    scale_shape_manual(values = manage_datasets_names(sort(unique(df_plot$regulon_evidence)), what = 'evidence2shape'), name = '') +
    mytheme + ylab('Coverage\n(number of TFs covered)') + xlab('Area under the PR curve (AUC)')
  if (legend_position == 'inside')  
    P = P + theme(legend.position = c(1,1), legend.justification = c(1.05, 1.05), legend.background = element_rect(size=0.5, linetype="solid", colour ="gray"), legend.text = element_text(size = 9), legend.title = element_blank())
  if (legend_position == 'bottom')  
    P = P +theme(legend.position = 'bottom', legend.background = element_rect(size=0.5, linetype="solid", colour ="gray"), legend.text = element_text(size = 9), legend.title = element_text(size = 10, face = 'bold'))
  if( ! is.null(y_scale) )
    P = P + scale_y_continuous(limits = y_scale)
  if( ! is.null(x_scale) )
    P = P + scale_x_continuous(limits = x_scale)
  return(list(P=P, aucs = df_plot))
}






plot_PR = function(ddf, y='NES', wrap_variable='regulon_dataset', 
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
  
  df_plot = dlply(ddf, c(wrap_variable, 'regulon_evidence'), function(x) {
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
  df_plot = melt(df_plot)
  df_plot[[wrap_variable]] = sapply(strsplit(df_plot$L1, split = "\\."), head, 1)
  df_plot$regulon_evidence = sapply(strsplit(df_plot$L1, split = "\\."), tail, 1)
  df_plot$wrap_variable = df_plot[, wrap_variable]
  
  P = ggplot(df_plot, aes(x=wrap_variable, y=value, color = wrap_variable)) + 
    geom_hline(yintercept = 0.5) +
    geom_boxplot(fill = 'white', alpha = 0, outlier.size = NA, outlier.color = NA, outlier.shape = NA) + 
    geom_quasirandom(size = 0.75, alpha =.3) +
    scale_y_continuous(limits = c(0.4, 1) ) +
    coord_flip() +
    scale_color_manual(values = line_colors, name = '') +
    geom_abline(intercept=1, slope=1, linetype="dashed") + xlab('') + ylab("Area under the PR curve (AUC)") + 
    mytheme + theme(legend.position = 'none')
  P
  return(P)
}


plot_rocs = function(ddf, y='NES', wrap_variable='regulon_dataset', 
                     balanced_accuracy = F,
                     regulon_evidence_filter=NULL, TFs_filter = NULL, 
                     line_colors = my_color_palette$EMBL){
  library(pROC)
  library(plotROC)
  ddf$wrap_variable = ddf[, wrap_variable]
  ddf$y = ddf[, y]
  if( ! is.null(regulon_evidence_filter) )
    ddf = subset(ddf, regulon_evidence %in% regulon_evidence_filter)
  if( ! is.null(TFs_filter) )
    ddf = subset(ddf, TF %in% TFs_filter)
  wrap_variable_filter = names(which(table(unique(ddf[, c('wrap_variable', 'is_TF_perturbed') ])$wrap_variable) == 2))
  ddf = subset(ddf, wrap_variable %in% wrap_variable_filter)
  
  df_plot = ddply(ddf, c(wrap_variable, 'regulon_evidence'), function(x) {
    observed = x$y
    expected = (x$is_TF_perturbed + 0)
    myroc = NULL
    try(myroc <- roc(predictor = observed, response = expected, smooth=T), silent = T)
    if( is.null(myroc) )
      myroc = roc(predictor = observed, response = expected, smooth=F)
    data.frame(specificity=myroc$specificities,
               sensitivity=myroc$sensitivities, 
               auc = compute_accuracy(observed, expected, method = 'aucROC', balanced = balanced_accuracy))
  })
  df_plot$wrap_variable = df_plot[, wrap_variable]
  df_plot$wrap_variable = paste(df_plot$wrap_variable, ' (auc = ', round(df_plot$auc, 2), ')', sep = '')
  
  P = ggplot(df_plot, aes(x=specificity, y=sensitivity, color = wrap_variable, group = wrap_variable)) + 
    geom_path() + 
    # geom_point(shape=3)  + 
    scale_x_reverse() + 
    scale_color_manual(values = line_colors, name = '') +
    geom_abline(intercept=1, slope=1, linetype="dashed") + xlab("1 - specificity") + ylab("sensitivity") + 
    mytheme + theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.text = element_text(size = 10))
  return(P)
}


plot_boxplots = function(ddf, y='NES', wrap_variable='regulon_dataset', regulon_evidence_filter=NULL, TFs_filter = NULL, p2_with_names = F){
  ddf$wrap_variable = ddf[, wrap_variable]
  ddf$y = ddf[, y]
  if( ! is.null(regulon_evidence_filter) )
    ddf = subset(ddf, regulon_evidence %in% regulon_evidence_filter)
  if( ! is.null(TFs_filter) )
    ddf = subset(ddf, TF %in% TFs_filter)
  wrap_variable_filter = names(which(table(unique(ddf[, c('wrap_variable', 'is_TF_perturbed') ])$wrap_variable) == 2))
  ddf = subset(ddf, wrap_variable %in% wrap_variable_filter)
  
  mddf = ddply(ddf, c('wrap_variable', 'regulon_evidence'), function(x){
    # message(x[1,])
    pval = NA
    try(pval <- t.test(subset(x, is_TF_perturbed)$y, subset(x, ! is_TF_perturbed)$y, alternative = 'less')$p.value, silent = T)
    return(pval)
    })
  names(mddf)[3] = 'pval'
  mddf$pval_lab = paste('p =', signif(mddf$pval, 3))
  mddf$effect_size = ddply(ddf, c('wrap_variable', 'regulon_evidence'), function(x) effectsize_cohensD(subset(x, is_TF_perturbed)$y, subset(x, ! is_TF_perturbed)$y) )$V1
  mddf$median = ddply(ddf, c('wrap_variable', 'regulon_evidence'), function(x) median(subset(x, is_TF_perturbed)$y, na.rm = T) )$V1
  mddf = mddf[ order(mddf$median), ]
  ddf$wrap_variable = factor(ddf$wrap_variable, levels = mddf$wrap_variable)
  mddf$wrap_variable = factor(mddf$wrap_variable, levels = mddf$wrap_variable)
  p1 = ggplot(ddf, aes(x=is_TF_perturbed, y = y, color = regulon_evidence)) + geom_boxplot() + 
    facet_wrap(~wrap_variable , ncol = 8) + 
    geom_text(data = mddf, aes(label=mddf$pval_lab, y= Inf, x= Inf), vjust = 1, hjust = 1, color = 'red') + mytheme + ylab(y)
  ddf$coverage = 'coverage'
  p2 = ggplot(subset(ddf, is_TF_perturbed), aes(x=wrap_variable, y = y, color = regulon_evidence, label = TF)) + 
    stat_boxplot(geom = "errorbar", width = 0.5) + 
    geom_boxplot(fill = 'white', alpha = 1, outlier.size = NA, outlier.color = NA, outlier.shape = NA) +
    geom_quasirandom(size = 1.5, alpha =.5, color = 'gray50') + 
    geom_boxplot(alpha = .2, outlier.size = NA, outlier.color = NA, outlier.shape = NA) +
    geom_hline(yintercept = 0.5, color = 'red') + coord_flip() + mytheme + ylab(y) + xlab(wrap_variable)  + scale_color_manual(values = my_color_palette$shiny)
  if( p2_with_names )
    p2 = p2 + geom_text_repel()
  p3 = ggplot(subset(ddf, is_TF_perturbed), aes(x=wrap_variable, fill = regulon_evidence)) + geom_bar() + 
    coord_flip() + mytheme + ylab('perturbations') + xlab('dataset') + scale_fill_manual(values = my_color_palette$shiny)
  # p4 = ggplot(mddf, aes(x=wrap_variable, fill = effect_size, y = -log10(pval))) + geom_bar(stat = 'identity') + facet_wrap(~regulon_evidence , ncol = 8) + 
    # coord_flip() + mytheme + xlab('dataset') + geom_hline(yintercept = -log10(0.05), color = 'red') + theme(legend.key.width = unit(1.5, "cm"))
  return(list(p1, p2, p3))
}

