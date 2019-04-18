create_contingency_table = function(genes, fungenes, background_genes){
  nonfungenes = setdiff(background_genes, fungenes)
  regulon_fungenes = intersect(fungenes, genes)
  regulon_nonfungenes = intersect(nonfungenes, genes)
  m = matrix(c(length(regulon_fungenes), length(fungenes),
               length(regulon_nonfungenes), length(nonfungenes)), ncol = 2, dimnames = list(c('in_geneset', 'out_geneset'), c('regulon', 'genome')) )
  return(m)  
}

enrichment = function(fun, genes, background_genes){
  funtable = create_contingency_table(genes, fun, background_genes)
  out = NULL
  if ( funtable[1,1] > 1 ){
    ft = fisher.test(funtable)
    pvalue = ft$p.value
    estimate = ft$estimate
    conf.int = ft$conf.int
    out = data.frame(pvalue= pvalue, odds.ratio = estimate, min_confint = conf.int[1], max_confint = conf.int[2],  stringsAsFactors = F)
  return(out)
  }
}

analyse_genesets = function(genes, geneset){
  if( length(genes) >= 5 ){
    message('provided ', length(genes), ' genes')
    background_genes = unique(unlist(geneset))
    genes = intersect(genes, background_genes)
    message('analysing ', length(genes), ' genes')
    tftests = lapply(geneset, enrichment, genes, background_genes)
    names(tftests) = names(geneset)
    tftests = Filter(Negate(is.null), tftests)
    if( length(tftests) > 0 ){
      mtftests = melt(tftests, id.vars = names(tftests[[1]]) )
      names(mtftests) = c(names(tftests[[1]]), 'geneset')
      mtftests$fdr = p.adjust(mtftests$pvalue, method = 'fdr')
      mtftests = mtftests[ order(mtftests$fdr), ]
      return(mtftests)
    }
  }
}

plot_enrichment = function(re, feature = 'TF_class'){
  if( ! is.null(feature) ){
    mdf = melt(re, id.vars = names(re[[1]]) )
    mdf = subset(mdf, L1 == feature)
  }else{
    mdf = melt(re, id.vars = names(re) )
  }
  mdf$log.odds.ratio = log(mdf$odds.ratio)
  mdf$geneset = gsub('_', ' ', mdf$geneset)
  mdf$geneset = factor(mdf$geneset, levels = rev(mdf$geneset))
  ggplot(subset(mdf, fdr < 0.05), aes(x = geneset, y= log.odds.ratio, color = fdr) ) +
    geom_point(stat = 'Identity', size = 5) +
    geom_errorbar(aes(ymin=log(min_confint), ymax=log(max_confint)), position = position_dodge(width = .5), width=0.3, size=0.3, color = 'black') +
    geom_hline(yintercept = 0) +
    scale_color_gradient(low = brewer.pal('RdBu', n = 8)[8], high = brewer.pal('RdBu', n = 8)[5], name = 'FDR') +
    mytheme + coord_flip() + 
    xlab('') + ylab('log( odds ratio )')  +
    theme(legend.key.width=unit(1.5,"cm"), legend.key.height = unit(0.3,"cm"))
}


plot_enrichment_grid = function(re){
  mdf = melt(re, id.vars = names(re[[1]]) )
  mdf$log.odds.ratio = log(mdf$odds.ratio)
  mdf$geneset = gsub('_', ' ', mdf$geneset)
  mdf$geneset
  mdf$L1 = gsub('at_least', '\u2265', mdf$L1)
  mdf$L1 = gsub('_', ' ', mdf$L1)
  ggplot(subset(mdf, fdr < 0.05), aes(x = geneset, y= log.odds.ratio, color = fdr) ) +
    geom_point(stat = 'Identity', size = 3) +
    geom_errorbar(aes(ymin=log(min_confint), ymax=log(max_confint)), position = position_dodge(width = .5), width=0.3, size=0.3, color = 'black') +
    geom_hline(yintercept = 0) +
    scale_color_gradient(low = brewer.pal('RdBu', n = 8)[8], high = brewer.pal('RdBu', n = 8)[5], name = 'FDR') +
    mytheme + coord_flip() + 
    xlab('') + ylab('log( odds ratio )')  +
    theme(legend.key.width=unit(1.5,"cm"), legend.key.height = unit(0.3,"cm")) +
    facet_wrap(~L1, nrow = 1)
}
