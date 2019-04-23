
source('code/utils.r')


TF_properties_enrichment = function(df, genesets){
  all_results = list()
  datasets2study = unique(df$regulon_dataset)
  for(dataset in datasets2study){
    message(dataset)
    
    perturbation_ranking = subset(df, regulon_dataset == dataset & TF == perturbed_TF)
    mySignature = perturbation_ranking$signature
    names(mySignature) = perturbation_ranking$TF
    mySignature = mySignature[order(mySignature, decreasing = T)]
    
    
    for ( feature in TFrole_features ) {
      feature_network = genesets[[feature]]
      
      # results = GSEAlm(Y = mySignature[ names(mySignature) %in% unlist(feature_network)  ], genesets = feature_network)
      results = GSEAlm(Y = mySignature, genesets = feature_network)
      if( ! is.null(results) ){
        results = results[ order(results$pval, decreasing = F) , ]
        results$padj = p.adjust(results$pval, method = 'fdr')
        results$regulon_evidence = perturbation_ranking$regulon_evidence[1]
        all_results[[feature]][[dataset]] = results
      }
    }
  }
  mdf = melt(all_results, id.vars = names(all_results[[1]][[1]]))
  return(mdf)
}


GSEAlm = function(Y, genesets){ # NOTE: method modified from Shahzad Asif 2010
  genesets = lapply(genesets, intersect, names(Y) )
  genesets = Filter(Negate(function(x) length(x) == 0  ), genesets)
  if( length(genesets) == 0 )
    return(NULL)
  gsm = geneset.list2binaryMat(gsl = genesets, Y = Y)
  gsm = gsm[ , ! duplicated(t(gsm)), drop=FALSE]
  gsm = gsm[ , apply(gsm, 2, function(x) any(x == '0') ), drop=FALSE]
  if( ncol(gsm) == 0 )
    return(NULL)
  if( sum(gsm == '0') == 0 )
    return(NULL)
  RESULT =  apply( gsm[ names(Y), , drop=FALSE], 2, function(X)  single_gsealm(Y, X = X))
  RESULT = melt(lapply(RESULT, function(re) data.frame(ES = re$score, pval = re$pval)  ), id.vars = c('ES', 'pval') )
  names(RESULT)[3] = 'feature'
  return(RESULT[, 3:1])
}


geneset.list2binaryMat = function(gsl, Y){
  mgs = melt(gsl)
  names(mgs) = c('gene', 'setname')
  bmgs = dcast(mgs, gene~setname, fill = 0, value.var = 'setname')
  rownames(bmgs) = as.character(bmgs$gene)
  bmgs = bmgs[,-1, drop=FALSE]
  genes2add = setdiff(names(Y), rownames(bmgs))
  bmgs = rbind(bmgs, 
               matrix(0, ncol = ncol(bmgs), nrow = length(genes2add), dimnames = list(genes2add, colnames(bmgs)) ) )
  return(bmgs)
}


single_gsealm = function(Y, X){ 
  names(X) = sapply(strsplit(names(X), split = "\\."), head, 1)
  fit = summary(lm(Y~X))
  rownames(fit$coefficients) = unlist(lapply( rownames(fit$coefficients) , function(x) strsplit(x, ']')[[1]][2] ))
  return(list(score = fit$coefficients[ -1, 'Estimate'], pval = fit$coefficients[ -1, 'Pr(>|t|)' ]))
}



