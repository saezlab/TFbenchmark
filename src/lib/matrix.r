


get_percentile = function(ob_nes, ra_nes){
  per = ecdf(ra_nes)(ob_nes)
}



effectsize_cohensD  =  function(x, y) { # Effect size
  lx  =  length(x)- 1
  ly  =  length(y)- 1
  md   =  abs(mean(x) - mean(y)) ## mean difference (numerator)
  csd  =  lx * var(x) + ly * var(y)
  csd  =  csd/(lx + ly)
  csd  =  sqrt(csd)                     ## common sd computation
  cd   =  md/csd                        ## cohen's d
}



average_duplicated_genes = function(X){
  genes = setdiff(rownames(X), '')
  mat = matrix(NA, ncol = ncol(X), nrow = length(genes), dimnames = list(genes, colnames(X)) )
  for( g in genes ){
    if( sum(rownames(X) == g) > 1  ){
      mat[g, ] = apply( X[ rownames(X) == g,  ], 2, mean)
    }else{
      mat[g, ] = X[ rownames(X) == g,  ]
    }
  }
  return(mat)
}