


retieve_aracne_network = function(regname, E, filter = F){
  message('- Id conversion ...')
  reg = get(regname)
  ensemblgenes_annot = unique(ensemblgenes_annot[, c('hgnc_symbol', 'entrezgene') ])
  reg = reg[ names(reg) %in% ensemblgenes_annot$entrezgene ]
  names(reg) = ensemblgenes_annot$hgnc_symbol[ match(names(reg), ensemblgenes_annot$entrezgene) ]
  reg = reg[order(names(reg))]
  dup = names(reg)[ duplicated(names(reg))]
  reg = reg[ ! names(reg) %in% dup ]
  reg = reg[ names(reg) %in% TFcensus ]
  pb = txtProgressBar(min = 0, max = length(names(reg)), style = 3)
  for( gene in names(reg)){
    setTxtProgressBar(pb, value = which(names(reg)==gene))
    reg[[gene]]$likelihood = reg[[gene]]$likelihood[ names(reg[[gene]]$tfmode) %in% ensemblgenes_annot$entrezgene ]
    reg[[gene]]$tfmode = reg[[gene]]$tfmode[ names(reg[[gene]]$tfmode) %in% ensemblgenes_annot$entrezgene ]
    if( is.null(reg[[gene]]) ){
      reg = reg[ ! names(gene) %in% gene]
    }else{
      names(reg[[gene]]$tfmode) = ensemblgenes_annot$hgnc_symbol[ match(names(reg[[gene]]$tfmode), ensemblgenes_annot$entrezgene) ]      
    }
  }
  close(pb)
  return(reg)
}


load_hugo_entrez_annotation = function(){
  get(load('/Volumes/GoogleDrive/My Drive/databases/ensemblgenes_annot.Rdata'))
}