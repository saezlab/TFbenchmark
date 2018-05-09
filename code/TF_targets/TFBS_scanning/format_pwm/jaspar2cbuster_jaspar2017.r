

format_matrix = function(ch){
  ch = gsub('\\[', '\\[ ', ch)
  ch = gsub(']', ' ]', ch)
  m = as.matrix(sapply(ch, function(x) as.integer(grep('[0-9]', unlist(strsplit(x, split = ' ')), value = T))) )
  colnames(m) = c('A', 'C', 'G', 'T')
  return(m)
}

jaspar = read.delim(file = 'data/TF_target_sources/TFBS_scanning/TFBS_models/jarpar2017/pfm_vertebrates.pfm', sep = '\t', stringsAsFactors = F, header = F)
entry_idx = which(jaspar$V2 != '')
start_idx = entry_idx + 1
end_idx = c((entry_idx - 1)[-1], nrow(jaspar))
jaspar_l = mapply(function(entry, st, en){ 
  format_matrix(jaspar[ st:en, 1])
  }, entry_idx, start_idx, end_idx)
names(jaspar_l) = paste(jaspar[entry_idx, 1], jaspar[entry_idx, 2], sep = '_')

for ( i  in  seq_along(jaspar_l) ){
  x = rbind(c(names(jaspar_l[i]), '', '', ''), jaspar_l[[i]])
  write.table(x, file = paste('data/TF_target_sources/TFBS_scanning/TFBS_models/jarpar2017/cbuster_individual/pfm_vertebrates.cbuster_', gsub('>', '', names(jaspar_l)[i]), sep = ''), sep = '\t', col.names = F, row.names = F, quote = F)
}
