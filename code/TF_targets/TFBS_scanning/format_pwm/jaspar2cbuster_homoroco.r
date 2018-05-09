

format_matrix = function(ch){
  m = as.matrix(sapply(ch, function(x) as.integer(grep('[0-9]', unlist(strsplit(x, split = '\t')), value = T))) )
  colnames(m) = c('A', 'C', 'G', 'T')
  return(m)
}

hocomoco = read.delim(file = 'data/TF_target_sources/TFBS_scanning/TFBS_models/hocomoco/HOCOMOCOv10_HUMAN_mono_jaspar_format.txt', sep = ',', stringsAsFactors = F, header = F)
entry_idx = grep('^>', hocomoco$V1)
start_idx = entry_idx + 1
end_idx = c((entry_idx - 1)[-1], nrow(hocomoco))
hocomoco_l = mapply(function(entry, st, en){ 
  format_matrix(hocomoco[ st:en, 1])
  }, entry_idx, start_idx, end_idx)
names(hocomoco_l) = hocomoco[entry_idx, 1]

for ( i  in  seq_along(hocomoco_l) ){
  x = rbind(c(names(hocomoco_l[i]), '', '', ''), hocomoco_l[[i]])
  write.table(x, file = paste('data/TF_target_sources/TFBS_scanning/TFBS_models/hocomoco/cbuster_individual/HOCOMOCOv10_HUMAN_mono_', gsub('>', '', names(hocomoco_l)[i]), sep = ''), sep = '\t', col.names = F, row.names = F, quote = F)
}
